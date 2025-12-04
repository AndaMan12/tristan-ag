from functools import partial
import numpy as np
from graphet import Data
from graphet.plugins import TristanV2
from tqdm import tqdm
import time
import os
import pickle
import multiprocessing


class RadiationReactionCalculator:
    def __init__(self, out_dir, input_file_name, bounds_x, bounds_y, 
                 interval, omegap0, wA_wp, weight_fac, unit_ch, CC, stride):
        """
        Initialize the radiation reaction calculator
        
        Parameters:
        -----------
        out_dir : str
            Output directory path
        input_file_name : str  
            Configuration file name
        bounds_x, bounds_y : list
            Simulation domain boundaries
        interval : float
            Time interval between frames
        omegap0, wA_wp : float
            Plasma frequency parameters
        weight_fac, unit_ch : float
            Particle weight and unit charge
        CC : float
            Speed of light
        stride : int
            Stride factor for field calculation
        """
        self.out_dir = out_dir
        self.input_file_name = input_file_name
        self.bounds_x = bounds_x
        self.bounds_y = bounds_y
        self.interval = interval
        self.omegap0 = omegap0
        self.wA_wp = wA_wp
        self.weight_fac = weight_fac
        self.unit_ch = unit_ch
        self.CC = CC
        self.stride = stride
        self.metric = np.diag([-1, 1, 1, 1])  # Minkowski metric
        
    def get_test_particle(self, frame, species, particle_id=None, region=None):
        """
        Get test particle data at specified frame
        
        Parameters:
        -----------
        frame : int
            Frame number
        species : int
            Particle species (1 for positrons, 2 for electrons)
        particle_id : int, optional
            Specific particle ID to track
        region : dict, optional
            Region to randomly select from: {'x_range': [min,max], 'y_range': [min,max]}
        
        Returns:
        --------
        dict : Test particle data
        """
        d = Data(TristanV2, steps=[frame], path=self.out_dir, 
                cfg_fname=self.input_file_name, params=True)
        
        particles = d.particles[species].sel(t=frame)
        
        if particle_id is not None:
            # Find specific particle
            mask = particles.idx.values == particle_id
            if not np.any(mask):
                raise ValueError(f"Particle ID {particle_id} not found in frame {frame}")
            idx = np.where(mask)[0][0]
        elif region is not None:
            # Random particle in region
            x_mask = ((particles.x.values >= region['x_range'][0]) & 
                     (particles.x.values <= region['x_range'][1]))
            y_mask = ((particles.y.values >= region['y_range'][0]) & 
                     (particles.y.values <= region['y_range'][1]))
            region_mask = x_mask & y_mask
            
            if not np.any(region_mask):
                raise ValueError("No particles found in specified region")
            
            valid_indices = np.where(region_mask)[0]
            print(particles.x.values[valid_indices])
            print(particles.y.values[valid_indices])
            print(particles.idx.values[valid_indices])
            print(valid_indices)
            idx = valid_indices[1] #7400000000273 #np.random.choice(valid_indices)
        else:
            # Random particle from entire domain
            idx = np.random.randint(0, len(particles.idx.values))
        
        # Extract particle data
        test_particle = {
            'id': particles.idx.values[idx],
            'x': particles.x.values[idx],
            'y': particles.y.values[idx], 
            'z': particles.z.values[idx],
            'u': particles.u.values[idx],
            'v': particles.v.values[idx],
            'w': particles.w.values[idx],
            'species': species,
            'mass': 1 * self.weight_fac,
            'charge': self.unit_ch * self.weight_fac if species == 1 else -self.unit_ch
        }
        
        return test_particle
    
    def delta_s2(self, x4_1, x4_2):
        """Calculate spacetime interval squared using Minkowski metric"""
        x_diff = x4_1 - x4_2
        return np.einsum('ij,jk,ik->i', x_diff, self.metric, x_diff)

    def calculate_field_contribution_from_frame(self, test_particle, observation_frame, past_frame):
        """
        Calculate E and B field contribution from a single past frame
        
        Parameters:
        -----------
        test_particle : dict
            Test particle data at observation frame
        observation_frame : int
            Frame when test particle is observed
        past_frame : int
            Past frame to calculate field contribution from
            
        Returns:
        --------
        tuple : (E_field_contribution, B_field_contribution) as 3D vectors
        """
        
        # Test particle 4-vector position at observation time
        test_pos_3d = np.array([test_particle['x'], test_particle['y'], test_particle['z']])
        x4_test = np.array([self.CC * observation_frame * self.interval, 
                           test_particle['x'], test_particle['y'], test_particle['z']])
        
        E_contribution = np.zeros(3)
        B_contribution = np.zeros(3)
        
        # Load particle data for this past frame and the next one
        try:
#             print("Field calculation...inside try")
            d_curr = Data(TristanV2, steps=[past_frame], path=self.out_dir,
                         cfg_fname=self.input_file_name, params=True)
            d_next = Data(TristanV2, steps=[past_frame + 1], path=self.out_dir,
                         cfg_fname=self.input_file_name, params=True)
#             print("Data loaded!")
        except:
            return E_contribution, B_contribution  # Return zeros if data not available
            print("Loading failed :( ")
        
        # Process both species
        for species in [1, 2]:
            particles_curr = d_curr.particles[species].sel(t=past_frame)
            particles_next = d_next.particles[species].sel(t=past_frame + 1)
            
            # Find common particle IDs
            current_ids = set(particles_curr.idx.values)
            next_ids = set(particles_next.idx.values)
            common_ids = list(current_ids.intersection(next_ids))
            
            if len(common_ids) == 0:
                continue
                
            # # Remove test particle from calculation if in same species
            # if test_particle['species'] == species:
            #     common_ids = [pid for pid in common_ids if pid != test_particle['id']]
            
            if len(common_ids) == 0:
                continue
            
            # Create index mappings
            curr_id_to_idx = {id_val: i for i, id_val in enumerate(particles_curr.idx.values)}
            next_id_to_idx = {id_val: i for i, id_val in enumerate(particles_next.idx.values)}
            
            curr_indices = np.array([curr_id_to_idx[idx] for idx in common_ids])
            next_indices = np.array([next_id_to_idx[idx] for idx in common_ids])
            
            # Construct 4-vector positions for light cone calculation
            curr_locs4 = np.column_stack((self.CC * past_frame * self.interval * np.ones_like(curr_indices),
                                         particles_curr.x.values[curr_indices],
                                         particles_curr.y.values[curr_indices],
                                         particles_curr.z.values[curr_indices]))
            
            next_locs4 = np.column_stack((self.CC * (past_frame + 1) * self.interval * np.ones_like(next_indices),
                                         particles_next.x.values[next_indices],
                                         particles_next.y.values[next_indices],
                                         particles_next.z.values[next_indices]))
            
            # Calculate spacetime intervals to test particle at observation time
            ds2_curr = self.delta_s2(curr_locs4, x4_test)
            ds2_next = self.delta_s2(next_locs4, x4_test)
            
            # Light cone intersection condition: ds2 changes sign between frames
            light_cone_mask = (ds2_curr * ds2_next) <= 0
            
            if not np.any(light_cone_mask):
                continue
            
            # Get indices of particles crossing light cone
            lc_indices = np.where(light_cone_mask)[0]
            lc_curr_indices = curr_indices[lc_indices]
            lc_next_indices = next_indices[lc_indices]
            
            # Get 3D positions of light cone particles
            lc_curr_pos = np.column_stack((particles_curr.x.values[lc_curr_indices],
                                         particles_curr.y.values[lc_curr_indices],
                                         particles_curr.z.values[lc_curr_indices]))
            lc_next_pos = np.column_stack((particles_next.x.values[lc_next_indices],
                                         particles_next.y.values[lc_next_indices],
                                         particles_next.z.values[lc_next_indices]))
            
            # Average positions and calculate distances
            lc_avg_pos = (lc_curr_pos + lc_next_pos) / 2
            lc_displacements = test_pos_3d - lc_avg_pos
            lc_distances = np.linalg.norm(lc_displacements, axis=1)
            
            # Calculate unit vectors from particles to test particle
            lc_unit_vectors = lc_displacements / lc_distances[:, np.newaxis]
            
            # Get 4-velocities of light cone particles
            lc_curr_4vel = np.column_stack((particles_curr.u.values[lc_curr_indices],
                                          particles_curr.v.values[lc_curr_indices],
                                          particles_curr.w.values[lc_curr_indices]))
            lc_next_4vel = np.column_stack((particles_next.u.values[lc_next_indices],
                                          particles_next.v.values[lc_next_indices],
                                          particles_next.w.values[lc_next_indices]))
            
            # Convert to beta (v/c)
            lc_curr_beta = lc_curr_4vel / np.sqrt(1 + np.linalg.norm(lc_curr_4vel, axis=1)**2)[:, np.newaxis]
            lc_next_beta = lc_next_4vel / np.sqrt(1 + np.linalg.norm(lc_next_4vel, axis=1)**2)[:, np.newaxis]
            lc_avg_beta = (lc_curr_beta + lc_next_beta) / 2
            lc_avg_gamma = 1 / np.sqrt(1 - np.linalg.norm(lc_avg_beta, axis = 1)**2)
            lc_beta_dot = (lc_next_beta - lc_curr_beta) / self.interval
            
            # Calculate electric field contribution (Liénard-Wiechert fields)
            charge_sign = 1 if species == 1 else -1
            E_acc_num = (charge_sign * self.weight_fac * self.unit_ch / self.CC * 
                    np.cross(lc_unit_vectors, np.cross(lc_unit_vectors - lc_avg_beta, lc_beta_dot)))
            E_acc_denom = ((1 - np.einsum('ij,ij->i', lc_unit_vectors, lc_avg_beta))**3 * lc_distances)
            
            E_vel_num = charge_sign * self.weight_fac * self.unit_ch * (lc_unit_vectors - lc_avg_beta)
            E_vel_denom = E_acc_denom * lc_distances * lc_avg_gamma**2 

            # Avoid division by zero
            valid_mask = (E_acc_denom != 0) & (E_vel_denom != 0)
            if np.any(valid_mask):
                E_fields = np.zeros_like(E_acc_num)
                E_fields[valid_mask] = E_acc_num[valid_mask] / E_acc_denom[valid_mask, np.newaxis] + E_vel_num[valid_mask] / E_vel_denom[valid_mask, np.newaxis]
                
                # Calculate magnetic field
                B_fields = np.cross(lc_unit_vectors, E_fields)
                
                # Sum contributions from this frame
                E_contribution += np.sum(E_fields, axis=0)
                B_contribution += np.sum(B_fields, axis=0)
        
        return E_contribution * self.stride, B_contribution * self.stride
    
    def calculate_forces(self, test_particle, E_field, B_field):
        """
        Calculate radiation reaction force using Faraday tensor
        
        Parameters:
        -----------
        test_particle : dict
            Test particle data
        E_field, B_field : array
            Electric and magnetic field vectors
            
        Returns:
        --------
        array : 4-force from radiation reaction
        """
        # Test particle 4-velocity
        u_4vel = np.array([test_particle['u'], test_particle['v'], test_particle['w']])
        gamma = np.sqrt(1 + np.linalg.norm(u_4vel)**2)
        u_mu = np.array([gamma, test_particle['u'], test_particle['v'], test_particle['w']]) # u^\mu
        
        # Construct Faraday tensor F_μν
        F_tensor = np.zeros((4, 4))
        F_tensor[0, 1:4] = E_field
        F_tensor[1:4, 0] = -E_field
        F_tensor[1, 2] = B_field[2]
        F_tensor[2, 1] = -B_field[2]
        F_tensor[1, 3] = -B_field[1]
        F_tensor[3, 1] = B_field[1]
        F_tensor[2, 3] = B_field[0]
        F_tensor[3, 2] = -B_field[0]

        F_contra_tensor = np.einsum('ij,jk,kl->il', self.metric, F_tensor, self.metric)
        
        # Calculate F_μν u^ν
        F_dot_u = np.einsum('ij,j->i', F_tensor, u_mu)
        
        # Calculate (F_αβ u^β)^2
        F_dot_u_squared = np.einsum('i,ij,j', F_dot_u, self.metric, F_dot_u)

        # Calculate F^αβ F_βν u^ν
        F_dot_F_dot_u = np.einsum('ij,j->i', F_contra_tensor, F_dot_u)
        
        # Radiation reaction 4-force: (2e^4/3m^2)[F^αβ F_βν u^ν + (F_μν u^ν)² u^α]        
        mass_squared = test_particle['mass']**2  # Normalized units
        prefactor = 2 * (test_particle['charge'])**4 / (3 * mass_squared)
        
        f_rad = prefactor * (F_dot_F_dot_u + F_dot_u_squared * u_mu)

        f_lorentz = test_particle['charge'] * np.einsum('ij,j->i', self.metric, F_dot_u)
        
        return f_rad, f_lorentz
    
    
    
    def distance_analysis(self, test_particle, frame, max_frames_back, normalize=True):
        """
        Efficiently analyze how radiation reaction force varies with retarded time/distance
        Uses incremental field calculation to avoid repeated computations
        
        Parameters:
        -----------
        test_particle : dict
            Test particle data
        frame : int
            Observation frame number
        max_frames_back : int
            Maximum number of past frames to include
        normalize : bool
            Whether to normalize by Lorentz force magnitude
            
        Returns:
        --------
        dict : Results containing forces vs retarded distance
        """
        frame_ranges = range(1, max_frames_back + 1)
        max_distances = np.array(frame_ranges) * self.interval * self.CC  # Convert frames to distance
        
        results = {
            'test_particle': test_particle,
            'frames_back': list(frame_ranges),
            'max_distances': max_distances,
            'f_rad': [],
            'f_lorentz': []            
        }
        
        # Initialize cumulative fields
        E_total = np.zeros(3)
        B_total = np.zeros(3)
        
        for frames_back in tqdm(frame_ranges, desc="Incremental retarded time analysis"):
            # Calculate the NEW contribution from the most distant frame
            past_frame = frame - frames_back
            
            if past_frame >= 0:  # Only if frame existsfrom multiprocessing.dummy import Pool as ThreadPool
                E_new, B_new = self.calculate_field_contribution_from_frame(
                    test_particle, frame, past_frame)
                
                # Add to cumulative total
                E_total += E_new
                B_total += B_new
            
            # Calculate forces using current cumulative fields
            f_rad, f_lorentz = self.calculate_forces(test_particle, E_total, B_total)            
            results['f_rad'] = f_rad
            results['f_lorentz'] = f_lorentz            
            
#             # Normalized ratio
#             if normalize and lorentz_mag > 0:
#                 results['normalized_ratios'].append(rad_mag / lorentz_mag)
#             else:
#                 results['normalized_ratios'].append(rad_mag)        
        
        return results
    
#     def distance_analysis_chunked(self, test_particle, frame, max_frames_back, num_chunks=5):
#         """
#         Coarse-grained parallel version of distance_analysis.
#         Each process handles a subset of past frames serially.

#         Parameters
#         ----------
#         test_particle : dict
#             Test particle data.
#         frame : int
#             Observation frame number.
#         max_frames_back : int
#             Maximum number of past frames to include.
#         num_chunks : int
#             Number of coarse parallel chunks (default 5).

#         Returns
#         -------
#         dict
#             Results with same structure as original distance_analysis2().
#         """

#         frame_ranges = list(range(1, max_frames_back + 1))
#         max_distances = np.array(frame_ranges) * self.interval * self.CC

#         results = {
#             "test_particle": test_particle,
#             "frames_back": frame_ranges,
#             "max_distances": max_distances,
#             "f_rad": [],
#             "f_lorentz": [],
#         }

#         past_frames = [frame - i for i in frame_ranges if frame - i >= 0]
#         if not past_frames:
#             return results

#         # Split past frames into chunks
#         chunks = np.array_split(past_frames, num_chunks)

#         calc_params = {
#             "out_dir": self.out_dir,
#             "input_file_name": self.input_file_name,
#             "bounds_x": self.bounds_x,
#             "bounds_y": self.bounds_y,
#             "interval": self.interval,
#             "omegap0": self.omegap0,
#             "wA_wp": self.wA_wp,
#             "weight_fac": self.weight_fac,
#             "unit_ch": self.unit_ch,
#             "CC": self.CC,
#             "stride": self.stride,
#         }

#         # Define worker function at top level of this module
#         def _coarse_worker(calc_params, test_particle, obs_frame, subset, output_file):
#             """
#             Serially process one subset of past frames and write partial results.
#             """
#             calc = RadiationReactionCalculator(**calc_params)
#             E_total = np.zeros(3)
#             B_total = np.zeros(3)
#             partial = []

#             for pf in subset:
#                 E_new, B_new = calc.calculate_field_contribution_from_frame(test_particle, obs_frame, pf)
#                 E_total += E_new
#                 B_total += B_new
#                 f_rad, f_lorentz = calc.calculate_forces(test_particle, E_total, B_total)
#                 partial.append((pf, f_rad, f_lorentz))

#             with open(output_file, "wb") as f:
#                 pickle.dump(partial, f)

#         # --- Launch parallel processes ---
#         jobs = []
#         for i, subset in enumerate(chunks):
#             if len(subset) == 0:
#                 continue
#             output_file = f"chunk_{i}.pkl"
#             args = (calc_params, test_particle, frame, subset.tolist(), output_file)
#             p = multiprocessing.Process(target=_coarse_worker, args=args)
#             p.start()
#             jobs.append(p)

#         for p in tqdm(jobs, desc="Processing chunks"):
#             p.join()

#         # --- Merge partial outputs ---
#         combined = []
#         for i, subset in enumerate(chunks):
#             output_file = f"chunk_{i}.pkl"
#             if os.path.exists(output_file):
#                 with open(output_file, "rb") as f:
#                     combined.extend(pickle.load(f))
#                 os.remove(output_file)

#         combined.sort(key=lambda r: r[0])  # sort by frame number

#         # Fill in same structure as before
#         for _, f_rad, f_lorentz in combined:
#             results["f_rad"].append(f_rad)
#             results["f_lorentz"].append(f_lorentz)

#         return results


# Example usage function
def run_radiation_reaction_analysis(out_dir, input_file_name, bounds_x, bounds_y,
                                   interval, omegap0, wA_wp, weight_fac, unit_ch, 
                                   CC, stride, frame, species=2, max_frames_back=20):
    """
    Run complete radiation reaction analysis with proper retarded time
    
    Parameters:
    -----------
    All simulation parameters as in original code
    frame : int
        Observation frame to analyze
    species : int  
        Particle species (1=positron, 2=electron)
    max_frames_back : int
        Maximum number of past frames to include (controls max retarded distance)
    """
    
    # Initialize calculator
    calc = RadiationReactionCalculator(
        out_dir, input_file_name, bounds_x, bounds_y,
        interval, omegap0, wA_wp, weight_fac, unit_ch, CC, stride
    )
    
    # Get random test particle
    test_particle = calc.get_test_particle(frame, species, region = {'x_range':bounds_x, 'y_range':bounds_y})
    print(f"Selected test particle: ID={test_particle['id']}, "
          f"pos=({test_particle['x']:.2f}, {test_particle['y']:.2f}, {test_particle['z']:.2f})")
    
    # Run retarded time analysis
    results = calc.distance_analysis_chunked(test_particle, frame, max_frames_back, num_chunks = 15)
    with open('rad_rxn_output.pkl', 'wb') as f:
        pickle.dump(results, f)
        f.close()

    print(f"Analysis complete. Max retarded distance: {results['max_distances'][-1]:.2f}")
    
    return results, test_particle




if __name__ == "__main__":
    
    out_dir = "/scratch/10446/anindya_12/tristan_mp_v2/vault/output_psi0.5_mul1_3_mul2_0.1_TT_1.9e-4_rad_drag/"
    input_file_name = out_dir + "temp_input_psi0.5_mul13_mul20.1_TT1.9e-4.in"
    
    # Parse input parameters to extract values
    import tristanVis.isolde as isolde
    input_params = isolde.parseInput(input_file_name)
    
    interval = input_params["output"]["interval"]
    lst_time = input_params["time"]["last"]
    grid_x = int(input_params["grid"]["mx0"])
    grid_y = int(input_params["grid"]["my0"])
    CC = input_params["algorithm"]["c"]
    COMP = input_params["plasma"]["c_omp"]
    ppc0 = input_params["plasma"]["ppc0"]
    stride = input_params["output"]["stride"]
    wA_wp = input_params["problem"]["wA_wp"]
    mode = input_params["problem"]["mode"]
    ramp_width = input_params["problem"]["ramp_width"]
    
    # Derived parameters
    unit_ch = CC**2 / (ppc0 * COMP**2)
    weight_fac = 2.0473518E-03  # From your notebook
    ds = wA_wp * 2 * np.pi / mode
    omegap0 = CC / ds
    init_x_boundary = int(5 * np.pi / mode) + 1
    fin_x_boundary = init_x_boundary + ramp_width * (2 * np.pi / mode)
    
    # Set analysis parameters (modify these as needed)
    bounds_x = [10220, 10250]  # e.g., [init_x_boundary, fin_x_boundary]
    bounds_y = [1040, 1070]  # e.g., [0, grid_y]
    frame = 2400  # e.g., 100 or int(lst_time // interval) // 2
    species = 2  # 1 for positrons, 2 for electrons
    max_frames_back = 150  # Maximum retarded time frames
    
#     print(f"Grid size: {grid_x} x {grid_y}")
#     print(f"Boundaries calculated: x=[{init_x_boundary}, {fin_x_boundary}]")
   
    
    results, test_particle = run_radiation_reaction_analysis(
        out_dir=out_dir,
        input_file_name=input_file_name, 
        bounds_x=bounds_x,
        bounds_y=bounds_y,
        interval=interval,
        omegap0=omegap0,
        wA_wp=wA_wp,
        weight_fac=weight_fac,
        unit_ch=unit_ch,
        CC=CC,
        stride=stride,
        frame=frame,
        species=species,
        max_frames_back=max_frames_back
    )
