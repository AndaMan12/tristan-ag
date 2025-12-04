# radiation_utils.py
from functools import partial
import numpy as np
from graphet import Data
from graphet.plugins import TristanV2
from tqdm import tqdm
import time
import os
import pickle
import multiprocessing
import tristanVis.isolde as isolde

class RadiationReactionCalculator:
    def __init__(self, out_dir, input_file_name, bounds_x, bounds_y, 
                 interval, omegap0, wA_wp, weight_fac, unit_ch, CC, stride):
        """
        RadiationReactionCalculator: cleaned and frame-centric friendly.

        Parameters are the same as your original script. Key fields:
          - out_dir, input_file_name: Tristan output config and path
          - bounds_x, bounds_y: region for test-particle selection
          - interval: time between frames (simulation output time)
          - weight_fac, unit_ch: weight and charge scaling used in your pipeline
          - CC: code units for speed of light
          - stride: multiplier applied to returned fields (kept to preserve original behavior)
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

        # Small softening length to avoid singular 1/r when r -> 0
        # You may adjust this (e.g. to grid cell size)
        self.softening = 1e-8

    # ------------------------------------------------------------------
    # Utility: small wrapper to fetch prtl/flds/spec files via isolde
    # ------------------------------------------------------------------
    @staticmethod
    def fetch_var_at_step(out_dir, var, step, rad=False):
        """
        Fetch variable file (prtl,flds,spec) from a directory and step index.
        Returns the parsed object from isolde.
        """
        filename = os.path.join(out_dir.rstrip('/'), var, f"{var}.tot.{step:05d}")
        if var == "prtl":
            return isolde.getParticles(filename)
        elif var == "flds":
            return isolde.getFields(filename)
        elif var == "spec":
            return isolde.getSpectra(filename, radiation=rad)
        else:
            raise ValueError(f"Unsupported var '{var}' in fetch_var_at_step")

    # ------------------------------------------------------------------
    # Frame loading helpers -- driver will prefer to call these and cache
    # ------------------------------------------------------------------
    def load_frame_pair(self, frame):
        """
        Load particle Data for 'frame' and 'frame+1' using graphet.Data(TristanV2).
        Returns (d_curr, d_next) or (None, None) if load fails.
        """
        try:
            d_curr = Data(TristanV2, steps=[frame], path=self.out_dir, cfg_fname=self.input_file_name, params=True)
            d_next = Data(TristanV2, steps=[frame + 1], path=self.out_dir, cfg_fname=self.input_file_name, params=True)
            return d_curr, d_next
        except Exception as e:
            print(f"[load_frame_pair] Failed to load frames {frame}/{frame+1}: {e}")
            return None, None

    def load_particles_from_frame(self, frame):
        """
        Convenience wrapper: returns the Data object for particles at given frame.
        (Equivalent to Data(...).particles but returns the Data object so the caller
         can access both species easily.)
        """
        try:
            d = Data(TristanV2, steps=[frame], path=self.out_dir, cfg_fname=self.input_file_name, params=True)
            return d
        except Exception as e:
            print(f"[load_particles_from_frame] Failed to load frame {frame}: {e}")
            return None

    # ------------------------------------------------------------------
    # Particle selection (fixes and optional density_dir)
    # ------------------------------------------------------------------
    def get_test_particle(self, frame, region=None, density_dir=None, percentile=25):
        """
        Select test particles in `region` at given `frame` using local density > percentile.
        If density_dir is provided it will be used to read flds files (useful when density
        diagnostics come from a different simulation run).
        
        Returns:
            selected_particles : list of dicts (same keys as in original code)
        """
        density_dir = density_dir or self.out_dir

        # Load Tristan particle Data for the observation frame
        try:
            d = Data(TristanV2, steps=[frame], path=self.out_dir, cfg_fname=self.input_file_name, params=True)
        except Exception as e:
            raise RuntimeError(f"Failed to load particle data for frame {frame}: {e}")

        selected_particles = []

        # Process species 1 and 2
        for species in (1, 2):
            particles = d.particles[species].sel(t=frame)

            # Read density field from density_dir/flds (user may provide different directory)
            dens_step = frame // 5 if frame >= 5 else 0  # your original code used frame//5
            try:
                flds = self.fetch_var_at_step(density_dir, "flds", dens_step)
            except Exception as e:
                raise RuntimeError(f"Failed to load fields from {density_dir} for step {dens_step}: {e}")

            dens_key = f"dens{species}"
            if dens_key not in flds:
                raise KeyError(f"Density key '{dens_key}' not found in fields file at step {dens_step}")

            dens_data = flds[dens_key]
            # dens_data may be 2D or 3D (z,y,x). Averaging over first axis (z) as in your original code
            dens_avg = np.mean(dens_data, axis=0)

            # Build mask for region (particle positions assumed to be in same coordinate system as dens grid indexing)
            if region is not None:
                x_mask = (particles.x.values >= region['x_range'][0]) & (particles.x.values <= region['x_range'][1])
                y_mask = (particles.y.values >= region['y_range'][0]) & (particles.y.values <= region['y_range'][1])
                region_mask = x_mask & y_mask
            else:
                region_mask = np.ones_like(particles.x.values, dtype=bool)

            # Determine threshold (percentile) on the sub-region density if region provided,
            # otherwise use global dens_avg percentile
            if region is not None:
                # region indices are being used as integer grid indices in original code:
                x_min, x_max = region['x_range']
                y_min, y_max = region['y_range']
                # guard against out-of-bounds slices
                Ny, Nx = dens_avg.shape
                x_min_i = max(0, int(x_min)); x_max_i = min(Nx, int(x_max))
                y_min_i = max(0, int(y_min)); y_max_i = min(Ny, int(y_max))
                if x_max_i <= x_min_i or y_max_i <= y_min_i:
                    # fallback to full domain
                    region_dens_vals = dens_avg.flatten()
                else:
                    region_dens_vals = dens_avg[y_min_i:y_max_i, x_min_i:x_max_i].flatten()
            else:
                region_dens_vals = dens_avg.flatten()

            if region_dens_vals.size == 0:
                dens_threshold = 0.0
            else:
                dens_threshold = np.percentile(region_dens_vals, percentile)

            # Map particle positions to grid indices (nearest integer). This assumes particle
            # coordinates are referenced in the same units as grid indices (as in your original).
            Ny, Nx = dens_avg.shape
            x_idx = np.clip(particles.x.values.astype(int), 0, Nx - 1)
            y_idx = np.clip(particles.y.values.astype(int), 0, Ny - 1)
            particle_local_dens = dens_avg[y_idx, x_idx]

            # Apply selection: local density > threshold and inside region mask
            high_dens_mask = (particle_local_dens > dens_threshold) & region_mask

            for idx in np.where(high_dens_mask)[0]:
                selected_particles.append({
                    'id': particles.idx.values[idx],
                    'x': particles.x.values[idx],
                    'y': particles.y.values[idx],
                    'z': particles.z.values[idx],
                    'u': particles.u.values[idx],
                    'v': particles.v.values[idx],
                    'w': particles.w.values[idx],
                    'species': species,
                    'mass': 1 * self.weight_fac,
                    'charge': (self.unit_ch * self.weight_fac if species == 1 else -self.unit_ch),
                    'density': float(particle_local_dens[idx])
                })

        print(f"[get_test_particle] Selected {len(selected_particles)} particles at frame {frame}")
        return selected_particles

    # ------------------------------------------------------------------
    # Minkowski interval helper
    # ------------------------------------------------------------------
    def delta_s2(self, x4_array, x4_single):
        """
        x4_array: (N,4) array of events (ct, x, y, z)
        x4_single: (4,) single event (ct, x, y, z)
        returns ds^2 for each row: (x - x0)^T * metric * (x - x0)
        """
        x_diff = x4_array - x4_single[np.newaxis, :]
        # metric diag(-1,1,1,1): ds^2 = -dt^2 + dx^2 + dy^2 + dz^2
        return -x_diff[:, 0]**2 + np.sum(x_diff[:, 1:]**2, axis=1)

    # ------------------------------------------------------------------
    # Main field contribution (now accepts preloaded Data objects)
    # ------------------------------------------------------------------
    def calculate_field_contribution_from_frame(self, test_particle, observation_frame, past_frame, d_curr=None, d_next=None):
        """
        Calculate E and B contribution from source particles whose worldlines
        cross the light cone between (past_frame) and (past_frame+1),
        evaluated at test_particle's spacetime point (observation_frame).

        Parameters:
          - test_particle: dict with particle fields at observation (keys: x,y,z,u,v,w,species,id,...)
          - observation_frame: int (frame index where test particle is observed)
          - past_frame: int (frame index to test emission from)
          - d_curr, d_next: optional Data objects preloaded for past_frame and past_frame+1
                             (if not provided, they will be loaded here)
        Returns:
          - (E_contribution (3,), B_contribution (3,))
        """
        # Prepare outputs
        E_contribution = np.zeros(3, dtype=float)
        B_contribution = np.zeros(3, dtype=float)

        # Test particle 3-position and 4-position at observation time
        test_pos_3d = np.array([test_particle['x'], test_particle['y'], test_particle['z']], dtype=float)
        t_obs = observation_frame * self.interval
        x4_test = np.array([self.CC * t_obs, test_pos_3d[0], test_pos_3d[1], test_pos_3d[2]], dtype=float)

        # Load frames if not provided
        if d_curr is None or d_next is None:
            try:
                d_curr = Data(TristanV2, steps=[past_frame], path=self.out_dir, cfg_fname=self.input_file_name, params=True)
                d_next = Data(TristanV2, steps=[past_frame + 1], path=self.out_dir, cfg_fname=self.input_file_name, params=True)
            except Exception as e:
                # If loading fails, return zero contributions
                print(f"[calculate_field_contribution_from_frame] Failed to load frames {past_frame}/{past_frame+1}: {e}")
                return E_contribution, B_contribution

        # Loop species
        for species in (1, 2):
            print("Doing species {}...".format(species))
            try:
                particles_curr = d_curr.particles[species].sel(t=past_frame)
                particles_next = d_next.particles[species].sel(t=past_frame + 1)
            except Exception:
                # if a species is missing in a snapshot, skip gracefully
                continue

            # Find common particle IDs present in both snapshots
            ids_curr = particles_curr.idx.values
            ids_next = particles_next.idx.values
            common_ids = np.intersect1d(ids_curr, ids_next, assume_unique=False)
            if common_ids.size == 0:
                continue

            # Create mapping from id->index for quick selection
            curr_id_to_idx = {val: i for i, val in enumerate(ids_curr)}
            next_id_to_idx = {val: i for i, val in enumerate(ids_next)}

            # Build index arrays for common ids
            curr_indices = np.array([curr_id_to_idx[val] for val in common_ids], dtype=int)
            next_indices = np.array([next_id_to_idx[val] for val in common_ids], dtype=int)

            # 4-positions at past_frame and past_frame+1
            t_curr = self.CC * past_frame * self.interval
            t_next = self.CC * (past_frame + 1) * self.interval
            curr_locs4 = np.column_stack((
                np.full(curr_indices.size, t_curr, dtype=float),
                particles_curr.x.values[curr_indices],
                particles_curr.y.values[curr_indices],
                particles_curr.z.values[curr_indices]
            ))
            next_locs4 = np.column_stack((
                np.full(next_indices.size, t_next, dtype=float),
                particles_next.x.values[next_indices],
                particles_next.y.values[next_indices],
                particles_next.z.values[next_indices]
            ))

            # Compute ds^2 to the observation event for all common particles at both snapshots
            ds2_curr = self.delta_s2(curr_locs4, x4_test)
            ds2_next = self.delta_s2(next_locs4, x4_test)
            
            print("ds2 s are calculated...")
            
            # Light-cone crossing mask: sign change or exactly zero
            light_cone_mask = (ds2_curr * ds2_next) <= 0
            if not np.any(light_cone_mask):
                continue

            # Indices of common_ids that cross the light cone
            lc_idx = np.where(light_cone_mask)[0]
            lc_curr_indices = curr_indices[lc_idx]
            lc_next_indices = next_indices[lc_idx]

            # Spatial positions at the two snapshots (N_lc x 3)
            lc_curr_pos = np.column_stack((
                particles_curr.x.values[lc_curr_indices],
                particles_curr.y.values[lc_curr_indices],
                particles_curr.z.values[lc_curr_indices]
            ))
            lc_next_pos = np.column_stack((
                particles_next.x.values[lc_next_indices],
                particles_next.y.values[lc_next_indices],
                particles_next.z.values[lc_next_indices]
            ))
            
            print("Positions of source particles sorted...")
            
            # Average source positions (linear interp approx)
            lc_avg_pos = 0.5 * (lc_curr_pos + lc_next_pos)
            # Displacement vectors from source -> observer
            lc_displacements = test_pos_3d[None, :] - lc_avg_pos  # shape (N_lc,3)

            # Distances and safe division
            lc_distances = np.linalg.norm(lc_displacements, axis=1)
            # Avoid division by zero
            very_small = lc_distances < self.softening
            if np.any(very_small):
                lc_distances[very_small] = self.softening
                lc_displacements[very_small] = lc_displacements[very_small] / (np.linalg.norm(lc_displacements[very_small], axis=1)[:, None] + 1e-30) * self.softening

            lc_unit_vectors = lc_displacements / lc_distances[:, None]

            # Get particle velocities (spatial components) at the two snapshots
            lc_curr_4vel = np.column_stack((
                particles_curr.u.values[lc_curr_indices],
                particles_curr.v.values[lc_curr_indices],
                particles_curr.w.values[lc_curr_indices]
            ))
            lc_next_4vel = np.column_stack((
                particles_next.u.values[lc_next_indices],
                particles_next.v.values[lc_next_indices],
                particles_next.w.values[lc_next_indices]
            ))

            # Convert to beta (v/c) -- this assumes u,v,w are spatial components of 4-velocity
            # Confirm convention in Tristan: if u/v/w are momentum-like, adjust accordingly.
            # Here we follow your original approach:
            def to_beta(spatial):
                norms = np.linalg.norm(spatial, axis=1)
                denom = np.sqrt(1.0 + norms**2)
                return spatial / denom[:, None]

            lc_curr_beta = to_beta(lc_curr_4vel)
            lc_next_beta = to_beta(lc_next_4vel)
            lc_avg_beta = 0.5 * (lc_curr_beta + lc_next_beta)
            # approximate time derivative of beta
            lc_beta_dot = (lc_next_beta - lc_curr_beta) / self.interval

            # gamma for averaged beta (used in E_vel denom)
            beta2 = np.sum(lc_avg_beta**2, axis=1)
            # guard: if beta2 >= 1 we clamp to slightly less than 1 to avoid nan gamma
            beta2 = np.minimum(beta2, 1.0 - 1e-12)
            lc_avg_gamma = 1.0 / np.sqrt(1.0 - beta2)
            
            print("Velocities of source sorted...")
            
            # Liénard–Wiechert contributions (vectorized)
            charge_sign = np.where(species == 1, 1.0, -1.0)
            qscale = charge_sign * self.weight_fac * self.unit_ch

            # Numerators and denominators (vectorized shapes)
            # E_acc term
            cross_inner = np.cross(lc_unit_vectors, np.cross(lc_unit_vectors - lc_avg_beta, lc_beta_dot))
            dot_n_beta = np.einsum('ij,ij->i', lc_unit_vectors, lc_avg_beta)
            # denom for acceleration term (note: shape (N_lc,))
            E_acc_denom = (1.0 - dot_n_beta)**3 * lc_distances
            # vel/near field term denom
            E_vel_denom = E_acc_denom * lc_distances * (lc_avg_gamma**2)

            # Avoid zero denominators
            valid_acc = np.abs(E_acc_denom) > 0
            valid_vel = np.abs(E_vel_denom) > 0
            if not np.any(valid_acc | valid_vel):
                continue

            # Prepare arrays
            E_fields = np.zeros((lc_unit_vectors.shape[0], 3), dtype=float)

            # acceleration term
            idx_acc = np.where(valid_acc)[0]
            if idx_acc.size > 0:
                E_fields[idx_acc] += (qscale / self.CC) * (cross_inner[idx_acc] / E_acc_denom[idx_acc][:, None])

            # velocity term
            idx_vel = np.where(valid_vel)[0]
            if idx_vel.size > 0:
                E_vel_num = qscale * (lc_unit_vectors[idx_vel] - lc_avg_beta[idx_vel])
                E_fields[idx_vel] += (E_vel_num / E_vel_denom[idx_vel][:, None])

            # Magnetic field from cross(n, E)
            B_fields = np.cross(lc_unit_vectors, E_fields)

            # Sum contributions from this species
            E_contribution += np.sum(E_fields, axis=0)
            B_contribution += np.sum(B_fields, axis=0)
            
            print("Fields computed ...")

        # apply stride scaling to match your previous code's behavior
        return E_contribution * self.stride, B_contribution * self.stride

    # ------------------------------------------------------------------
    # Force computation (unchanged physics, with small safety guards)
    # ------------------------------------------------------------------
    def calculate_forces(self, test_particle, E_field, B_field):
        """
        Calculate radiation reaction 4-force and Lorentz force.

        Returns:
            f_rad (4-vector), f_lorentz (4-vector)
        Note: This uses the same algebra as your original script. Double-check
        that test_particle u,v,w convention matches expected 4-velocity components.
        """
        # spatial u vector (these are the spatial components of 4-velocity in your code)
        u_sp = np.array([test_particle['u'], test_particle['v'], test_particle['w']], dtype=float)
        gamma = np.sqrt(1.0 + np.sum(u_sp**2))
        u_mu = np.array([gamma, u_sp[0], u_sp[1], u_sp[2]], dtype=float)

        # Build Faraday tensor (covariant components)
        F = np.zeros((4, 4), dtype=float)
        F[0, 1:4] = E_field
        F[1:4, 0] = -E_field
        # B-field entries (Cartesian)
        F[1, 2] = B_field[2]
        F[2, 1] = -B_field[2]
        F[1, 3] = -B_field[1]
        F[3, 1] = B_field[1]
        F[2, 3] = B_field[0]
        F[3, 2] = -B_field[0]

        # Raise indices using metric diag(-1,1,1,1)
        F_contra = np.einsum('ij,jk,kl->il', self.metric, F, self.metric)

        # F_{μν} u^ν
        F_dot_u = np.einsum('ij,j->i', F, u_mu)

        # (F_{αβ} u^β)^2 using metric
        F_dot_u_squared = np.einsum('i,ij,j', F_dot_u, self.metric, F_dot_u)

        # F^αβ F_βν u^ν
        F_dot_F_dot_u = np.einsum('ij,j->i', F_contra, F_dot_u)

        mass_squared = (test_particle['mass'])**2
        # Prevent divide by zero mass
        if mass_squared == 0:
            mass_squared = 1e-30

        prefactor = 2.0 * (test_particle['charge'])**4 / (3.0 * mass_squared)
        f_rad = prefactor * (F_dot_F_dot_u + F_dot_u_squared * u_mu)

        # Lorentz 4-force: q * g^{μα} F_{αβ} u^β  (we follow your original expression)
        f_lorentz = test_particle['charge'] * np.einsum('ij,j->i', self.metric, F_dot_u)

        return f_rad, f_lorentz






# from functools import partial
# import numpy as np
# from graphet import Data
# from graphet.plugins import TristanV2
# from tqdm import tqdm
# import time
# import os
# import pickle
# import multiprocessing
# import tristanVis.isolde as isolde

# class RadiationReactionCalculator:
#     def __init__(self, out_dir, input_file_name, bounds_x, bounds_y, 
#                  interval, omegap0, wA_wp, weight_fac, unit_ch, CC, stride):
#         """
#         Initialize the radiation reaction calculator
        
#         Parameters:
#         -----------
#         out_dir : str
#             Output directory path
#         input_file_name : str  
#             Configuration file name
#         bounds_x, bounds_y : list
#             Simulation domain boundaries
#         interval : float
#             Time interval between frames
#         omegap0, wA_wp : float
#             Plasma frequency parameters
#         weight_fac, unit_ch : float
#             Particle weight and unit charge
#         CC : float
#             Speed of light
#         stride : int
#             Stride factor for field calculation
#         """
#         self.out_dir = out_dir
#         self.input_file_name = input_file_name
#         self.bounds_x = bounds_x
#         self.bounds_y = bounds_y
#         self.interval = interval
#         self.omegap0 = omegap0
#         self.wA_wp = wA_wp
#         self.weight_fac = weight_fac
#         self.unit_ch = unit_ch
#         self.CC = CC
#         self.stride = stride
#         self.metric = np.diag([-1, 1, 1, 1])  # Minkowski metric
        
# #     def get_test_particle(self, frame, species, particle_id=None, region=None):
# #         """
# #         Get test particle data at specified frame
        
# #         Parameters:
# #         -----------
# #         frame : int
# #             Frame number
# #         species : int
# #             Particle species (1 for positrons, 2 for electrons)
# #         particle_id : int, optional
# #             Specific particle ID to track
# #         region : dict, optional
# #             Region to randomly select from: {'x_range': [min,max], 'y_range': [min,max]}
        
# #         Returns:
# #         --------
# #         dict : Test particle data
# #         """
# #         d = Data(TristanV2, steps=[frame], path=self.out_dir, 
# #                 cfg_fname=self.input_file_name, params=True)
        
# #         particles = d.particles[species].sel(t=frame)
        
# #         if particle_id is not None:
# #             # Find specific particle
# #             mask = particles.idx.values == particle_id
# #             if not np.any(mask):
# #                 raise ValueError(f"Particle ID {particle_id} not found in frame {frame}")
# #             idx = np.where(mask)[0][0]
# #         elif region is not None:
# #             # Random particle in region
# #             x_mask = ((particles.x.values >= region['x_range'][0]) & 
# #                      (particles.x.values <= region['x_range'][1]))
# #             y_mask = ((particles.y.values >= region['y_range'][0]) & 
# #                      (particles.y.values <= region['y_range'][1]))
# #             region_mask = x_mask & y_mask
            
# #             if not np.any(region_mask):
# #                 raise ValueError("No particles found in specified region")
            
# #             valid_indices = np.where(region_mask)[0]
# # #             print(particles.x.values[valid_indices])
# # #             print(particles.y.values[valid_indices])
# # #             print(particles.idx.values[valid_indices])
# # #             print(valid_indices)
# #             idx = valid_indices[1] #7400000000273 #np.random.choice(valid_indices)
# #         else:
# #             # Random particle from entire domain
# #             idx = np.random.randint(0, len(particles.idx.values))
        
# #         # Extract particle data
# #         test_particle = {
# #             'id': particles.idx.values[idx],
# #             'x': particles.x.values[idx],
# #             'y': particles.y.values[idx], 
# #             'z': particles.z.values[idx],
# #             'u': particles.u.values[idx],
# #             'v': particles.v.values[idx],
# #             'w': particles.w.values[idx],
# #             'species': species,
# #             'mass': 1 * self.weight_fac,
# #             'charge': self.unit_ch * self.weight_fac if species == 1 else -self.unit_ch
# #         }
        
# #         return test_particle
    

#     def get_test_particle(self, frame, region=None):
#         """
#         Get test particle data for all particles in regions where
#         local density > 25th percentile of density in the region.

#         Parameters
#         ----------
#         frame : int
#             Frame number
#         species : int
#             Particle species (1 for positrons, 2 for electrons)
#         region : dict, optional
#             Region to select from: {'x_range': [min,max], 'y_range': [min,max]}

#         Returns
#         -------
#         list of dict
#             List of test particle data dictionaries
#         """
#         def fetch_var_at_step(out_dir, var, step, rad = False):
#             """
#             This wrapper function fetches either particle or field data at any output step (!NOT! simulation time step) 
#             from an output directory. Obviously.
#             """
#             filename = out_dir + var + '/' + var + '.tot.%05i'%step
#             if var == "prtl":

#                 return(isolde.getParticles(filename))
#             elif var == "flds":

#                 return(isolde.getFields(filename))

#             elif var == "spec":        
#                 return(isolde.getSpectra(filename, radiation = rad))

#             else:
#                 print("Not supported yet!")
#                 return False

#         # Load particle data
#         d = Data(TristanV2, steps=[frame], path=self.out_dir, 
#                  cfg_fname=self.input_file_name, params=True)
                
#         selected_particles = []

#         # Load density data
#         for species in [1, 2]:
#             particles = d.particles[species].sel(t=frame)
#             dens_key = "dens1" if species == 1 else "dens2"
#             dens_data = fetch_var_at_step("/scratch/10446/anindya_12/tristan_mp_v2/vault/output_psi0.5_mul1_3_mul2_0.1_TT_1.9e-4_rad_drag/", "flds", frame//5)[dens_key]
#             print("Density of species {} is loaded".format(species))
#             dens_avg = np.average(dens_data, axis=0)  # average over z or time axis

#             # Define region mask (if region is provided)
#             if region is not None:
#                 x_mask = ((particles.x.values >= region['x_range'][0]) & 
#                           (particles.x.values <= region['x_range'][1]))
#                 y_mask = ((particles.y.values >= region['y_range'][0]) & 
#                           (particles.y.values <= region['y_range'][1]))
#                 region_mask = x_mask & y_mask
#             else:
#                 region_mask = np.ones_like(particles.x.values, dtype=bool)

#             # Compute density threshold (25th percentile)
#             region_dens = dens_avg
#             if region is not None:
#                 # Assuming dens_avg has grid coordinates that can be indexed
#                 x_min, x_max = region['x_range']
#                 y_min, y_max = region['y_range']
#                 region_dens = dens_avg[y_min:y_max, x_min:x_max]  # adjust indexing as needed

#             dens_threshold = np.percentile(region_dens, 25)

#             # Interpolate local density for each particle
#             # Assuming dens_avg.shape corresponds to simulation grid (Ny, Nx)
#             # and particle positions are in normalized or grid units
#             Ny, Nx = dens_avg.shape
#             x_idx = np.clip(particles.x.values.astype(int), 0, Nx - 1)
#             y_idx = np.clip(particles.y.values.astype(int), 0, Ny - 1)

#             # Local density at each particle position
#             particle_local_dens = dens_avg[y_idx, x_idx]
            
#             # Filter particles based on density > 25th percentile
#             high_dens_mask = (particle_local_dens > dens_threshold) & region_mask

#             # Select and store particles            
            
#             for idx in np.where(high_dens_mask)[0]:
#                 selected_particles.append({
#                     'id': particles.idx.values[idx],
#                     'x': particles.x.values[idx],
#                     'y': particles.y.values[idx],
#                     'z': particles.z.values[idx],
#                     'u': particles.u.values[idx],
#                     'v': particles.v.values[idx],
#                     'w': particles.w.values[idx],
#                     'species': species,
#                     'mass': 1 * self.weight_fac,
#                     'charge': self.unit_ch * self.weight_fac if species == 1 else -self.unit_ch,
#                     'density': float(particle_local_dens[idx])  # store local density value
#                 })
#         print("{} particles selected".format(len(selected_particles)))
#         return selected_particles

    
#     def delta_s2(self, x4_1, x4_2):
#         """Calculate spacetime interval squared using Minkowski metric"""
#         x_diff = x4_1 - x4_2
#         return np.einsum('ij,jk,ik->i', x_diff, self.metric, x_diff)

#     def calculate_field_contribution_from_frame(self, test_particle, observation_frame, past_frame):
#         """
#         Calculate E and B field contribution from a single past frame
        
#         Parameters:
#         -----------
#         test_particle : dict
#             Test particle data at observation frame
#         observation_frame : int
#             Frame when test particle is observed
#         past_frame : int
#             Past frame to calculate field contribution from
            
#         Returns:
#         --------
#         tuple : (E_field_contribution, B_field_contribution) as 3D vectors
#         """
        
#         # Test particle 4-vector position at observation time
#         test_pos_3d = np.array([test_particle['x'], test_particle['y'], test_particle['z']])
#         x4_test = np.array([self.CC * observation_frame * self.interval, 
#                            test_particle['x'], test_particle['y'], test_particle['z']])
        
#         E_contribution = np.zeros(3)
#         B_contribution = np.zeros(3)
        
#         # Load particle data for this past frame and the next one
#         try:
# #             print("Field calculation...inside try")
#             d_curr = Data(TristanV2, steps=[past_frame], path=self.out_dir,
#                          cfg_fname=self.input_file_name, params=True)
#             d_next = Data(TristanV2, steps=[past_frame + 1], path=self.out_dir,
#                          cfg_fname=self.input_file_name, params=True)
# #             print("Data loaded!")
#         except:
#             return E_contribution, B_contribution  # Return zeros if data not available
#             print("Loading failed :( ")
        
#         # Process both species
#         for species in [1, 2]:
#             particles_curr = d_curr.particles[species].sel(t=past_frame)
#             particles_next = d_next.particles[species].sel(t=past_frame + 1)
            
#             # Find common particle IDs
#             current_ids = set(particles_curr.idx.values)
#             next_ids = set(particles_next.idx.values)
#             common_ids = list(current_ids.intersection(next_ids))
            
#             if len(common_ids) == 0:
#                 continue
                
#             # # Remove test particle from calculation if in same species
#             # if test_particle['species'] == species:
#             #     common_ids = [pid for pid in common_ids if pid != test_particle['id']]
            
#             if len(common_ids) == 0:
#                 continue
            
#             # Create index mappings
#             curr_id_to_idx = {id_val: i for i, id_val in enumerate(particles_curr.idx.values)}
#             next_id_to_idx = {id_val: i for i, id_val in enumerate(particles_next.idx.values)}
            
#             curr_indices = np.array([curr_id_to_idx[idx] for idx in common_ids])
#             next_indices = np.array([next_id_to_idx[idx] for idx in common_ids])
            
#             # Construct 4-vector positions for light cone calculation
#             curr_locs4 = np.column_stack((self.CC * past_frame * self.interval * np.ones_like(curr_indices),
#                                          particles_curr.x.values[curr_indices],
#                                          particles_curr.y.values[curr_indices],
#                                          particles_curr.z.values[curr_indices]))
            
#             next_locs4 = np.column_stack((self.CC * (past_frame + 1) * self.interval * np.ones_like(next_indices),
#                                          particles_next.x.values[next_indices],
#                                          particles_next.y.values[next_indices],
#                                          particles_next.z.values[next_indices]))
            
#             # Calculate spacetime intervals to test particle at observation time
#             ds2_curr = self.delta_s2(curr_locs4, x4_test)
#             ds2_next = self.delta_s2(next_locs4, x4_test)
            
#             # Light cone intersection condition: ds2 changes sign between frames
#             light_cone_mask = (ds2_curr * ds2_next) <= 0
            
#             if not np.any(light_cone_mask):
#                 continue
            
#             # Get indices of particles crossing light cone
#             lc_indices = np.where(light_cone_mask)[0]
#             lc_curr_indices = curr_indices[lc_indices]
#             lc_next_indices = next_indices[lc_indices]
            
#             # Get 3D positions of light cone particles
#             lc_curr_pos = np.column_stack((particles_curr.x.values[lc_curr_indices],
#                                          particles_curr.y.values[lc_curr_indices],
#                                          particles_curr.z.values[lc_curr_indices]))
#             lc_next_pos = np.column_stack((particles_next.x.values[lc_next_indices],
#                                          particles_next.y.values[lc_next_indices],
#                                          particles_next.z.values[lc_next_indices]))
            
#             # Average positions and calculate distances
#             lc_avg_pos = (lc_curr_pos + lc_next_pos) / 2
#             lc_displacements = test_pos_3d - lc_avg_pos
#             lc_distances = np.linalg.norm(lc_displacements, axis=1)
            
#             # Calculate unit vectors from particles to test particle
#             lc_unit_vectors = lc_displacements / lc_distances[:, np.newaxis]
            
#             # Get 4-velocities of light cone particles
#             lc_curr_4vel = np.column_stack((particles_curr.u.values[lc_curr_indices],
#                                           particles_curr.v.values[lc_curr_indices],
#                                           particles_curr.w.values[lc_curr_indices]))
#             lc_next_4vel = np.column_stack((particles_next.u.values[lc_next_indices],
#                                           particles_next.v.values[lc_next_indices],
#                                           particles_next.w.values[lc_next_indices]))
            
#             # Convert to beta (v/c)
#             lc_curr_beta = lc_curr_4vel / np.sqrt(1 + np.linalg.norm(lc_curr_4vel, axis=1)**2)[:, np.newaxis]
#             lc_next_beta = lc_next_4vel / np.sqrt(1 + np.linalg.norm(lc_next_4vel, axis=1)**2)[:, np.newaxis]
#             lc_avg_beta = (lc_curr_beta + lc_next_beta) / 2
#             lc_avg_gamma = 1 / np.sqrt(1 - np.linalg.norm(lc_avg_beta, axis = 1)**2)
#             lc_beta_dot = (lc_next_beta - lc_curr_beta) / self.interval
            
#             # Calculate electric field contribution (Liénard-Wiechert fields)
#             charge_sign = 1 if species == 1 else -1
#             E_acc_num = (charge_sign * self.weight_fac * self.unit_ch / self.CC * 
#                     np.cross(lc_unit_vectors, np.cross(lc_unit_vectors - lc_avg_beta, lc_beta_dot)))
#             E_acc_denom = ((1 - np.einsum('ij,ij->i', lc_unit_vectors, lc_avg_beta))**3 * lc_distances)
            
#             E_vel_num = charge_sign * self.weight_fac * self.unit_ch * (lc_unit_vectors - lc_avg_beta)
#             E_vel_denom = E_acc_denom * lc_distances * lc_avg_gamma**2 

#             # Avoid division by zero
#             valid_mask = (E_acc_denom != 0) & (E_vel_denom != 0)
#             if np.any(valid_mask):
#                 E_fields = np.zeros_like(E_acc_num)
#                 E_fields[valid_mask] = E_acc_num[valid_mask] / E_acc_denom[valid_mask, np.newaxis] + E_vel_num[valid_mask] / E_vel_denom[valid_mask, np.newaxis]
                
#                 # Calculate magnetic field
#                 B_fields = np.cross(lc_unit_vectors, E_fields)
                
#                 # Sum contributions from this frame
#                 E_contribution += np.sum(E_fields, axis=0)
#                 B_contribution += np.sum(B_fields, axis=0)
        
#         return E_contribution * self.stride, B_contribution * self.stride
    
#     def calculate_forces(self, test_particle, E_field, B_field):
#         """
#         Calculate radiation reaction force using Faraday tensor
        
#         Parameters:
#         -----------
#         test_particle : dict
#             Test particle data
#         E_field, B_field : array
#             Electric and magnetic field vectors
            
#         Returns:
#         --------
#         array : 4-force from radiation reaction
#         """
#         # Test particle 4-velocity
#         u_4vel = np.array([test_particle['u'], test_particle['v'], test_particle['w']])
#         gamma = np.sqrt(1 + np.linalg.norm(u_4vel)**2)
#         u_mu = np.array([gamma, test_particle['u'], test_particle['v'], test_particle['w']]) # u^\mu
        
#         # Construct Faraday tensor F_μν
#         F_tensor = np.zeros((4, 4))
#         F_tensor[0, 1:4] = E_field
#         F_tensor[1:4, 0] = -E_field
#         F_tensor[1, 2] = B_field[2]
#         F_tensor[2, 1] = -B_field[2]
#         F_tensor[1, 3] = -B_field[1]
#         F_tensor[3, 1] = B_field[1]
#         F_tensor[2, 3] = B_field[0]
#         F_tensor[3, 2] = -B_field[0]

#         F_contra_tensor = np.einsum('ij,jk,kl->il', self.metric, F_tensor, self.metric)
        
#         # Calculate F_μν u^ν
#         F_dot_u = np.einsum('ij,j->i', F_tensor, u_mu)
        
#         # Calculate (F_αβ u^β)^2
#         F_dot_u_squared = np.einsum('i,ij,j', F_dot_u, self.metric, F_dot_u)

#         # Calculate F^αβ F_βν u^ν
#         F_dot_F_dot_u = np.einsum('ij,j->i', F_contra_tensor, F_dot_u)
        
#         # Radiation reaction 4-force: (2e^4/3m^2)[F^αβ F_βν u^ν + (F_μν u^ν)² u^α]        
#         mass_squared = test_particle['mass']**2  # Normalized units
#         prefactor = 2 * (test_particle['charge'])**4 / (3 * mass_squared)
        
#         f_rad = prefactor * (F_dot_F_dot_u + F_dot_u_squared * u_mu)

#         f_lorentz = test_particle['charge'] * np.einsum('ij,j->i', self.metric, F_dot_u)
        
#         return f_rad, f_lorentz
    


# # Example placeholder, updated by wrapper
# sim_params = {'out_dir': '/scratch/10446/anindya_12/tristan_mp_v2/vault/output_psi0.5_mul1_3_mul2_0.1_TT_1.9e-4_rad_drag_prtl_maxxed/', 'input_file_name': '/scratch/10446/anindya_12/tristan_mp_v2/vault/output_psi0.5_mul1_3_mul2_0.1_TT_1.9e-4_rad_drag_prtl_maxxed/temp_input_psi0.5_mul13_mul20.1_TT1.9e-4.in', 'bounds_x': [10220, 10280], 'bounds_y': [1020, 1050], 'interval': 2.0, 'omegap0': 0.022281692032865345, 'wA_wp': 0.01, 'weight_fac': 0.0020473518, 'unit_ch': 0.028248018259818944, 'CC': 0.35, 'stride': 10.0}
