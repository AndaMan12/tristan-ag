import numpy as np
import scipy as sp
import tristanVis.isolde as isolde
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from functools import partial
import warnings
from numba import jit
import os
from mpi4py import MPI

warnings.filterwarnings('ignore')

from graphet import Data
from graphet.plugins import TristanV2

# Set n_processes based on SLURM_CPUS_PER_TASK
if 'SLURM_CPUS_PER_TASK' in os.environ:
    n_cpus = int(os.environ['SLURM_CPUS_PER_TASK'])
else:
    n_cpus = cpu_count()

out_dir = "../vault/output_psi0.5_mul1_3_mul2_0.1_TT_1.9e-2_rad_drag_temp/"

# Input file reading
input_file_name = out_dir + "temp_input_psi0.5_mul13_mul20.1_TT1.9e-2.in"
input_params = isolde.parseInput(input_file_name)
m1 = input_params["particles"]["m1"]
m2 = input_params["particles"]["m2"]
interval = input_params["output"]["interval"]
lst_time = input_params["time"]["last"]
grid_x = int(input_params["grid"]["mx0"])
grid_y = int(input_params["grid"]["my0"])
Nsteps = int(lst_time // interval)
CC = input_params["algorithm"]["c"]
COMP = input_params["plasma"]["c_omp"]
ppc0 = input_params["plasma"]["ppc0"]
SIGMA = input_params["plasma"]["sigma"]
stride = input_params["output"]["stride"]
B_norm = CC**2 * SIGMA**0.5 / COMP
unit_ch = CC**2 / (ppc0 * COMP**2)
weight_fac = 1.9316248E-03
wA_wp = input_params["problem"]["wA_wp"]
B_0 = input_params["problem"]["B_0"]
B_amp = input_params["problem"]["B_amplitude"]
spread = input_params["problem"]["spread"]
mode = input_params["problem"]["mode"]
psi = input_params["problem"]["psi"]
mult1 = input_params["problem"]["multiplicity_1"]
mult2 = input_params["problem"]["multiplicity_2"]
TT = input_params["problem"]["temperature"]
ramp_width = input_params["problem"]["ramp_width"]
init_x_boundary = int(5 * np.pi / mode) + 1
fin_x_boundary = init_x_boundary + ramp_width * (2 * np.pi / mode)
ds = wA_wp * 2 * np.pi / mode
omegap0 = CC / ds

# Helper functions
bounds_x = [fin_x_boundary, grid_x]
bounds_y = [0, grid_y]
bounds_z = [0, 0]

def farthest_pt(x_o):
    farthest_dict = {
        "UFL": np.array([max(bounds_x), min(bounds_y), min(bounds_z)]),
        "UFR": np.array([min(bounds_x), min(bounds_y), min(bounds_z)]),
        "UBL": np.array([max(bounds_x), max(bounds_y), min(bounds_z)]),
        "UBR": np.array([min(bounds_x), max(bounds_y), min(bounds_z)]),
        "DFL": np.array([max(bounds_x), min(bounds_y), max(bounds_z)]),
        "DFR": np.array([min(bounds_x), min(bounds_y), max(bounds_z)]),
        "DBL": np.array([max(bounds_x), max(bounds_y), max(bounds_z)]),
        "DBR": np.array([min(bounds_x), max(bounds_y), max(bounds_z)])
    }
    sector_code = ""
    if x_o[2] > sum(bounds_z)/2:
        sector_code += "U"
    else:
        sector_code += "D"
    if x_o[1] > sum(bounds_y)/2:
        sector_code += "F"
    else:
        sector_code += "B"
    if x_o[0] > sum(bounds_x)/2:
        sector_code += "R"
    else:
        sector_code += "L"
    return farthest_dict[sector_code]

def nearest_pt(x_o):
    nearest_dict = {
        "UFL": np.array([min(bounds_x), max(bounds_y), max(bounds_z)]),
        "UFC": np.array([x_o[0]       , max(bounds_y), max(bounds_z)]),
        "UFR": np.array([max(bounds_x), max(bounds_y), max(bounds_z)]),
        "UCL": np.array([min(bounds_x), x_o[1]       , max(bounds_z)]),
        "UCC": np.array([x_o[0]       , x_o[1]       , max(bounds_z)]),
        "UCR": np.array([max(bounds_x), x_o[1]       , max(bounds_z)]),
        "UBL": np.array([min(bounds_x), min(bounds_y), max(bounds_z)]),
        "UBC": np.array([x_o[0]       , min(bounds_y), max(bounds_z)]),
        "UBR": np.array([max(bounds_x), min(bounds_y), max(bounds_z)]),

        "CFL": np.array([min(bounds_x), max(bounds_y), x_o[2]       ]),
        "CFC": np.array([x_o[0]       , max(bounds_y), x_o[2]       ]),
        "CFR": np.array([max(bounds_x), max(bounds_y), x_o[2]       ]),
        "CCL": np.array([min(bounds_x), x_o[1]       , x_o[2]       ]),
        "CCC": np.array([x_o[0]       , x_o[1]       , x_o[2]       ]),
        "CCR": np.array([max(bounds_x), x_o[1]       , x_o[2]       ]),
        "CBL": np.array([min(bounds_x), min(bounds_y), x_o[2]       ]),
        "CBC": np.array([x_o[0]       , min(bounds_y), x_o[2]       ]),
        "CBR": np.array([max(bounds_x), min(bounds_y), x_o[2]       ]),

        "DFL": np.array([min(bounds_x), max(bounds_y), min(bounds_z)]),
        "DFC": np.array([x_o[0]       , max(bounds_y), min(bounds_z)]),
        "DFR": np.array([max(bounds_x), max(bounds_y), min(bounds_z)]),
        "DCL": np.array([min(bounds_x), x_o[1]       , min(bounds_z)]),
        "DCC": np.array([x_o[0]       , x_o[1]       , min(bounds_z)]),
        "DCR": np.array([max(bounds_x), x_o[1]       , min(bounds_z)]),
        "DBL": np.array([min(bounds_x), min(bounds_y), min(bounds_z)]),
        "DBC": np.array([x_o[0]       , min(bounds_y), min(bounds_z)]),
        "DBR": np.array([max(bounds_x), min(bounds_y), min(bounds_z)])
    }
    sector_code = ""
    if x_o[2] > max(bounds_z):
        sector_code += "U"
    elif x_o[2] < min(bounds_z):
        sector_code += "D"
    else:
        sector_code += "C"
    if x_o[1] > max(bounds_y):
        sector_code += "F"
    elif x_o[1] < min(bounds_y):
        sector_code += "B"
    else:
        sector_code += "C"
    if x_o[0] > max(bounds_x):
        sector_code += "R"
    elif x_o[0] < min(bounds_x):
        sector_code += "L"
    else:
        sector_code += "C"
    return nearest_dict[sector_code]

metric = np.diag(np.array([-1, 1, 1, 1], dtype=np.float64))

def delta_s2(x4_1, x4_2):
    x_diff = x4_1 - x4_2
    return np.sum(x_diff * (metric @ x_diff.T).T, axis=1)

@jit(nopython=True)
def compute_fields_numba(dirs, beta_av, beta_dot, dists, charge_unit, weight_factor, CC):
    n = dirs.shape[0]
    E_fields = np.zeros((n, 3), dtype=np.float64)
    B_fields = np.zeros((n, 3), dtype=np.float64)
    field_const = weight_factor * charge_unit / (4 * np.pi) / CC
    for i in range(n):
        n_minus_beta = dirs[i] - beta_av[i]
        cross_inner = np.cross(n_minus_beta, beta_dot[i])
        E_fields[i] = field_const * np.cross(dirs[i], cross_inner)
        dot_product = dirs[i][0] * beta_av[i][0] + dirs[i][1] * beta_av[i][1] + dirs[i][2] * beta_av[i][2]
        denom = (1 - dot_product)**3 * dists[i]
        if abs(denom) > 1e-15:
            E_fields[i] /= denom
        else:
            E_fields[i] = 0.0
        B_fields[i] = np.cross(dirs[i], E_fields[i])
    return np.sum(E_fields, axis=0), np.sum(B_fields, axis=0)

def process_particle_type_optimized(particles_curr, particles_next, frame, interval, CC, x4_obs, 
                                  omegap0, wA_wp, charge_unit, weight_factor):
    curr_ids = particles_curr.idx.values
    next_ids = particles_next.idx.values
    curr_mask = np.isin(curr_ids, next_ids)
    if not np.any(curr_mask):
        return np.zeros(3, dtype=np.float64), np.zeros(3, dtype=np.float64)
    curr_x = particles_curr.x.values[curr_mask].astype(np.float64)
    curr_y = particles_curr.y.values[curr_mask].astype(np.float64)
    curr_z = particles_curr.y.values[curr_mask].astype(np.float64)
    curr_u = particles_curr.u.values[curr_mask].astype(np.float64)
    curr_v = particles_curr.v.values[curr_mask].astype(np.float64)
    curr_w = particles_curr.w.values[curr_mask].astype(np.float64)
    next_mask = np.isin(next_ids, curr_ids[curr_mask])
    next_x = particles_next.x.values[next_mask].astype(np.float64)
    next_y = particles_next.y.values[next_mask].astype(np.float64)
    next_z = particles_next.z.values[next_mask].astype(np.float64)
    next_u = particles_next.u.values[next_mask].astype(np.float64)
    next_v = particles_next.v.values[next_mask].astype(np.float64)
    next_w = particles_next.w.values[next_mask].astype(np.float64)
    n_particles = len(curr_x)
    curr_locs4 = np.column_stack((
        CC * frame * interval * np.ones(n_particles, dtype=np.float64),
        curr_x, curr_y, curr_z
    ))
    next_locs4 = np.column_stack((
        CC * (frame + 1) * interval * np.ones(n_particles, dtype=np.float64),
        next_x, next_y, next_z
    ))
    ds2_curr = delta_s2(curr_locs4, x4_obs)
    ds2_next = delta_s2(next_locs4, x4_obs)
    light_cone_mask = (ds2_curr * ds2_next) <= 0
    if not np.any(light_cone_mask):
        return np.zeros(3, dtype=np.float64), np.zeros(3, dtype=np.float64)
    lc_curr_locs3 = curr_locs4[light_cone_mask, 1:]
    lc_next_locs3 = next_locs4[light_cone_mask, 1:]
    lc_curr_4vel = np.column_stack((curr_u[light_cone_mask], curr_v[light_cone_mask], curr_w[light_cone_mask])).astype(np.float64)
    lc_next_4vel = np.column_stack((next_u[light_cone_mask], next_v[light_cone_mask], next_w[light_cone_mask])).astype(np.float64)
    lc_locs3_av = 0.5 * (lc_curr_locs3 + lc_next_locs3)
    lc_disps = x4_obs[1:].astype(np.float64) - lc_locs3_av
    lc_dists = np.sqrt(np.sum(lc_disps**2, axis=1))
    valid_dists = lc_dists > 1e-10
    if not np.any(valid_dists):
        return np.zeros(3, dtype=np.float64), np.zeros(3, dtype=np.float64)
    lc_disps = lc_disps[valid_dists]
    lc_dists = lc_dists[valid_dists]
    lc_curr_4vel = lc_curr_4vel[valid_dists]
    lc_next_4vel = lc_next_4vel[valid_dists]
    lc_dirs = lc_disps / lc_dists[:, np.newaxis]
    lc_curr_beta = lc_curr_4vel / np.sqrt(1 + np.sum(lc_curr_4vel**2, axis=1))[:, np.newaxis]
    lc_next_beta = lc_next_4vel / np.sqrt(1 + np.sum(lc_next_4vel**2, axis=1))[:, np.newaxis]
    lc_beta_av = 0.5 * (lc_curr_beta + lc_next_beta)
    lc_beta_dot = (lc_next_beta - lc_curr_beta) / (interval * omegap0 * wA_wp)
    lc_dists_dim_less = lc_dists * omegap0 * wA_wp / CC
    E_fields, B_fields = compute_fields_numba(lc_dirs, lc_beta_av, lc_beta_dot, lc_dists_dim_less, 
                                             charge_unit, weight_factor, CC)
    return E_fields, B_fields

def process_frame_pair(frame_data, x4_obs, interval, CC, omegap0, wA_wp, unit_ch, weight_fac):
    try:
        frame, particles_p_curr, particles_e_curr, particles_p_next, particles_e_next = frame_data
        if particles_p_curr is None or particles_e_curr is None:
            return np.zeros(3), np.zeros(3)
        E_frame = np.zeros(3, dtype=np.float64)
        B_frame = np.zeros(3, dtype=np.float64)
        if len(particles_p_curr.idx.values) > 0 and len(particles_p_next.idx.values) > 0:
            E_p, B_p = process_particle_type_optimized(
                particles_p_curr, particles_p_next, frame, interval, CC, x4_obs, 
                omegap0, wA_wp, unit_ch, weight_fac
            )
            E_frame += E_p
            B_frame += B_p
        if len(particles_e_curr.idx.values) > 0 and len(particles_e_next.idx.values) > 0:
            E_e, B_e = process_particle_type_optimized(
                particles_e_curr, particles_e_next, frame, interval, CC, x4_obs, 
                omegap0, wA_wp, -unit_ch, weight_fac
            )
            E_frame += E_e
            B_frame += B_e
        return E_frame, B_frame
    except Exception as e:
        print(f"Error processing frame {frame}: {e}")
        return np.zeros(3, dtype=np.float64), np.zeros(3, dtype=np.float64)

def load_frame_data(frame_range, out_dir, input_file_name):
    frame_data = {}
    print("Pre-loading frame data...")
    for frame in tqdm(frame_range, desc="Loading frames"):
        try:
            d = Data(
                TristanV2,
                steps=[frame],
                path=out_dir,    
                cfg_fname=input_file_name,
                params=True,
            )
            frame_data[frame] = {
                'particles_p': d.particles[1].sel(t=frame),
                'particles_e': d.particles[2].sel(t=frame)
            }
        except Exception as e:
            print(f"Warning: Could not load frame {frame}: {e}")
            frame_data[frame] = {'particles_p': None, 'particles_e': None}
    return frame_data

def compute_radiation_fields_optimized(out_dir, input_file_name, frame_f, frame_i, 
                                     interval, CC, omegap0, wA_wp, unit_ch, weight_fac, 
                                     x4_obs, stride, n_processes=None, use_preloading=True):
    frames = list(range(frame_f, frame_i))
    if n_processes is None:
        if 'SLURM_CPUS_PER_TASK' in os.environ:
            n_cpus = int(os.environ['SLURM_CPUS_PER_TASK'])
        else:
            n_cpus = cpu_count()
        n_processes = min(n_cpus, len(frames))
    print(f"Processing {len(frames)} frames using {n_processes} processes")
    if use_preloading and len(frames) < 100:
        frame_data_dict = load_frame_data(frames, out_dir, input_file_name)
        frame_pairs = []
        for i, frame in enumerate(frames[:-1]):
            curr_data = frame_data_dict[frame]
            next_data = frame_data_dict[frame + 1] if frame + 1 in frame_data_dict else frame_data_dict[frames[i + 1]]
            frame_pairs.append((
                frame,
                curr_data['particles_p'],
                curr_data['particles_e'],
                next_data['particles_p'],
                next_data['particles_e']
            ))
        process_func = partial(
            process_frame_pair,
            x4_obs=x4_obs,
            interval=interval,
            CC=CC,
            omegap0=omegap0,
            wA_wp=wA_wp,
            unit_ch=unit_ch,
            weight_fac=weight_fac
        )
        with Pool(processes=n_processes) as pool:
            results = pool.map(process_func, frame_pairs)
    else:
        results = process_frames_sequentially_parallel(
            frames, out_dir, input_file_name, x4_obs, interval, CC, 
            omegap0, wA_wp, unit_ch, weight_fac, n_processes
        )
    E_tot = np.zeros(3, dtype=np.float64)
    B_tot = np.zeros(3, dtype=np.float64)
    for E_frame, B_frame in results:
        E_tot += E_frame
        B_tot += B_frame
    E_tot *= stride
    B_tot *= stride
    return E_tot, B_tot

def process_frames_sequentially_parallel(frames, out_dir, input_file_name, x4_obs, interval, CC, 
                                       omegap0, wA_wp, unit_ch, weight_fac, n_processes):
    batch_size = max(1, len(frames) // n_processes)
    results = []
    d_init = Data(TristanV2, steps=[frames[0]], path=out_dir, cfg_fname=input_file_name, params=True)
    particles_p_curr = d_init.particles[1].sel(t=frames[0])
    particles_e_curr = d_init.particles[2].sel(t=frames[0])
    for frame in tqdm(frames[1:], desc="Processing frames"):
        try:
            d = Data(TristanV2, steps=[frame], path=out_dir, cfg_fname=input_file_name, params=True)
            particles_p_next = d.particles[1].sel(t=frame)
            particles_e_next = d.particles[2].sel(t=frame)
            frame_data = (frame, particles_p_curr, particles_e_curr, particles_p_next, particles_e_next)
            E_frame, B_frame = process_frame_pair(
                frame_data, x4_obs, interval, CC, omegap0, wA_wp, unit_ch, weight_fac
            )
            results.append((E_frame, B_frame))
            particles_p_curr = particles_p_next
            particles_e_curr = particles_e_next
        except Exception as e:
            print(f"Error processing frame {frame}: {e}")
            results.append((np.zeros(3, dtype=np.float64), np.zeros(3, dtype=np.float64)))
    return results

def main_radiation_calculation(r, theta, phi, n_processes=None):
    x_origin = np.array([sum(bounds_x)/2, sum(bounds_y)/2, 0], dtype=np.float64)
    x_obs = x_origin + r * np.array([np.sin(theta) * np.cos(phi), 
                                        np.sin(theta) * np.sin(phi), 
                                        np.cos(theta)], 
                                        dtype=np.float64)
    x_near = nearest_pt(x_obs)
    x_far = farthest_pt(x_obs)
    t_obs_start = np.linalg.norm(x_far - x_obs) / CC
    t_sim = lst_time
    t_obs_end = t_sim + np.linalg.norm(x_near - x_obs) / CC
    t_par = 0.5
    t_obs = (1 - t_par) * t_obs_start + t_par * t_obs_end
    x4_obs = np.array([CC * t_obs,] + list(x_obs), dtype=np.float64)
    t_contact_i = t_obs - np.linalg.norm(x_near - x_obs) / CC
    t_contact_f = t_obs - np.linalg.norm(x_far - x_obs) / CC
    # print("x_obs = ", x_obs)
    # print("t_obs = ", t_obs)    
    # print("t_contact_i = ", t_contact_i)
    # print("t_contact_f = ", t_contact_f)
    frame_i = int(t_contact_i // interval + 1)
    frame_f = int(t_contact_f // interval - 1)
    E_tot, B_tot = compute_radiation_fields_optimized(
        out_dir=out_dir,
        input_file_name=input_file_name,
        frame_f=frame_f,
        frame_i=frame_i,
        interval=interval,
        CC=CC,
        omegap0=omegap0,
        wA_wp=wA_wp,
        unit_ch=unit_ch,
        weight_fac=weight_fac,
        x4_obs=x4_obs,
        stride=stride,
        n_processes=n_processes,
        use_preloading=False
    )    
    Poynting_vec = CC / (4*np.pi) * np.cross(E_tot, B_tot)
    print("Poynting vector:", Poynting_vec, " from x_obs = ", x_obs)
    return Poynting_vec


# Spherical integration setup
r = 10 * grid_x  # Example radius, adjust as per your script
n = 47  # Lebedev quadrature order
x, w = sp.integrate.lebedev_rule(n)
x = x.T
theta = np.arccos(x[:, 2])
phi = np.mod(np.arctan2(x[:, 1], x[:, 0]), 2 * np.pi)

# MPI setup
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

total_points = len(w)
points_per_process = total_points // size
remainder = total_points % size

start = rank * points_per_process + min(rank, remainder)
end = start + points_per_process + (1 if rank < remainder else 0)

# Initialize local arrays for data collection
theta_local = theta[start:end]
phi_local = phi[start:end]
F_r_local = np.zeros(len(theta_local))
local_flux = 0.0

# Compute flux and collect Poynting vector data
for idx, i in enumerate(range(start, end)):
    V_i = main_radiation_calculation(r, theta[i], phi[i], n_processes=n_cpus)
    F_r_local[idx] = np.dot(V_i, x[i])  # Radial component of Poynting vector
    local_flux += w[i] * F_r_local[idx]

# Gather data to rank 0
all_flux = comm.gather(local_flux, root=0)
theta_all = comm.gather(theta_local, roots=0)
phi_all = comm.gather(phi_local, root=0)
F_r_all = comm.gather(F_r_local, root=0)

# Process and save data on rank 0
if rank == 0:
    flux = sum(all_flux)
    flux *= (r * omegap0 * wA_wp / CC)**2  # Adjust scaling as per your script
    print(f"Total flux: {flux}")
    
    # Concatenate gathered arrays
    theta_all = np.concatenate(theta_all)
    phi_all = np.concatenate(phi_all)
    F_r_all = np.concatenate(F_r_all)
    
    # Save to .npz file for visualization
    np.savez('poynting_data_{}.npz'.format(n), theta=theta_all, phi=phi_all, poynting_radial=F_r_all)
    print("Poynting vector data saved to 'poynting_data_{}.npz'".format(n))
else:
    flux = None