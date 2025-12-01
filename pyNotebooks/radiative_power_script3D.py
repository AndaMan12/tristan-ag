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
import pandas as pd


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

# Define extended bounds for 3D duplication
# Original 2D simulation bounds
original_bounds_y = [0, grid_y]
original_bounds_z = [0, 0]  # 2D simulation has no z extent

# Extended bounds for 3D analysis - MODIFY THESE AS NEEDED
extended_bounds_y = [0, 2 * np.pi / mode]  # Extend y by an Alfven wavelength
extended_bounds_z = [0, 2 * np.pi / mode]  # Extend z by an Alfven wavelength

# Helper functions
bounds_x = [fin_x_boundary, grid_x]
bounds_y = extended_bounds_y  # Use extended bounds
bounds_z = extended_bounds_z  # Use extended bounds

# def duplicate_particles_3d(particles, original_y_size, original_z_size, extended_y_bounds, extended_z_bounds):
#     """
#     Duplicate 2D particle data to create 3D periodic structure
#     """
#     if particles is None or len(particles.x.values) == 0:
#         return particles
    
#     # Get original particle data
#     orig_x = particles.x.values
#     orig_y = particles.y.values
#     orig_z = particles.z.values if hasattr(particles, 'z') else np.zeros_like(orig_x)
#     orig_u = particles.u.values
#     orig_v = particles.v.values
#     orig_w = particles.w.values
#     orig_idx = particles.idx.values
    
#     # Calculate duplication parameters
#     y_period = original_y_size
#     z_period = original_z_size  # Use same period for z as y for simplicity
    
#     # Number of copies in each direction
#     n_y_copies = int(np.ceil((extended_y_bounds[1] - extended_y_bounds[0]) / y_period))
#     n_z_copies = int(np.ceil((extended_z_bounds[1] - extended_z_bounds[0]) / z_period))
    
#     # Lists to store duplicated data
#     all_x, all_y, all_z = [], [], []
#     all_u, all_v, all_w = [], [], []
#     all_idx = []
    
#     idx_offset = 0
    
#     # Duplicate across y and z
#     for i_z in range(n_z_copies):
#         z_offset = i_z * z_period + extended_z_bounds[0]
#         for i_y in range(n_y_copies):
#             y_offset = i_y * y_period + extended_y_bounds[0]
            
#             # New positions
#             new_x = orig_x.copy()
#             new_y = orig_y + y_offset
#             new_z = orig_z + z_offset
            
#             # Filter to keep only particles within extended bounds
#             valid_mask = ((new_y >= extended_y_bounds[0]) & (new_y <= extended_y_bounds[1]) &
#                          (new_z >= extended_z_bounds[0]) & (new_z <= extended_z_bounds[1]))
            
#             if np.any(valid_mask):
#                 all_x.append(new_x[valid_mask])
#                 all_y.append(new_y[valid_mask])
#                 all_z.append(new_z[valid_mask])
#                 all_u.append(orig_u[valid_mask])
#                 all_v.append(orig_v[valid_mask])
#                 all_w.append(orig_w[valid_mask])
#                 # Create unique indices for duplicated particles
#                 all_idx.append(orig_idx[valid_mask] + idx_offset)
#                 idx_offset += len(orig_idx)
    
#     if not all_x:
#         return particles
    
#     # Concatenate all duplicated data
#     particles_dict = {
#         'x': np.concatenate(all_x),
#         'y': np.concatenate(all_y),
#         'z': np.concatenate(all_z),
#         'u': np.concatenate(all_u),
#         'v': np.concatenate(all_v),
#         'w': np.concatenate(all_w),
#         'idx': np.concatenate(all_idx)
#     }
    
#     # Create a new particles object with the duplicated data
#     # Note: This assumes particles has a similar structure - may need adjustment
#     import pandas as pd
#     duplicated_particles = pd.DataFrame(particles_dict)
    
#     # Copy attributes from original particles object
#     for attr in dir(particles):
#         if not attr.startswith('_') and attr not in particles_dict:
#             try:
#                 setattr(duplicated_particles, attr, getattr(particles, attr))
#             except:
#                 pass
    
#     return duplicated_particles

"""
Parallelised + shared memory duplication
"""


def process_duplication_chunk(chunk_data, orig_x, orig_y, orig_z, orig_u, orig_v, orig_w, orig_idx, 
                             y_period, z_period, extended_y_bounds, extended_z_bounds):
    """
    Process a chunk of (i_z, i_y) combinations for particle duplication
    """
    chunk_results = {
        'x': [], 'y': [], 'z': [], 'u': [], 'v': [], 'w': [], 'idx': []
    }
    
    for i_z, i_y in chunk_data:
        z_offset = i_z * z_period + extended_z_bounds[0]
        y_offset = i_y * y_period + extended_y_bounds[0]
        
        # New positions
        new_x = orig_x.copy()
        new_y = orig_y + y_offset
        new_z = orig_z + z_offset
        
        # Filter to keep only particles within extended bounds
        valid_mask = ((new_y >= extended_y_bounds[0]) & (new_y <= extended_y_bounds[1]) &
                     (new_z >= extended_z_bounds[0]) & (new_z <= extended_z_bounds[1]))
        
        if np.any(valid_mask):
            chunk_results['x'].append(new_x[valid_mask])
            chunk_results['y'].append(new_y[valid_mask])
            chunk_results['z'].append(new_z[valid_mask])
            chunk_results['u'].append(orig_u[valid_mask])
            chunk_results['v'].append(orig_v[valid_mask])
            chunk_results['w'].append(orig_w[valid_mask])
            # Create unique indices - we'll handle global uniqueness later
            chunk_results['idx'].append(orig_idx[valid_mask])
    
    return chunk_results

def create_chunks(n_z_copies, n_y_copies, n_processes):
    """
    Create balanced chunks of (i_z, i_y) combinations for parallel processing
    """
    # Create all (i_z, i_y) combinations
    all_combinations = [(i_z, i_y) for i_z in range(n_z_copies) for i_y in range(n_y_copies)]
    
    # Split into chunks
    chunk_size = len(all_combinations) // n_processes
    if chunk_size == 0:
        chunk_size = 1
    
    chunks = []
    for i in range(0, len(all_combinations), chunk_size):
        chunks.append(all_combinations[i:i + chunk_size])
    
    return chunks

# Usage example - replace your existing function call:
# OLD: particles_curr_3d = duplicate_particles_3d(particles_curr, grid_y, 1, extended_bounds_y, extended_bounds_z)
# NEW: particles_curr_3d = duplicate_particles_3d_parallel(particles_curr, grid_y, 1, extended_bounds_y, extended_bounds_z)

def duplicate_particles_3d_parallel(particles, original_y_size, original_z_size, 
                                   extended_y_bounds, extended_z_bounds, n_processes=None):
    """
    Duplicate 2D particle data to create 3D periodic structure using parallel processing
    """
    if particles is None or len(particles.x.values) == 0:
        return particles
    
    if n_processes is None:
        n_processes = min(n_cpus, 8)  # Cap at 8 to avoid too much overhead
    
    # Get original particle data
    orig_x = particles.x.values
    orig_y = particles.y.values
    orig_z = particles.z.values if hasattr(particles, 'z') else np.zeros_like(orig_x)
    orig_u = particles.u.values
    orig_v = particles.v.values
    orig_w = particles.w.values
    orig_idx = particles.idx.values
    
    # Calculate duplication parameters
    y_period = original_y_size
    z_period = original_z_size  # Use same period for z as y for simplicity
    
    # Number of copies in each direction
    n_y_copies = int(np.ceil((extended_y_bounds[1] - extended_y_bounds[0]) / y_period))
    n_z_copies = int(np.ceil((extended_z_bounds[1] - extended_z_bounds[0]) / z_period))
    
    total_combinations = n_z_copies * n_y_copies
    print(f"Processing {total_combinations} combinations using {n_processes} processes")
    
    # Create chunks for parallel processing
    chunks = create_chunks(n_z_copies, n_y_copies, n_processes)
    
    # Create partial function with fixed arguments
    process_func = partial(
        process_duplication_chunk,
        orig_x=orig_x, orig_y=orig_y, orig_z=orig_z,
        orig_u=orig_u, orig_v=orig_v, orig_w=orig_w, orig_idx=orig_idx,
        y_period=y_period, z_period=z_period,
        extended_y_bounds=extended_y_bounds, extended_z_bounds=extended_z_bounds
    )
    
    # Process chunks in parallel
    with Pool(n_processes) as pool:
        chunk_results = list(tqdm(
            pool.imap(process_func, chunks),
            total=len(chunks),
            desc="Duplicating particles"
        ))
    
    # Combine results from all chunks
    all_x, all_y, all_z = [], [], []
    all_u, all_v, all_w = [], [], []
    all_idx = []
    
    idx_offset = 0
    for chunk_result in chunk_results:
        for i in range(len(chunk_result['x'])):
            all_x.append(chunk_result['x'][i])
            all_y.append(chunk_result['y'][i])
            all_z.append(chunk_result['z'][i])
            all_u.append(chunk_result['u'][i])
            all_v.append(chunk_result['v'][i])
            all_w.append(chunk_result['w'][i])
            # Create globally unique indices
            all_idx.append(chunk_result['idx'][i] + idx_offset)
            idx_offset += len(orig_idx)
    
    if not all_x:
        return particles
    
    # Concatenate all duplicated data
    particles_dict = {
        'x': np.concatenate(all_x),
        'y': np.concatenate(all_y),
        'z': np.concatenate(all_z),
        'u': np.concatenate(all_u),
        'v': np.concatenate(all_v),
        'w': np.concatenate(all_w),
        'idx': np.concatenate(all_idx)
    }
    
    # Create a new particles object with the duplicated data
    duplicated_particles = pd.DataFrame(particles_dict)
    
    # Copy attributes from original particles object
    for attr in dir(particles):
        if not attr.startswith('_') and attr not in particles_dict:
            try:
                setattr(duplicated_particles, attr, getattr(particles, attr))
            except:
                pass
    
    return duplicated_particles

# Alternative version using shared memory for very large datasets
def duplicate_particles_3d_shared_memory(particles, original_y_size, original_z_size, 
                                        extended_y_bounds, extended_z_bounds, n_processes=None):
    """
    Version using shared memory - useful for very large particle datasets
    """
    try:
        from multiprocessing import shared_memory
    except ImportError:
        print("Shared memory not available, falling back to regular parallel version")
        return duplicate_particles_3d_parallel(particles, original_y_size, original_z_size, 
                                             extended_y_bounds, extended_z_bounds, n_processes)
    
    if particles is None or len(particles.x.values) == 0:
        return particles
    
    if n_processes is None:
        n_processes = min(n_cpus, 8)
    
    # Get original particle data
    orig_x = particles.x.values.astype(np.float64)
    orig_y = particles.y.values.astype(np.float64)
    orig_z = particles.z.values.astype(np.float64) if hasattr(particles, 'z') else np.zeros_like(orig_x, dtype=np.float64)
    orig_u = particles.u.values.astype(np.float64)
    orig_v = particles.v.values.astype(np.float64)
    orig_w = particles.w.values.astype(np.float64)
    orig_idx = particles.idx.values.astype(np.int64)
    
    # Create shared memory arrays
    shm_x = shared_memory.SharedMemory(create=True, size=orig_x.nbytes)
    shm_y = shared_memory.SharedMemory(create=True, size=orig_y.nbytes)
    shm_z = shared_memory.SharedMemory(create=True, size=orig_z.nbytes)
    shm_u = shared_memory.SharedMemory(create=True, size=orig_u.nbytes)
    shm_v = shared_memory.SharedMemory(create=True, size=orig_v.nbytes)
    shm_w = shared_memory.SharedMemory(create=True, size=orig_w.nbytes)
    shm_idx = shared_memory.SharedMemory(create=True, size=orig_idx.nbytes)
    
    # Copy data to shared memory
    shared_x = np.ndarray(orig_x.shape, dtype=np.float64, buffer=shm_x.buf)
    shared_y = np.ndarray(orig_y.shape, dtype=np.float64, buffer=shm_y.buf)
    shared_z = np.ndarray(orig_z.shape, dtype=np.float64, buffer=shm_z.buf)
    shared_u = np.ndarray(orig_u.shape, dtype=np.float64, buffer=shm_u.buf)
    shared_v = np.ndarray(orig_v.shape, dtype=np.float64, buffer=shm_v.buf)
    shared_w = np.ndarray(orig_w.shape, dtype=np.float64, buffer=shm_w.buf)
    shared_idx = np.ndarray(orig_idx.shape, dtype=np.int64, buffer=shm_idx.buf)
    
    shared_x[:] = orig_x
    shared_y[:] = orig_y
    shared_z[:] = orig_z
    shared_u[:] = orig_u
    shared_v[:] = orig_v
    shared_w[:] = orig_w
    shared_idx[:] = orig_idx
    
    try:
        # Calculate duplication parameters
        y_period = original_y_size
        z_period = original_z_size
        
        n_y_copies = int(np.ceil((extended_y_bounds[1] - extended_y_bounds[0]) / y_period))
        n_z_copies = int(np.ceil((extended_z_bounds[1] - extended_z_bounds[0]) / z_period))
        
        # Create chunks and process
        chunks = create_chunks(n_z_copies, n_y_copies, n_processes)
        
        # For shared memory version, we'd need to modify the processing function
        # This is more complex and may not provide significant benefits for most use cases
        # For now, fall back to regular parallel version
        result = duplicate_particles_3d_parallel(particles, original_y_size, original_z_size, 
                                                extended_y_bounds, extended_z_bounds, n_processes)
        
    finally:
        # Clean up shared memory
        shm_x.close()
        shm_x.unlink()
        shm_y.close()
        shm_y.unlink()
        shm_z.close()
        shm_z.unlink()
        shm_u.close()
        shm_u.unlink()
        shm_v.close()
        shm_v.unlink()
        shm_w.close()
        shm_w.unlink()
        shm_idx.close()
        shm_idx.unlink()
    
    return result

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
    # Apply 3D duplication to both current and next frame particles
    # particles_curr_3d = duplicate_particles_3d(particles_curr, grid_y, 1, extended_bounds_y, extended_bounds_z)
    # particles_next_3d = duplicate_particles_3d(particles_next, grid_y, 1, extended_bounds_y, extended_bounds_z)
    particles_curr_3d = duplicate_particles_3d_parallel(particles_curr, grid_y, 1, extended_bounds_y, extended_bounds_z)
    particles_next_3d = duplicate_particles_3d_parallel(particles_next, grid_y, 1, extended_bounds_y, extended_bounds_z)
    
    curr_ids = particles_curr_3d.idx.values
    next_ids = particles_next_3d.idx.values
    curr_mask = np.isin(curr_ids, next_ids)
    if not np.any(curr_mask):
        return np.zeros(3, dtype=np.float64), np.zeros(3, dtype=np.float64)
    
    curr_x = particles_curr_3d.x.values[curr_mask].astype(np.float64)
    curr_y = particles_curr_3d.y.values[curr_mask].astype(np.float64)
    curr_z = particles_curr_3d.z.values[curr_mask].astype(np.float64)
    curr_u = particles_curr_3d.u.values[curr_mask].astype(np.float64)
    curr_v = particles_curr_3d.v.values[curr_mask].astype(np.float64)
    curr_w = particles_curr_3d.w.values[curr_mask].astype(np.float64)
    
    next_mask = np.isin(next_ids, curr_ids[curr_mask])
    next_x = particles_next_3d.x.values[next_mask].astype(np.float64)
    next_y = particles_next_3d.y.values[next_mask].astype(np.float64)
    next_z = particles_next_3d.z.values[next_mask].astype(np.float64)
    next_u = particles_next_3d.u.values[next_mask].astype(np.float64)
    next_v = particles_next_3d.v.values[next_mask].astype(np.float64)
    next_w = particles_next_3d.w.values[next_mask].astype(np.float64)
    
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

def process_frames_sequentially(frames, out_dir, input_file_name, x4_obs, interval, CC, 
                               omegap0, wA_wp, unit_ch, weight_fac):
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

def compute_radiation_fields_optimized(out_dir, input_file_name, frame_f, frame_i, 
                                     interval, CC, omegap0, wA_wp, unit_ch, weight_fac, 
                                     x4_obs, stride):
    frames = list(range(frame_f, frame_i))
    print(f"Processing {len(frames)} frames")
    
    results = process_frames_sequentially(
        frames, out_dir, input_file_name, x4_obs, interval, CC, 
        omegap0, wA_wp, unit_ch, weight_fac
    )
    
    E_tot = np.zeros(3, dtype=np.float64)
    B_tot = np.zeros(3, dtype=np.float64)
    for E_frame, B_frame in results:
        E_tot += E_frame
        B_tot += B_frame
    
    E_tot *= stride
    B_tot *= stride
    
    return E_tot, B_tot

def main_radiation_calculation(r, theta, phi):
    x_origin = np.array([sum(bounds_x)/2, sum(bounds_y)/2, sum(bounds_z)/2], dtype=np.float64)
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
        stride=stride
    )    
    
    Poynting_vec = CC / (4*np.pi) * np.cross(E_tot, B_tot)
    print("Poynting vector:", Poynting_vec, " from x_obs = ", x_obs)
    return Poynting_vec

# Spherical integration setup
r = 10 * max(bounds_x[1], bounds_y[1], bounds_z[1])  # Use extended bounds for radius
n = 3  # Lebedev quadrature order
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
    V_i = main_radiation_calculation(r, theta[i], phi[i])
    F_r_local[idx] = np.dot(V_i, x[i])  # Radial component of Poynting vector
    local_flux += w[i] * F_r_local[idx]

# Gather data to rank 0
all_flux = comm.gather(local_flux, root=0)
theta_all = comm.gather(theta_local, root=0)
phi_all = comm.gather(phi_local, root=0)
F_r_all = comm.gather(F_r_local, root=0)

# Process and save data on rank 0
if rank == 0:
    flux = sum(all_flux)
    
    # Calculate volume scaling factor for per-unit-volume normalization
    original_volume = (bounds_x[1] - bounds_x[0]) * (original_bounds_y[1] - original_bounds_y[0]) * 1  # 2D has unit z
    extended_volume = (bounds_x[1] - bounds_x[0]) * (bounds_y[1] - bounds_y[0]) * (bounds_z[1] - bounds_z[0])
    volume_ratio = extended_volume / original_volume
    
    flux *= (r * omegap0 * wA_wp / CC)**2 / volume_ratio  # Normalize by volume ratio
    print(f"Total flux per unit volume: {flux}")
    print(f"Volume scaling factor: {volume_ratio}")
    
    # Concatenate gathered arrays
    theta_all = np.concatenate(theta_all)
    phi_all = np.concatenate(phi_all)
    F_r_all = np.concatenate(F_r_all)
    
    # Save to .npz file for visualization
    np.savez('poynting_data_{}.npz'.format(n), theta=theta_all, phi=phi_all, poynting_radial=F_r_all)
    print("Poynting vector data saved to 'poynting_data_{}.npz'".format(n))
else:
    flux = None