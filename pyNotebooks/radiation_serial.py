import os
import pickle
import numpy as np
from tqdm import tqdm
from radiation_utils import RadiationReactionCalculator
import tristanVis.isolde as isolde

# ================================================================
# 1. Setup parameters
# ================================================================
out_dir = "/scratch/10446/anindya_12/tristan_mp_v2/vault/output_psi0.5_mul1_3_mul2_0.1_TT_1.9e-4_rad_drag_prtl_maxxed/"
input_file_name = os.path.join(out_dir, "temp_input_psi0.5_mul13_mul20.1_TT1.9e-4.in")

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
weight_fac = 2.0473518E-03
ds = wA_wp * 2 * np.pi / mode
omegap0 = CC / ds

bounds_x = [10220, 10280]
bounds_y = [1020, 1050]
obs_frame = 2400 * 5
max_frames_back = 150 * 5

calc_params = {
    "out_dir": out_dir,
    "input_file_name": input_file_name,
    "bounds_x": bounds_x,
    "bounds_y": bounds_y,
    "interval": interval,
    "omegap0": omegap0,
    "wA_wp": wA_wp,
    "weight_fac": weight_fac,
    "unit_ch": unit_ch,
    "CC": CC,
    "stride": stride
}

# ================================================================
# 2. Prepare directories and calculator
# ================================================================
os.makedirs("partial_frames", exist_ok=True)
calc = RadiationReactionCalculator(**calc_params)

# ================================================================
# 3. Load or select test particles
# ================================================================
if os.path.exists("test_particles.pkl"):
    with open("test_particles.pkl", "rb") as f:
        test_particles = pickle.load(f)
    print(f"Loaded {len(test_particles)} test particles from file.")
else:
    # NOTE: if you want to use a different density directory, pass density_dir via get_test_particle
    test_particles = calc.get_test_particle(obs_frame, region={"x_range": bounds_x, "y_range": bounds_y}, 
                                            density_dir = "/scratch/10446/anindya_12/tristan_mp_v2/vault/output_psi0.5_mul1_3_mul2_0.1_TT_1.9e-4_rad_drag/")
    if not isinstance(test_particles, list):
        test_particles = [test_particles]
    with open("test_particles.pkl", "wb") as f:
        pickle.dump(test_particles, f)
    print(f"Saved {len(test_particles)} test particles to test_particles.pkl")

# ================================================================
# 4. Loop over past frames (serial)
# ================================================================
past_frames = list(range(obs_frame - max_frames_back, obs_frame))


# Initialize first pair
pf0 = past_frames[0]
d_curr, d_next = calc.load_frame_pair(pf0)
if d_curr is None or d_next is None:
    raise RuntimeError(f"Failed to load initial frame pair {pf0}/{pf0+1}")

print("First frame pair loaded!")
    
for i, pf in enumerate(past_frames):
    print("pf Loop started!")
    
    out_file = f"partial_frames/frame_{pf:06d}.pkl"
    if os.path.exists(out_file):
        # advance the sliding window before continuing
        if i < len(past_frames) - 1:
            # load the next frame pair incrementally
            next_pf = past_frames[i + 1]
            new_next = calc.load_particles_from_frame(next_pf + 1)
            d_curr, d_next = d_next, new_next
        continue

    if d_curr is None or d_next is None:
        print(f"Skipping frame {pf}: failed to load frame pair")
        continue

    results = []
    print("Empty results made")
    for tp in test_particles:
        print("Particle id {} starting:".format(tp["id"]))
        try:
            E_new, B_new = calc.calculate_field_contribution_from_frame(
                tp, obs_frame, pf, d_curr=d_curr, d_next=d_next
            )
        except Exception as e:
            print(f"Error at frame {pf} for particle {tp['id']}: {e}")
            continue

        results.append({
            "test_id": tp["id"],
            "frame": pf,
            "E": E_new,
            "B": B_new
        })
        print("Particle id {} ... done!".format(tp["id"]))

    with open(out_file + ".tmp", "wb") as f:
        pickle.dump(results, f)
    os.replace(out_file + ".tmp", out_file)

    # Prepare next iteration
    if i < len(past_frames) - 1:
        next_pf = past_frames[i + 1]
        try:
            new_next = calc.load_particles_from_frame(next_pf + 1)
            d_curr, d_next = d_next, new_next
        except Exception as e:
            print(f"Warning: failed to load frame {next_pf + 1}: {e}")
            d_curr, d_next = None, None

print("All frame contributions computed.")

# ================================================================
# 5. Aggregate results
# ================================================================
E_totals = {tp["id"]: np.zeros(3) for tp in test_particles}
B_totals = {tp["id"]: np.zeros(3) for tp in test_particles}

frame_files = sorted(os.listdir("partial_frames"))
for ff in tqdm(frame_files, desc="Aggregating"):
    with open(os.path.join("partial_frames", ff), "rb") as f:
        frame_results = pickle.load(f)
    for res in frame_results:
        tid = res["test_id"]
        E_totals[tid] += res["E"]
        B_totals[tid] += res["B"]

# Compute final forces
results = {}
for tp in test_particles:
    tid = tp["id"]
    f_rad, f_lorentz = calc.calculate_forces(tp, E_totals[tid], B_totals[tid])
    results[tid] = {
        "E_total": E_totals[tid],
        "B_total": B_totals[tid],
        "f_rad": f_rad,
        "f_lorentz": f_lorentz
    }

with open("final_radiation_results.pkl", "wb") as f:
    pickle.dump(results, f)

print("Final results saved to final_radiation_results.pkl")
