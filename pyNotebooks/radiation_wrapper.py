import os
import pickle
import numpy as np
import subprocess
from tqdm import tqdm
import multiprocessing
from radiation_utils import RadiationReactionCalculator, sim_params
import tristanVis.isolde as isolde

def run_field_worker(calc_params, test_particle, obs_frame, frame_chunks):
    """
    Launch multiple subprocesses of field_contribution_worker.py in parallel.
    (Already parallelized across frame chunks.)
    """
    particle_id = test_particle["id"]
    particle_file = f"test_particle_{particle_id}.pkl"

    with open(particle_file, "wb") as f:
        pickle.dump(test_particle, f)

    # Update global sim_params file so that worker sees correct values
    with open("radiation_utils.py", "r") as f:
        lines = f.readlines()
    with open("radiation_utils.py", "w") as f:
        for line in lines:
            if line.strip().startswith("sim_params ="):
                f.write("sim_params = " + repr(calc_params) + "\n")
            else:
                f.write(line)

    processes = []
    output_files = []

    for i, (start, end) in enumerate(frame_chunks):
        output_file = f"field_chunk_{particle_id}_{i}.pkl"
        output_files.append(output_file)
        cmd = [
            "python3", "field_contribution_worker.py",
            "--out_dir", calc_params["out_dir"],
            "--input_file", calc_params["input_file_name"],
            "--frame_obs", str(obs_frame),
            "--frame_start", str(start),
            "--frame_end", str(end),
            "--particle_file", particle_file,
            "--output_file", output_file
        ]
        p = subprocess.Popen(cmd)
        processes.append(p)

    # Wait for all to finish
    for p in tqdm(processes, desc=f"Running field workers (particle {particle_id})"):
        p.wait()

    # Collect results
    all_results = []
    for f in output_files:
        if os.path.exists(f):
            with open(f, "rb") as fp:
                all_results.extend(pickle.load(fp))
            os.remove(f)

    os.remove(particle_file)
    all_results.sort(key=lambda r: r["past_frame"], reverse=True)
    return all_results


def process_single_particle(calc_params, frame, max_frames_back, num_chunks, test_particle):
    """
    Process radiation reaction analysis for a single test particle.
    """
    calc = RadiationReactionCalculator(**calc_params)

    past_frames = list(range(frame - max_frames_back, frame))
    chunks = np.array_split(past_frames, num_chunks)
    frame_chunks = [(int(ch[0]), int(ch[-1])) for ch in chunks if len(ch) > 0]

    field_data = run_field_worker(calc_params, test_particle, frame, frame_chunks)

    # Combine fields cumulatively and compute forces
    E_total = np.zeros(3)
    B_total = np.zeros(3)
    results = {
        "test_particle": test_particle,
        "frames_back": [],
        "max_distances": [],
        "f_rad": [],
        "f_lorentz": []
    }

    for entry in field_data:
        E_total += entry["E"]
        B_total += entry["B"]
        f_rad, f_lorentz = calc.calculate_forces(test_particle, E_total, B_total)
        results["frames_back"].append(entry["past_frame"])
        results["max_distances"].append(calc.CC * calc.interval * (frame - entry["past_frame"]))
        results["f_rad"].append(f_rad)
        results["f_lorentz"].append(f_lorentz)

    return test_particle["id"], results


def run_radiation_analysis(calc_params, frame, max_frames_back, num_chunks=12, parallel=True):
    """
    Run radiation analysis for all test particles (not just one).
    Optionally parallelize across particles.
    """
    calc = RadiationReactionCalculator(**calc_params)

    # Get all test particles in the region
    test_particles = calc.get_test_particle(frame, region={"x_range": calc_params["bounds_x"],"y_range":calc_params["bounds_y"]})

    if not isinstance(test_particles, list):
        test_particles = [test_particles]  # backward compatibility

    print(f"Selected {len(test_particles)} test particles for analysis.")

    all_results = {}

    if parallel and len(test_particles) > 1:
        # Outer-level parallelization over particles
        with multiprocessing.Pool(processes=min(len(test_particles), os.cpu_count())) as pool:
            tasks = [
                pool.apply_async(
                    process_single_particle,
                    (calc_params, frame, max_frames_back, num_chunks, p)
                )
                for p in test_particles
            ]
            for task in tqdm(tasks, desc="Processing all particles"):
                pid, res = task.get()
                all_results[pid] = res
    else:
        # Serial fallback
        for p in tqdm(test_particles, desc="Processing particles"):
            pid, res = process_single_particle(calc_params, frame,  max_frames_back, num_chunks, p)
            all_results[pid] = res

    # Save results
    with open("final_radiation_results.pkl", "wb") as f:
        pickle.dump(all_results, f)

    print(f"Analysis complete for {len(test_particles)} particles.")
    print("Results saved to final_radiation_results.pkl")

    return all_results


if __name__ == "__main__":
    import tristanVis.isolde as isolde

    out_dir = "/scratch/10446/anindya_12/tristan_mp_v2/vault/output_psi0.5_mul1_3_mul2_0.1_TT_1.9e-4_rad_drag_prtl_maxxed/"
    input_file_name = out_dir + "temp_input_psi0.5_mul13_mul20.1_TT1.9e-4.in"

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

    unit_ch = CC**2 / (ppc0 * COMP**2)
    weight_fac = 2.0473518E-03
    ds = wA_wp * 2 * np.pi / mode
    omegap0 = CC / ds
    init_x_boundary = int(5 * np.pi / mode) + 1
    fin_x_boundary = init_x_boundary + ramp_width * (2 * np.pi / mode)

    bounds_x = [10220, 10280]
    bounds_y = [1020, 1050]
    frame = 2400 * 5    
    max_frames_back = 150 * 5
    num_chunks = 20

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

    run_radiation_analysis(calc_params, frame,  max_frames_back, num_chunks, parallel=False)





# import os
# import pickle
# import numpy as np
# import subprocess
# from tqdm import tqdm
# import multiprocessing
# from radiation_utils import RadiationReactionCalculator, sim_params


# def run_field_worker(calc_params, test_particle, obs_frame, frame_chunks):
#     """
#     Launch multiple subprocesses of field_contribution_worker.py in parallel.
#     """
#     with open("test_particle.pkl", "wb") as f:
#         pickle.dump(test_particle, f)

#     # Update global sim_params file so that worker sees correct values
#     with open("radiation_utils.py", "r") as f:
#         lines = f.readlines()
#     with open("radiation_utils.py", "w") as f:
#         for line in lines:
#             if line.strip().startswith("sim_params ="):
#                 f.write("sim_params = " + repr(calc_params) + "\n")
#                 print("sim_params = " + repr(calc_params) + "\n")
#             else:
#                 f.write(line)

#     processes = []
#     output_files = []

#     for i, (start, end) in enumerate(frame_chunks):
#         output_file = f"field_chunk_{i}.pkl"
#         output_files.append(output_file)
#         cmd = [
#             "python3", "field_contribution_worker.py",
#             "--out_dir", calc_params["out_dir"],
#             "--input_file", calc_params["input_file_name"],
#             "--frame_obs", str(obs_frame),
#             "--frame_start", str(start),
#             "--frame_end", str(end),
#             "--particle_file", "test_particle.pkl",
#             "--output_file", output_file
#         ]
#         p = subprocess.Popen(cmd)
#         processes.append(p)

#     # Wait for all to finish
#     for p in tqdm(processes, desc="Running field workers"):
#         p.wait()
# #     tqdm(processes, desc="Running field workers")

#     # Collect results
#     all_results = []
#     for f in output_files:
#         if os.path.exists(f):
#             with open(f, "rb") as fp:
#                 all_results.extend(pickle.load(fp))
#             os.remove(f)

#     all_results.sort(key=lambda r: r["past_frame"], reverse=True)
#     return all_results


# def run_radiation_analysis(calc_params, frame, species, max_frames_back, num_chunks=8):
#     calc = RadiationReactionCalculator(**calc_params)

#     # Pick test particle
#     test_particle = calc.get_test_particle(frame, species, region={
#         "x_range": calc_params["bounds_x"],
#         "y_range": calc_params["bounds_y"]
#     })

#     print(f"Selected test particle ID={test_particle['id']}")

#     # Divide past frames into chunks
#     past_frames = list(range(frame - max_frames_back, frame))
#     chunks = np.array_split(past_frames, num_chunks)
#     frame_chunks = [(int(ch[0]), int(ch[-1])) for ch in chunks if len(ch) > 0]

#     field_data = run_field_worker(calc_params, test_particle, frame, frame_chunks)

#     # Combine fields cumulatively and compute forces
#     E_total = np.zeros(3)
#     B_total = np.zeros(3)
#     results = {
#         "test_particle": test_particle,
#         "frames_back": [],
#         "max_distances": [],
#         "f_rad": [],
#         "f_lorentz": []
#     }

#     for entry in tqdm(field_data, desc="Combining results"):
#         E_total += entry["E"]
#         B_total += entry["B"]
#         f_rad, f_lorentz = calc.calculate_forces(test_particle, E_total, B_total)
#         results["frames_back"].append(entry["past_frame"])
#         results["max_distances"].append(calc.CC * calc.interval * (frame - entry["past_frame"]))
#         results["f_rad"].append(f_rad)
#         results["f_lorentz"].append(f_lorentz)

#     with open("final_radiation_results.pkl", "wb") as f:
#         pickle.dump(results, f)

#     print("Analysis complete. Results saved to final_radiation_results.pkl")
#     return results


# if __name__ == "__main__":
#     # Example parameter setup — replace with your actual values
#     out_dir = "/scratch/10446/anindya_12/tristan_mp_v2/vault/output_psi0.5_mul1_3_mul2_0.1_TT_1.9e-4_rad_drag_prtl_maxxed/"
#     input_file_name = out_dir + "temp_input_psi0.5_mul13_mul20.1_TT1.9e-4.in"
    
#     # Parse input parameters to extract values
#     import tristanVis.isolde as isolde
#     input_params = isolde.parseInput(input_file_name)
    
#     interval = input_params["output"]["interval"]
#     lst_time = input_params["time"]["last"]
#     grid_x = int(input_params["grid"]["mx0"])
#     grid_y = int(input_params["grid"]["my0"])
#     CC = input_params["algorithm"]["c"]
#     COMP = input_params["plasma"]["c_omp"]
#     ppc0 = input_params["plasma"]["ppc0"]
#     stride = input_params["output"]["stride"]
#     wA_wp = input_params["problem"]["wA_wp"]
#     mode = input_params["problem"]["mode"]
#     ramp_width = input_params["problem"]["ramp_width"]
    
#     # Derived parameters
#     unit_ch = CC**2 / (ppc0 * COMP**2)
#     weight_fac = 2.0473518E-03  # From your notebook
#     ds = wA_wp * 2 * np.pi / mode
#     omegap0 = CC / ds
#     init_x_boundary = int(5 * np.pi / mode) + 1
#     fin_x_boundary = init_x_boundary + ramp_width * (2 * np.pi / mode)
    
#     # Set analysis parameters (modify these as needed)
#     bounds_x = [10220, 10280]  # e.g., [init_x_boundary, fin_x_boundary]
#     bounds_y = [1020, 1050]  # e.g., [0, grid_y]
#     frame = 2400*5  # e.g., 100 or int(lst_time // interval) // 2
#     species = 2  # 1 for positrons, 2 for electrons
#     max_frames_back = 150*5  # Maximum retarded time frames
#     calc_params = {
#         "out_dir": out_dir,
#         "input_file_name": input_file_name,
#         "bounds_x": bounds_x,
#         "bounds_y": bounds_y,
#         "interval": interval,
#         "omegap0": omegap0,
#         "wA_wp": wA_wp,
#         "weight_fac": weight_fac,
#         "unit_ch": unit_ch,
#         "CC": CC,
#         "stride": stride
#     }

#     frame = 2400*5
#     species = 2
#     max_frames_back = 150*5
#     num_chunks = 20

#     run_radiation_analysis(calc_params, frame, species, max_frames_back, num_chunks)

