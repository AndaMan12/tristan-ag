import multiprocessing
import numpy as np
import pickle
import os
from tqdm import tqdm

def distance_analysis_chunked(self, test_particle, frame, max_frames_back, num_chunks=5):
    frame_ranges = np.arange(1, max_frames_back + 1)
    past_frames = [frame - i for i in frame_ranges if frame - i >= 0]
    chunks = np.array_split(past_frames, num_chunks)

    calc_params = {
        "out_dir": self.out_dir,
        "input_file_name": self.input_file_name,
        "bounds_x": self.bounds_x,
        "bounds_y": self.bounds_y,
        "interval": self.interval,
        "omegap0": self.omegap0,
        "wA_wp": self.wA_wp,
        "weight_fac": self.weight_fac,
        "unit_ch": self.unit_ch,
        "CC": self.CC,
        "stride": self.stride,
    }

    jobs = []
    for i, subset in enumerate(chunks):
        if len(subset) == 0:
            continue
        output_file = f"chunk_{i}.pkl"
        args = (calc_params, test_particle, frame, subset.tolist(), output_file)
        p = multiprocessing.Process(target=run_chunk, args=args)
        p.start()
        jobs.append(p)

    for p in tqdm(jobs, desc="Waiting for chunk processes"):
        p.join()

    # Merge partial results
    all_results = []
    for i, subset in enumerate(chunks):
        output_file = f"chunk_{i}.pkl"
        if os.path.exists(output_file):
            with open(output_file, "rb") as f:
                all_results.extend(pickle.load(f))
            os.remove(output_file)

    # Sort by frame number
    all_results.sort(key=lambda r: r[0])
    frames_back, E_s, B_s = zip(*all_results)
    
    return {"frames_back": frames_back, "f_rad": f_rad, "f_lorentz": f_lorentz}
