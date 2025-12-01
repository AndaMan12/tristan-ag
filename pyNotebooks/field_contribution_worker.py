import argparse
import pickle
import numpy as np
from graphet import Data
from graphet.plugins import TristanV2
from radiation_utils import RadiationReactionCalculator


def main():
    parser = argparse.ArgumentParser(description="Compute field contributions for a test particle over frame range.")
    parser.add_argument("--out_dir", required=True, help="Path to output directory.")
    parser.add_argument("--input_file", required=True, help="Input configuration file.")
    parser.add_argument("--frame_obs", type=int, required=True, help="Observation frame number.")
    parser.add_argument("--frame_start", type=int, required=True, help="First past frame to compute.")
    parser.add_argument("--frame_end", type=int, required=True, help="Last past frame to compute (inclusive).")
    parser.add_argument("--particle_file", required=True, help="Path to .pkl file containing test particle dict.")
    parser.add_argument("--output_file", required=True, help="Path to save .pkl output (field results).")

    args = parser.parse_args()

    # Load test particle
    with open(args.particle_file, "rb") as f:
        test_particle = pickle.load(f)

    # Load simulation parameters (from the same file as in wrapper)
    from radiation_utils import sim_params
    calc = RadiationReactionCalculator(**sim_params)

    field_results = []

    for past_frame in range(args.frame_start, args.frame_end + 1):
        E, B = calc.calculate_field_contribution_from_frame(test_particle, args.frame_obs, past_frame)
        field_results.append({
            "past_frame": past_frame,
            "E": E,
            "B": B
        })

    with open(args.output_file, "wb") as f:
        pickle.dump(field_results, f)

    print(f"Saved field results for frames {args.frame_start}–{args.frame_end} to {args.output_file}")


if __name__ == "__main__":
    main()

