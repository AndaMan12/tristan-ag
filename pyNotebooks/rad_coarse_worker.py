# coarse_worker.py
import pickle
from rad_rxn import RadiationReactionCalculator

def run_chunk(calc_params, test_particle, frame, frame_subset, output_file):
    calc = RadiationReactionCalculator(**calc_params)
    
    results = []
    for pf in frame_subset:
        E_new, B_new = calc.calculate_field_contribution_from_frame(test_particle, frame, pf)        
        results.append((frame, E_new, B_new))

    with open(output_file, "wb") as f:
        pickle.dump(results, f)
