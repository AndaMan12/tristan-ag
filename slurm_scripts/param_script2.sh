#!/bin/bash

#SBATCH --job-name="Alfven_run_sweep"
#SBATCH --mail-type=all
#SBATCH --mail-user=anindyaguria@iisc.ac.in
#SBATCH -p normal
#SBATCH -t 48:00:00  # dd-hh:mm:ss
#SBATCH -N 4
#SBATCH -n 80
#SBATCH --output=%x-%j.log

# Load necessary modules
# module load zlib/
# module load szip/
# module load mpich/
# module load hdf5/parallel
module restore tristan_modules

# Define parameter values
PSI_VALUES=(0.5 )        # Values for psi
MULT1_VALUES=(3 )        # Values for multiplicities
MULT2_VALUES=(0.1 )
TT_VALUES=(1.9e-2 )      # Values for temperatures

# Path to the master input file
MASTER_INPUT=/work/10446/anindya_12/ls6/tristan_mp_v2/inputs/inputAG.2d_EM_wave_embed_mw
PYTHON_SCRIPT=/work/10446/anindya_12/ls6/tristan_mp_v2/pyNotebooks/COMP_calc.py  # Ensure this script is in the correct directory

for psi in "${PSI_VALUES[@]}"; do
  # Compute cos(psi) in radians for movwingam
  MOVWINGAM_TRUE=$(python3 -c "import math; print(math.cos($psi))")
  
  for mul1 in "${MULT1_VALUES[@]}"; do
    for mul2 in "${MULT2_VALUES[@]}"; do
      for TT in "${TT_VALUES[@]}"; do

        # Create a temporary copy of the input file
        NEW_INPUT="temp_input_psi${psi}_mul1${mul1}_mul2${mul2}_TT${TT}.in"

        # Modify psi, movwingam and other parameters in the input file
        sed -e "s/^\(\s*psi\s*=\s*\).*/\1$psi/" \
            -e "s/^\(\s*movwingam\s*=\s*\).*/\1$MOVWINGAM_TRUE/" \
            -e "s/^\(\s*multiplicity_1\s*=\s*\).*/\1$mul1/" \
            -e "s/^\(\s*multiplicity_2\s*=\s*\).*/\1$mul2/" \
            -e "s/^\(\s*temperature\s*=\s*\).*/\1$TT/" \
            "$MASTER_INPUT" > "$NEW_INPUT"

        # Run the Python script to compute the correct c_omp value
        COMP_TRUE=$(python3 "$PYTHON_SCRIPT" -i "$NEW_INPUT" | tail -1)
        echo "$COMP_TRUE"

        # Update the c_omp value in the input file
        sed -i "s/^\(\s*c_omp\s*=\s*\).*/\1$COMP_TRUE/" "$NEW_INPUT"

        # Run the simulation
        RUN_TS=$(date +%Y%m%d_%H%M%S)
        LOGFILE="out_TSI_EM_psi${psi}_mul1${mul1}_mul2${mul2}_TT${TT}_${RUN_TS}.log"

        # srun /work/10446/anindya_12/ls6/tristan_mp_v2/bin/tristan-mp2d -i "$NEW_INPUT" > "$LOGFILE"

        # Move the results to the vault
        # OUTPUT_TAG="/scratch/10446/anindya_12/ls6/tristan_mp_v2/vault/output_psi${psi}_mul1_${mul1}_mul2_${mul2}_TT_${TT}_mov_win"
        # mv "$NEW_INPUT" output
        # mv "$LOGFILE" output
        # mv output "$OUTPUT_TAG"

      done
    done
  done
done
