#!/bin/bash

#SBATCH --job-name="BTSI_run_sweep"
#SBATCH --mail-type=END
#SBATCH --mail-user=anindyaguria@iisc.ac.in
#SBATCH -p debug
#SBATCH -t 20-24:00:00  # dd-hh:mm:ss
#SBATCH -n 10
#SBATCH --output=%x-%j.log

# Load necessary modules
module load zlib/
module load szip/
module load mpich/
module load hdf5/parallel


# Define parameter values
SHIFT_GAMMA_VALUES=( 1.020621 )  # Values for shift gamma (beta 0.5 and 0.2)

TT_VALUES=(0 ) # Values for temperatures

# Path to the master input file
MASTER_INPUT=../inputs/inputAG.2dtwostream_thermal_Bfield

for TT in "${TT_VALUES[@]}"; do
  for gamma in "${SHIFT_GAMMA_VALUES[@]}"; do    
  # Create a temporary copy of the input file
  NEW_INPUT="temp_input_gamma${gamma}_TT${TT}_noB.in"

  # Modify psi and multiplicity values in the input file
  sed -e "s/^\(\s*temperature\s*=\s*\).*/\1$TT/" \
      -e "s/^\(\s*shift_gamma\s*=\s*\).*/\1$gamma/" \
      "$MASTER_INPUT" > "$NEW_INPUT"

  
  # Run the simulation
  RUN_TS=$(date +%Y%m%d_%H%M%S)
  LOGFILE="out_B_TSI_${gamma}_TT${TT}_${RUN_TS}_noB.log"

  srun ../bin/tristan-mp2dBTSI_2D -i "$NEW_INPUT" > "$LOGFILE"

  # Move the results to the vault
  OUTPUT_TAG="../vault/output_gamma${gamma}_TT_${TT}_noB_TSI"
  mv "$NEW_INPUT" output
  mv "$LOGFILE" output
  mv output "$OUTPUT_TAG"
  done  
done

