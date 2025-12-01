#!/bin/bash

#SBATCH --job-name="Alfven_run_sweep"
#SBATCH --mail-type=END
#SBATCH --mail-user=anindyaguria@iisc.ac.in
#SBATCH -p debug
#SBATCH -t 20-24:00:00  # dd-hh:mm:ss
#SBATCH -n 220
#SBATCH --output=%x-%j.log

# Load necessary modules
module load zlib/
module load szip/
module load mpich/
module load hdf5/parallel


# Define parameter values
PSI_VALUES=(0.5 )  # Values for psi
MULT1_VALUES=(3 )                  # Values for multiplicities
MULT2_VALUES=(0.1)

# Path to the master input file
MASTER_INPUT=../inputs/inputAG.2d_EM_wave_source
PYTHON_SCRIPT=../pyNotebooks/COMP_calc_src.py  # Ensure this script is in the correct directory

for psi in "${PSI_VALUES[@]}"; do
  for mul1 in "${MULT1_VALUES[@]}"; do
    for mul2 in "${MULT2_VALUES[@]}"; do
    # Create a temporary copy of the input file
     NEW_INPUT="temp_input_psi${psi}_mul1${mul1}_mul2${mul2}.in"

    # Modify psi and multiplicity values in the input file
    sed -e "s/^\(\s*psi\s*=\s*\).*/\1$psi/" \
        -e "s/^\(\s*multiplicity_1\s*=\s*\).*/\1$mul1/" \
        -e "s/^\(\s*multiplicity_2\s*=\s*\).*/\1$mul2/" \
        "$MASTER_INPUT" > "$NEW_INPUT"

    # Run the Python script to compute the correct c_omp value
    COMP_TRUE=$(python3 "$PYTHON_SCRIPT" -i "$NEW_INPUT" | tail -1)

    echo $COMP_TRUE

    # Update the c_omp value in the input file
    sed -i "s/^\(\s*c_omp\s*=\s*\).*/\1$COMP_TRUE/" "$NEW_INPUT"

    # Run the simulation
    RUN_TS=$(date +%Y%m%d_%H%M%S)
    LOGFILE="out_TSI_EM_psi${psi}_mul1${mul1}_mul2${mul2}_${RUN_TS}.log"

    srun ../bin/tristan-mp2d -i "$NEW_INPUT" > "$LOGFILE"

    # Move the results to the vault
    OUTPUT_TAG="../vault/output_psi${psi}_mul1_${mul1}_mul2_${mul2}_src"
    mv "$NEW_INPUT" output
    mv "$LOGFILE" output
    mv output "$OUTPUT_TAG"
    done
  done
done

