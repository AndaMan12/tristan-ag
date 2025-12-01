#!/bin/bash


#SBATCH --job-name="TSI_run"
#SBATCH --mail-type=END         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=anindyaguria@iisc.ac.in    # Where to send mail.  Set this to your email address
#SBATCH -p debug
#SBATCH -t 20-24:00:00  #dd-hh:mm:ss
#SBATCH -n 25
#SBATCH --output=%x-%j.log


module load zlib/
module load szip/
module load mpich/
module load hdf5/parallel


srun ../bin/tristan-mp1dTSI1D -i ../inputs/inputAG.1dtwostream_thermal > outTSI.log
