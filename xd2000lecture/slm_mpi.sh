#!/bin/bash
#SBATCH --job-name=KH
#SBATCH --partition=M-debug
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=30G
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.out
#SBATCH --hint=nomultithread

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
cd ${SLURM_SUBMIT_DIR}
module list

srun ./kh_mpi.x
