#!/bin/bash
#SBATCH --job-name=KHmpi
#SBATCH --partition=M-large-t
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=1G
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.out
#SBATCH --time=00:15:00
#SBATCH --hint=nomultithread
#SBATCH --reservation=beginner

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
cd ${SLURM_SUBMIT_DIR}
module list

srun ./kh_mpi.x
