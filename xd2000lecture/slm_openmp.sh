#!/bin/bash
#SBATCH --job-name=KHomp
#SBATCH --partition=M-debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=1G
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.out
#SBATCH --time=00:15:00
#SBATCH --hint=nomultithread
#SBATCH --reservation=beginner

cd ${SLURM_SUBMIT_DIR}
module list

echo $SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

export OMP_PROC_BIND=true
export OMP_PLACES=cores                                                         

./kh_openmp.x
