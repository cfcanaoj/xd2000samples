#!/bin/bash
#SBATCH --job-name=KH
#SBATCH --partition=M-debug
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.out
#SBATCH --hint=nomultithread

cd ${SLURM_SUBMIT_DIR}
module list

echo $SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# gnu in CPE
export OMP_PROC_BIND=True
export OMP_PLACES=cores

srun -c${SLURM_CPUS_PER_TASK} ./kh_hybrid.x
