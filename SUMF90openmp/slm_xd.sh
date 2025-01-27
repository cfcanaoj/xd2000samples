#!/bin/bash
#SBATCH --job-name=SUM
#SBATCH --partition=M-test-cfca
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH -o slmlog%J.out
#SBATCH --hint=nomultithread


cd ${SLURM_SUBMIT_DIR}

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export SLURM_CPU_BIND=ldoms
#export SLURM_CPU_BIND=rank_ldom

date >& log.$SLURM_JOB_ID
time srun ./sum.x >> log.$SLURM_JOB_ID
#time srun vtune -collect hotspots -r=vtout$SLURM_JOB_ID  ./sum.x >> log.$SLURM_JOB_ID 
#time srun vtune -collect threading -r=vtout$SLURM_JOB_ID  ./sum.x >> log.$SLURM_JOB_ID 

date >> log.$SLURM_JOB_ID
