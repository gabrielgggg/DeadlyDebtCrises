#!/bin/bash
#SBATCH --account=ACCOUNT_HERE_IF_NEEDED
#SBATCH -p PARTITION_HERE
#SBATCH --nodes=14
#SBATCH --ntasks-per-node=28
#SBATCH --mem=0 
#SBATCH -t 0-04:00:00
#SBATCH -J baseline
#SBATCH -o out.txt
#SBATCH -e err.txt

set -xe

# adjust module load or compiler config, if needed
module load intel/2021.3.0

export MV2_ENABLE_AFFINITY=0
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

echo "$SLURM_JOB_NUM_NODES nodes with $SLURM_CPUS_ON_NODE processors each." > log.txt

cd $SLURM_SUBMIT_DIR
date >> log.txt
make clean &>> log.txt
make -j$SLURM_CPUS_ON_NODE &>> log.txt

date >> log.txt
srun -n $SLURM_JOB_NUM_NODES -c $SLURM_CPUS_ON_NODE ./bin/pdef &>> log.txt

date >> log.txt
echo "End." &>> ./log.txt

