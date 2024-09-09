#!/bin/bash
#SBATCH --job-name=blastn
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --nodelist=node08
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --error=./server/stderr.log
#SBATCH --output=./server/stdout.log

# Read positional arguments
DIR=$1
DB=$2

# Print SLURM Environment variables
echo "Running blastn with the following parameters:"
echo "SLURM_JOB_NAME: $SLURM_JOB_NAME"
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo "SLURM_JOB_NUM_NODES: $SLURM_JOB_NUM_NODES"
echo "SLURM_JOB_NODELIST: $SLURM_JOB_NODELIST"
echo "SLURM_NTASKS_PER_NODE: $SLURM_NTASKS_PER_NODE"
echo "SLURM_CPUS_PER_TASK: $SLURM_CPUS_PER_TASK"
echo "SLURM_CPUS_ON_NODE: $SLURM_CPUS_ON_NODE"

TOTAL_NTASKS=$((SLURM_NTASKS_PER_NODE * SLURM_JOB_NUM_NODES))
echo "TOTAL_NTASKS: $TOTAL_NTASKS"

TOTAL_CPUS_ALLOCATED=$((TOTAL_NTASKS * SLURM_CPUS_PER_TASK))
echo "TOTAL_CPUS_ALLOCATED: $TOTAL_CPUS_ALLOCATED"
echo

# Automatically set N_THREADS
N_THREADS=$TOTAL_CPUS_ALLOCATED

# Print command-line variables
echo "Path to directory with .ab1 files: $DIR"
echo "Database: $DB"
echo "Number of threads: $N_THREADS"
echo

# Run program
python main.py -d $DIR -db $DB -nt $N_THREADS  

# Move output to directory with data
mv ./server/stdout.log $DIR
mv ./server/stderr.log $DIR 
