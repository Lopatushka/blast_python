#!/bin/bash
#SBATCH --job-name=blast
#SBATCH --nodes=1
#SBATCH --nodelist=node08
#SBATCH --ntasks-per-node=32
#SBATCH --partition=debug
#SBATCH --error=./server/stderr.log
#SBATCH --output=./server/stdout.log

# Get arguments
DIR=$1
DB=$2

# Run program
python main.py -d $DIR -db $DB -nt 32  

# Move output to directory with data
mv ./server/stdout.log $DIR
mv ./server/stderr.log $DIR 
