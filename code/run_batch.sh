#!/bin/bash
#SBATCH -p short-serial
#SBATCH -o job_output/%j.out
#SBATCH -e job_output/%j.err
#SBATCH -t 24:00:00
#SBATCH --mem=64000  # 64M needed for HadGEM3-MM, else 32M usually sufficient
#SBATCH --array=0-45

conda activate heatstress
./run_utci.py ${SLURM_ARRAY_TASK_ID}

## Usage:
## make sure you are in the correct directory: from the repo base:
## $ cd jasmin-scripts

## then
## $ sbatch run_batch.sh 
