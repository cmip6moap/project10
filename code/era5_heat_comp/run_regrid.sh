#!/bin/bash
#SBATCH -p short-serial
#SBATCH -o job_output/%j.out
#SBATCH -e job_output/%j.err
#SBATCH -t 24:00:00
#SBATCH --mem=32000 

conda activate heatstress
./regrid_models.py

## Usage:
## make sure you are in the correct directory: from the repo base:
## $ cd code

## then
## $ sbatch run_regrid.sh 
