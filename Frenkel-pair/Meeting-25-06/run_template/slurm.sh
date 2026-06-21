#!/bin/bash
#SBATCH -J jobFP_XXX.sh
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=30G
#SBATCH --output=/dev/null

echo $(hostname)
kmc base.py XXX

