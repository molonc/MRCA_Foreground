#!/bin/bash
#SBATCH --job-name=fg_head
#SBATCH --partition=upgrade
#SBATCH --account=aparicio_lab
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

echo "Running head node fg job"
nextflow main.nf -resume
