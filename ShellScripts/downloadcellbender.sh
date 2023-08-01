#!/bin/bash
#SBATCH --account=def-dmberman
#SBATCH --job-name=cellBender
#SBATCH --nodes=1                # number of Nodes
#SBATCH --ntasks-per-node=4      
#SBATCH --mem=20G
#SBATCH --time 02:00:00
#SBATCH --output=/home/garvena/projects/def-dmberman/garvena/singlecellseq/SlurmOut/cellBender.%J.out
#SBATCH --error=/home/garvena/projects/def-dmberman/garvena/singlecellseq/SlurmOut/cellBender.%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=#

module load StdEnv/2020
module load apptainer

apptainer pull cellbender.sif docker://us.gcr.io/broad-dsde-methods/cellbender:latest
