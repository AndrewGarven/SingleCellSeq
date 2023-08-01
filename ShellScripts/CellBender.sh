#!/bin/bash
#SBATCH --account=def-dmberman
#SBATCH --job-name=cellBender
#SBATCH --nodes=1                # number of Nodes
#SBATCH --gpus-per-node=a100:2
#SBATCH --ntasks-per-node=32       
#SBATCH --mem=127000M
#SBATCH --time 08:00:00
#SBATCH --output=/home/garvena/projects/def-dmberman/garvena/singlecellseq/SlurmOut/cellBender.%J.out
#SBATCH --error=/home/garvena/projects/def-dmberman/garvena/singlecellseq/SlurmOut/cellBender.%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=#

module load apptainer 

# Set the root directory where the subdirectories with raw_feature_bc_matrix.h5 files are located
input_directory="/home/garvena/projects/def-dmberman/garvena/singlecellseq/Data/CellRanger/Lai/"
output_directory="/home/garvena/projects/def-dmberman/garvena/singlecellseq/Data/CellBender/Lai/"

# Iterate through each subdirectory and process the raw_feature_bc_matrix.h5 file
for subdir in "$input_directory"/*; do
  # Check if the raw_feature_bc_matrix.h5 file exists in the current subdirectory
  if [ -f "$subdir/outs/raw_feature_bc_matrix.h5" ]; then
    sample=$(basename "$subdir")
    
    # Process the raw_feature_bc_matrix.h5 file in the current subdirectory
    apptainer exec -C -W ${SLURM_TMPDIR} --nv \
      --bind "$input_directory":"/mnt/input" \
      --bind "$output_directory":"/mnt/output" \
      cellbender.sif cellbender remove-background \
      --input "/mnt/input/$(basename "$subdir")/outs/raw_feature_bc_matrix.h5" \
      --output "/mnt/output/$sample.h5" \
      --cuda \
      --expected-cells 8000 \
      --total-droplets-included 20000 \
      --fpr 0.01 \
      --epochs 150
  fi
done
