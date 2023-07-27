#!/bin/bash
#SBATCH --account=def-dmberman
#SBATCH -J test_cellranger
#SBATCH --export=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --signal=2
#SBATCH --no-requeue
### Alternatively: --ntasks=1 --cpus-per-task={NUM_THREADS}
###   Consult with your cluster administrators to find the combination that
###   works best for single-node, multi-threaded applications on your system.
#SBATCH --mem=250G
#SBATCH --time 02:00:00
#SBATCH --output=/home/garvena/projects/def-dmberman/garvena/singlecellseq/SlurmOut/test_cellranger.%J.out
#SBATCH --error=/home/garvena/projects/def-dmberman/garvena/singlecellseq/SlurmOut/test_cellranger.%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=#
# Replace "your_directory_path" with the actual directory path you want to search

cd /home/garvena/projects/def-dmberman/garvena/singlecellseq/Data/CellRanger/Lai

d = '/home/garvena/projects/def-dmberman/garvena/singlecellseq/Data/InputData/Lai/SRR9897622/'

cellranger count --id=$(basename "$d") \
                 --transcriptome=/home/garvena/projects/def-dmberman/garvena/singlecellseq/Data/genome-ref/refdata-gex-GRCh38-2020-A \
                 --fastqs=$d \
                 --chemistry=SC3Pv2
