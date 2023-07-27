#!/bin/bash

#SBATCH --account=def-dmberman
#SBATCH --job-name=fastqc
#SBATCH --qos=privileged
#SBATCH --nodes=12                # number of Nodes
#SBATCH --tasks-per-node=4        # number of MPI processes per node
#SBATCH --mem 16g
#SBATCH --time 24:00:00
#SBATCH --output=/home/garvena/projects/def-dmberman/garvena/singlecellseq/SlurmOut/fastqc.%J.out
#SBATCH --error=/home/garvena/projects/def-dmberman/garvena/singlecellseq/SlurmOut/fastqc.%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=# add your email address

module load fastqc

fastqc ~/projects/def-dmberman/garvena/singlecellseq/Data/InputData/Lai/*/*.fastq.gz \
       -o ~/projects/def-dmberman/garvena/singlecellseq/Data/fastQC -f fastq -t 12