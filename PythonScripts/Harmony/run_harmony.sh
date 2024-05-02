#!/bin/bash
#SBATCH --account=def-dmberman
#SBATCH --job-name=scanorama
#SBATCH --qos=privileged
#SBATCH --nodes=1                # number of Nodes
#SBATCH --tasks-per-node=1        # number of MPI processes per node
#SBATCH --mem 10G
#SBATCH --time 00:30:00
#SBATCH --output=/home/garvena/projects/def-dmberman/garvena/singlecellseq/SlurmOut/scanorama.%J.out
#SBATCH --error=/home/garvena/projects/def-dmberman/garvena/singlecellseq/SlurmOut/scanorama.%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=#

module load StdEnv/2023

module load python/3.11
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip

pip install --no-index harmony_pytorch
pip install --no-index scanpy

python /home/garvena/projects/def-dmberman/garvena/singlecellseq/PythonScripts/harmony/slurm_harmony.py