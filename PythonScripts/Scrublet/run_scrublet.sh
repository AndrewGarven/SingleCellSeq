#!/bin/bash
#SBATCH --account=def-dmberman
#SBATCH --job-name=scrublet
#SBATCH --qos=privileged
#SBATCH --nodes=1                # number of Nodes
#SBATCH --tasks-per-node=1        # number of MPI processes per node
#SBATCH --mem 50G
#SBATCH --time 10:00:00
#SBATCH --output=/home/garvena/projects/def-dmberman/garvena/singlecellseq/SlurmOut/scrublet.%J.out
#SBATCH --error=/home/garvena/projects/def-dmberman/garvena/singlecellseq/SlurmOut/scrublet.%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=#

module load StdEnv/2023

module load python/3.10
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip

pip install --no-index scanpy
pip install --no-index scikit-image

# Run your Python script
python /home/garvena/projects/def-dmberman/garvena/singlecellseq/PythonScripts/slurm_scrublet.py