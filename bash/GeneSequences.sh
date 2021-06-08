#!/bin/bash

#SBATCH --job-name=G4RNAScreener
#SBATCH --mem=16G
#SBATCH --account=def-jpviroid
#SBATCH --time=06:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user anais.vannutelli@usherbrooke.ca

conda activate smake

python3 scripts/GeneSequences.py -sp $1
