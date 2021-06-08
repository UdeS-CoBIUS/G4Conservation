#!/bin/bash

#SBATCH --job-name=getFigures
#SBATCH --mem=32G
#SBATCH --account=def-jpviroid
#SBATCH --time=06:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user anais.vannutelli@usherbrooke.ca

python /home/vana2406/scratch/G4Conservation/scripts/getMainDensities.py -sp $1

python /home/vana2406/scratch/G4Conservation/scripts/getDataFig.py -sp $1
