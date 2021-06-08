#!/bin/bash

#BATCH --job-name=generateSequences
#SBATCH --mem=38G
#SBATCH --account=def-jpviroid
#SBATCH --time=06:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user anais.vannutelli@usherbrooke.ca

python /home/vana2406/scratch/G4Conservation/scripts/generateSequences.py -sp $1
