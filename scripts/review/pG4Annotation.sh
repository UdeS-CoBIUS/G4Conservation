#!/bin/bash

#SBATCH --job-name=G4Annotation
#SBATCH --mem=32G
#SBATCH --account=def-jpviroid
#SBATCH --time=3:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user anais.vannutelli@usherbrooke.ca

python /home/vana2406/scratch/G4Conservation/scripts/review/pG4Annotation.py -sp $1 -r $2
