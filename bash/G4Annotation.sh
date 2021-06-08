#!/bin/bash

#SBATCH --job-name=G4Annotation
#SBATCH --mem=32G
#SBATCH --account=def-jpviroid
#SBATCH --time=016:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user anais.vannutelli@usherbrooke.ca

python /home/vana2406/scratch/G4Conservation/scripts/G4Annotation.py -sp $1 -W $2
