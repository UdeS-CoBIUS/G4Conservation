#!/bin/bash

#SBATCH --job-name=G4AnnGene
#SBATCH --mem=32G
#SBATCH --account=def-jpviroid
#SBATCH --time=016:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user anais.vannutelli@usherbrooke.ca

python /home/vana2406/scratch/G4Conservation/scripts/G4AnnotationGene.py -sp $1
