#!/bin/bash

#SBATCH --job-name=G4RNAScreener
#SBATCH --mem=16G
#SBATCH --account=def-jpviroid
#SBATCH --time=02:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user anais.vannutelli@usherbrooke.ca

conda activate smake

python3 scripts/fasta_enumerator.py --fasta $1 --ntLimit 5000000 \
        --large False --output $2
