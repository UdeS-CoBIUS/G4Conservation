#!/bin/bash

#SBATCH --job-name=fasta_enum
#SBATCH --mem=4G
#SBATCH --account=def-jpviroid
#SBATCH --time=00:30:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user anais.vannutelli@usherbrooke.ca

conda activate smake

python3 /home/vana2406/scratch/G4Conservation/scripts/fasta_enumerator.py --fasta $1 --ntLimit 5000000 \
        --large False --output $2
