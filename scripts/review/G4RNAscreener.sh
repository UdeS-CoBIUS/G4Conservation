#!/bin/bash

#SBATCH --job-name=G4RNAScreener
#SBATCH --mem=32G
#SBATCH --array=0
#SBATCH --account=def-jpviroid
#SBATCH --time=02:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user anais.vannutelli@usherbrooke.ca

source /home/vana2406/software/vpythonG4RNAScreener/bin/activate

/home/vana2406/software/g4rna_screener/screen.py $1/$3.fa \
-a /home/vana2406/software/g4rna_screener/G4RNA_2016-11-07.pkl -w 60 -s 10 \
-c description cGcC G4H G4NN sequence start end -e \
> $2/$3.csv
