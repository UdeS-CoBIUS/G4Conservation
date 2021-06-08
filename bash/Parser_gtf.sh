#!/bin/bash

#SBATCH --job-name=ParserGTF
#SBATCH --mem=32G
#SBATCH --account=def-jpviroid
#SBATCH --time=0-16:00:00

source /home/vana2406/software/vpythonG4RNAScreener/bin/activate

python /home/vana2406/scratch/G4Conservation/scripts/Parser_gtf.py -sp $1 -chr $2
