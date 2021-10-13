#!/bin/bash

#SBATCH --job-name=G4RNAScreener
#SBATCH --account=def-jpviroid
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user anais.vannutelli@usherbrooke.ca

startTime=`date +%s`
list_file=$(ls $1)

for sp in $(ls ./reviewShuffle/);
do
  for repro in $(ls ./reviewShuffle/$sp/);
  do
    for file in $(ls ./reviewShuffle/$sp/$repro/SplitFile/ | cut -d '.' -f 1);
    do
      jobs_nb=$(squeue -u vana2406 | wc -l)
      if [ $jobs_nb -le 850 ]; then
          echo $file
          sbatch ~/scratch/G4Conservation/scripts/G4RNAscreener.sh reviewShuffle/$sp/$repro/SplitFile reviewShuffle/$sp/$repro/CSVFile $file
          sleep 0.1
      else
          sleep 42
      fi;
    done
  done
done
