#!/bin/bash

#SBATCH --job-name=G4RNAScreenerCentro
#SBATCH --account=def-jpviroid
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user anais.vannutelli@usherbrooke.ca

startTime=`date +%s`

for sp in $(ls ./reviewTRCentro/);
do
  for repro in $(seq 1 10);
  do
    for file in $(ls ./reviewTRCentro/$sp/Repro$repro/SplitFile/ | grep -v 'centro' | cut -d '.' -f 1);
    do
      jobs_nb=$(squeue -u vana2406 | wc -l)
      if [ $jobs_nb -le 850 ]; then
          echo $file $sp $repro
          sbatch ~/scratch/G4Conservation/scripts/G4RNAscreener.sh reviewTRCentro/$sp/Repro$repro/SplitFile reviewTRCentro/$sp/Repro$repro/CSVFile $file
          sleep 0.1
      else
          sleep 42
      fi;
    done
  done
done

echo 'Shuffle done, now WT'

for sp in $(ls ./reviewTRCentro/);
do
  for file in $(ls ./reviewTRCentro/$sp/SplitFile/ | grep -v 'centro' | cut -d '.' -f 1);
  do
    jobs_nb=$(squeue -u vana2406 | wc -l)
    if [ $jobs_nb -le 850 ]; then
        echo $file $sp $repro
        sbatch ~/scratch/G4Conservation/scripts/G4RNAscreener.sh reviewTRCentro/$sp/SplitFile reviewTRCentro/$sp/CSVFile $file
        sleep 0.1
    else
        sleep 42
    fi;
  done
done
