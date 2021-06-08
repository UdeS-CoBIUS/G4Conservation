#!/bin/bash

#SBATCH --job-name=G4RNAScreenerSubmit
#SBATCH --mem=4G
#SBATCH --account=def-jpviroid
#SBATCH --time=24:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user anais.vannutelli@usherbrooke.ca

source /home/vana2406/software/vpythonG4RNAScreener/bin/activate

list_sp=$(ls data/)
list_spMissing=()

for sp in $list_sp;
do
  nb_csv=$(ls /home/vana2406/scratch/G4Conservation/data/$sp/CSVFile/ | wc -l)
  nb_split=$(ls /home/vana2406/scratch/G4Conservation/data/$sp/SplitFile/ | wc -l)
  if [[ $nb_csv -lt $nb_split ]];
  then
    list_spMissing+="$sp "
   fi;
done

for sp in $list_spMissing;
do
  list_csv=$(ls /home/vana2406/scratch/G4Conservation/data/$sp/CSVFile/ | cut -d '.' -f 1)
  list_split=$(ls /home/vana2406/scratch/G4Conservation/data/$sp/SplitFile/ |  cut -d '.' -f 1)
  diff <(echo "$list_split") <(echo "$list_csv") > test
  for file in $(cat test | grep 'Sequences' | grep '<' | cut -d ' ' -f 2| sort -u);
  do
    jobs_nb=$(squeue -u vana2406 | wc -l)
    if [ $jobs_nb -le 990 ];
    then
      echo "/home/vana2406/software/g4rna_screener/screen.py data/$sp/SplitFile/$file.fa \
      -a /home/vana2406/software/g4rna_screener/G4RNA_2016-11-07.pkl -w 60 -s 10 \
      -c description cGcC G4H G4NN sequence start end -e \
      >  data/$sp/CSVFile/$file.csv"
      sbatch scripts/G4RNAscreener.sh data/$sp/SplitFile data/$sp/CSVFile $file
      sleep 0.1
    else
      sleep 42
    fi
  done
done
