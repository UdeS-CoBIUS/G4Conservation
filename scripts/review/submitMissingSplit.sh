#!/bin/bash

#SBATCH --job-name=G4RNAScreenerSubmit
#SBATCH --account=def-jpviroid
#SBATCH --time=4:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user anais.vannutelli@usherbrooke.ca

source /home/vana2406/software/vpythonG4RNAScreener/bin/activate

list_sp=$(ls reviewShuffle/)
list_spMissing=()
list_reproMissing=()

for sp in $list_sp;
do
  for repro in $(ls reviewShuffle/$sp | grep -v 'pG4_Shuffle');
  do
    nb_csv=$(ls /home/vana2406/scratch/G4Conservation/reviewShuffle/$sp/$repro/CSVFile/ | wc -l)
    nb_split=$(ls /home/vana2406/scratch/G4Conservation/reviewShuffle/$sp/$repro/SplitFile/ | wc -l)
    if [[ $nb_csv -lt $nb_split ]];
    then
      list_spMissing+="$sp "
      list_reproMissing+="$repro "
     fi;
  done
done

list_spMissing=$(for i in ${list_spMissing[@]}; do echo $i; done | sort -u)
list_reproMissing=$(for i in ${list_reproMissing[@]}; do echo $i; done | sort -u)

for sp in $list_spMissing;
do
  for repro in $list_reproMissing;
  do
    list_csv=$(ls /home/vana2406/scratch/G4Conservation/reviewShuffle/$sp/$repro/CSVFile/ | cut -d '.' -f 1)
    list_split=$(ls /home/vana2406/scratch/G4Conservation/reviewShuffle/$sp/$repro/SplitFile/ |  cut -d '.' -f 1)
    diff <(echo "$list_split") <(echo "$list_csv") > test
    for file in $(cat test | grep 'Sequences' | grep '<' | cut -d ' ' -f 2| sort -u);
    do
      # echo $file $repro $sp
      jobs_nb=$(squeue -u vana2406 | wc -l)
      if [ $jobs_nb -le 990 ];
      then
        echo $sp $repro $file
        sbatch scripts/G4RNAscreener.sh reviewShuffle/$sp/$repro/SplitFile reviewShuffle/$sp/$repro/CSVFile $file
        sleep 0.1
      else
        sleep 42
      fi
    done
  done
done
