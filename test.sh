#!/bin/bash 

reference_file=""
control_file=""
test_file=""

helpFunction()
{
   echo ""
   echo "$0 usage: [-r REFERENCE] [-c CONTROL] [-t TEST]"
   echo -e "\t-r reference text file with names"
   echo -e "\t-c Location of the control fastq file"
   echo -e "\t-t Location of the test fastq file"
   exit 1 # Exit script after printing help
}

while getopts "r:c:t:" opt
do
   case "$opt" in
      r ) reference_file="$OPTARG" ;;
      c ) control_file="$OPTARG" ;;
      t ) test_file="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
# if [ -r "$reference_file" ] || [ -c "$control_file" ] || [ -t "$test_file" ]
# then
#    echo "Must specify all required parameters.";
#    helpFunction
# fi

# Begin script in case all parameters are correct
#echo $reference_file
#echo $control_file
#echo $test_file

#cat transcripts.txt | xargs -P 4 -n 2 ./dana.sh $test_file $control_file
cat $reference_file | xargs -P 8 -n 1 ./adeno-m6A-KO1-P1.sh $test_file $control_file

#mkdir $1
