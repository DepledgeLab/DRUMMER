#!/bin/bash
set -e
set -u
set -o pipefail

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"


usage(){
echo "
DRUMMER (2020 - present)

Description:  Detection of RNA modifications through comparative profiling of basecall 
error rates in nanopore direct RNA sequencing data

https://github.com/DepledgeLab/DRUMMER

Usage: drummer.sh -r [FASTA] -u|-n [TARGETS] -c [CONTROL] -t [TEST] -o [OUTPUT] -m [RUNMODE] (OPTIONS)

Required flags:
-r 		fasta format reference genome (exome) or transcriptome (isoforms)

-u 		list of transcripts (isoform) to be examined (single column or XX-column format)
OR
-n		name of genome (exome) - must match fasta file header

-c 		sorted.bam file - control (RNA modification(s) present)
-t 		sorted.bam file - treatment (RNA modification(s) absent)
-o 		output directory
-m 		runmode (exome|isoform)

Optional flags:
-y		specify odds ratio requirement (default >= 1.5)
-z		specify adjusted p_value (G-test) requirement (default<= 0.05)

DRUMMER was written by Jonathan S. Abebe & Daniel P. Depledge.
If you encounter any problems with DRUMMER, please log them as an issue at the GitHub page https://github.com/DepledgeLab/DRUMMER
"
}

### Function for showing usage
if [[ ( $# == "--help") ||  $# == "-h" ||  $# == 0 ]] 
    then 
    usage
    exit 0
fi


### GETOPTS
#PTSTRING contains the option letters to be recognized; if a letter is followed by a colon, the option is expected to have an argument

while getopts ":r:n:c:t:o:m:u:y:z:" opt; do
   case $opt in
      r) reference_file="$OPTARG";;
      n) name="$OPTARG";;
      c) control_file="$OPTARG";;
      t) test_file="$OPTARG" ;;
      u) transcripts="$OPTARG" ;;
      o) output_dir="$OPTARG" ;;
      y) odds="$OPTARG" ;;
      z) padj="$OPTARG" ;;
      m) runmode="$OPTARG";;
      ?) usage ;; # Print usage in case parameter is non-existent
   esac
done


echo "$reference_file"

## CHECK THAT MANDATORY FLAGS ARE PROVIDED
## Note this should provide compatibility with old BASH versions (i.e. < v4.2)

if [ ! -z "reference_file" ]; then
echo "-r [flag] is set"
else
echo "-r [flag] is required"
exit
fi

if [ ! -z "test_file" ]; then
echo "-t [flag] is set"
else
echo "-t [flag] is required"
exit
fi

if [ ! -z "control_file" ]; then
echo "-c [flag] is set"
else
echo "-c [flag] is required"
exit
fi

if [ ! -z "output_dir" ]; then
echo "-o [flag] is set"
else
echo "-o [flag] is required"
exit
fi


### APPLY VALUES TO OPTIONAL PARAMETERS OR REVERT TO DEFAULTS
if [ ! -z "odds" ]; then
  echo "using user-specified value for Odds Ratio Test"
else
  odds=1.5
  echo "using default value of 1.5 for Odds Ratio Test"
fi

if [ ! -z "padj" ]; then
  echo "using user-specified value for padj cutoff"
else
  padj="0.05"
  echo "using default value of 0.05 for padj cutoff"
fi



### DETERMINE RUN PATH (EXOME VS ISOFORM)
if [ "$runmode" == "exome" ]; then
  if [ ! -z "name" ]; then
     echo "-m set to exome"
  else	
     echo "-m exome requires -n [flag]"
     exit
  fi



#echo "Hallelujah"
"$DIR"/core/drummer-core-exome.sh $reference_file $name $test_file $control_file $output_dir $odds $padj
#xargs -P 8 -n 1 ./core/drummer-core-exome.sh $reference_file $name $test_file $control_file $output_dir $log2fc $odds $padj

elif [ "$runmode" == "isoform" ]; then
  if [ ! -z "transcripts" ]; then
     echo "-m set to isoform"
  else
     echo "-m isoform requires -u [flag]"
     exit
  fi




#echo "Hallelujah"
#xargs -P 8 -n 1 ./core/drummer-core-isoform.sh $reference_file $transcripts $test_file $control_file $output_dir $log2fc $odds $padj
"$DIR"/core/drummer-core-isoform.sh $reference_file $transcripts $test_file $control_file $output_dir $odds $padj
else
echo "ERROR:: User must specify runmode as -m exome|isoform"
echo ""
usage
fi










### TEST RUN
# ./drummer.sh -r /gpfs/data/depledgelab/reference-genomes/Adenovirus-Ad5.fasta -n Ad5 -c /gpfs/scratch/depled01/m6A/Ad5_topstrand.M3P1.sorted.bam -t /gpfs/scratch/depled01/m6A/Ad5_topstrand.M3KO1.sorted.bam -o /gpfs/scratch/depled01/DRUMTEST2 -m exome



