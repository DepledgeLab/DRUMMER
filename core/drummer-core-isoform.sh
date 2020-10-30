#!/bin/bash 

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"

echo "$DIR"

### ISOFORM LEVEL ANALYSIS

transcriptome_file=$1 ### Input transcriptome multi-FASTA file
list=$2 ### list of transcript IDs - must match FASTA header
test_file=$3 ### INPUT CTRL (RNA MOD PRESENT) SORTED BAM FILE
control_file=$4  ### INPUT CTRL (RNA MOD ABSENT) SORTED BAM FILE
output_dir=$5
log2fc=$6
odds=$7
padj=$8
m6A_status=$9 
fraction_diff=${10} 
visualization=${11}

echo $transcriptome_file $list $test_file $control_file $output_dir 

### CHECK IF OUTPUT DIRECTORIES EXIST, OTHERWISE MAKE IT
if [ -d "$output_dir" ]; then
    echo "Directory /path/to/dir exists."
else
    mkdir "$output_dir"
fi

if [ -d "$output_dir"/map ]; then
    echo "Directory /path/to/dir exists."
else
    mkdir "$output_dir"/map
fi

if [ -d "$output_dir"/bam_readcount ]; then
    echo "Directory /path/to/dir exists."
else
    mkdir "$output_dir"/bam_readcount
fi

#if [ -d "$output_dir"/gTest ]; then
#    echo "Directory /path/to/dir exists."
#else
#    mkdir "$output_dir"/gTest
#fi

if [ -d "$output_dir"/transcripts ]; then
    echo "Directory /path/to/dir exists."
else
    mkdir "$output_dir"/transcripts
fi

if [ -d "$output_dir"/complete_analysis ]; then
    echo "Directory /path/to/dir exists."
else
    mkdir "$output_dir"/complete_analysis
fi


### ADD FUNCTION TO SPLIT MULTI-FASTA INTO INDIVIDUAL TRANSCRIPTS

cat $transcriptome_file | awk '{ if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")} print $0 > "'$output_dir'/transcripts/"filename }'



### LOOP FUNCTION REQUIRED BELOW????
#FOR testing yes. FOR final version I do not think so, should handled by drummer.sh when xargs is called, this should hand parallelization as well.
#(i.e.) cat $reference_file | xargs -P 8 -n 1 ./all_test.sh $test_file $control_file

#for id in $(echo $list | cut -f1 -d$'\t'); do


while IFS=$'\t' read id remainder; 
do

echo $id

### ONE LINER TO DETERMINE SEQUENCE LENGTH
#bam_readcount=/gpfs/data/tsirigoslab/home/ja3539/Nanopore/OPEN/bam-readcount/bam-readcount/build/bin/bam-readcount

length=$(awk '/^>/ {if (seqlen){print seqlen}; seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' "$output_dir"/transcripts/"$id".fa)

#### UNIQUE MAP FILTER METHOD NOTE: REQUIRES -alt.sorted.bam
samtools view -b "$test_file" "$id":1-"$length" -o "$output_dir"/map/"$id".UNMOD.bam
samtools sort -o "$output_dir"/map/"$id".UNMOD.sorted.bam "$output_dir"/map/"$id".UNMOD.bam
samtools index "$output_dir"/map/"$id".UNMOD.sorted.bam

samtools view -b "$control_file" "$id":1-"$length" -o "$output_dir"/map/"$id".MOD.bam
samtools sort -o "$output_dir"/map/"$id".MOD.sorted.bam "$output_dir"/map/"$id".MOD.bam
samtools index "$output_dir"/map/"$id".MOD.sorted.bam

"$DIR"/../modules/bam-readcount -f "$output_dir"/transcripts/"$id".fa "$output_dir"/map/"$id".UNMOD.sorted.bam "$id" > "$output_dir"/bam_readcount/"$id".UNMOD.bamreadcount.txt
"$DIR"/../modules/bam-readcount -f "$output_dir"/transcripts/"$id".fa "$output_dir"/map/"$id".MOD.sorted.bam "$id" > "$output_dir"/bam_readcount/"$id".MOD.bamreadcount.txt

#done < $list

#### JONATHAN TO CHECK BELOW

input_bamreadcounts="$output_dir"/bam_readcount/$id
python3 "$DIR"/../modules/readcount_filter.py -i $input_bamreadcounts -o $output_dir

merged_transcripts="$output_dir"/merged/$id.*

python3 "$DIR"/../modules/odds_ratio.py -i $merged_transcripts -o $output_dir

motif_transcripts="$output_dir"/odds_ratio/$id.*

python3 "$DIR"/../modules/motif_information.py -i $motif_transcripts -o $output_dir -m $m6A_status

gtest_transcripts="$output_dir"/motif_information/$id.*

python3 "$DIR"/../modules/Gtest.py -i $gtest_transcripts -o $output_dir

candidate_transcripts="$output_dir"/gTest/$id.*

python3 "$DIR"/../modules/find_candidates.py -i $candidate_transcripts -r $odds -p $padj -d $fraction_diff -o $output_dir/complete_analysis/$id.complete.txt

python3 "$DIR"/../modules/genomic_locations.py -i $output_dir/complete_analysis/$id.complete.txt -t $list -o $output_dir/complete_analysis/$id.complete.txt

if [ $visualization = "True" ]
then
if [ $m6A_status == "True" ]
then
python3 "$DIR"/../modules/candidate_visualization.py -i $output_dir/complete_analysis/$id.complete.txt -o $output_dir
fi
fi 

##Location of transcript file?
#python3 "$DIR"/../modules/genomic_locations.py -i $output_dir/$id.complete.txt -t $list -o $output_dir/$id.complete.txt
done < $list

python3 "$DIR"/../extras/summary.py -i $output_dir/complete_analysis/ -o $output_dir/summary.txt -m m6A_status
#input_bamreadcounts=bam_readcount/$name
#echo $transcript_name
#echo $control_file
#echo $test_file

#Creates a directory called filter and subdir of transcript name
#python3 ../modules/readcount_filter.py -i $input_bamreadcounts

#paste -d "\t" "$input_bamreadcounts".TEST.filtered.txt "$input_bamreadcounts".TEST.filtered.txt > "$input_bamreadcounts".merged.filtered.txt

#odds_ratio / make new directory with odds ratio added to individual reads
#merged_transcripts=merged/$transcript_name.*
#merged_transcripts=merged/$name.*
#echo $transcript_name
#echo $merged_transcripts
#python3 ../modules/odds_ratio.py -i $merged_transcripts

#motif_information / make new directory with motif information added to individual reads
#python3 motif_information.py -i odds_ratio/AdPol.merged.filtered.odds_ratio.txt

#motif_transcripts=odds_ratio/$name.*
#motif_transcripts=odds_ratio/$transcript_name.*
#echo $motif_transcripts
#python3 ../modules/motif_information.py -i $motif_transcripts

#Gtest/ make new directory of gtest values added to individual reads
#gtest_transcripts=motif_information/$transcript_name.*
#gtest_transcripts=motif_information/$transcript_name.*
#python3 create_output_file.py -i _ -o gTest
#echo $gtest_transcripts
# Rscript ../modules/Gtest.R $gtest_transcripts "$output_dir"/gTest/$transcript_name.merged.odds_ratio.motif_information.gtest.csv
#python3 ../modules/Gtest.py -i $gtest_transcripts

#candidate_transcripts=gTest/$transcript_name.*
#python3 ../modules/find_candidates.py -i $candidate_transcripts -r $odds -l $log2fc -p $padj -o $output_dir/$transcript_name.complete.txt


