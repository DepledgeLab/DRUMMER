#!/bin/bash 

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"

echo "$DIR"


#### EXOME LEVEL ANALYSIS

genome_file=$1 ### Input genome FASTA file
name=$2 ### working ID - must match FASTA header
test_file=$3 ### INPUT CTRL (RNA MOD PRESENT) SORTED BAM FILE
control_file=$4  ### INPUT CTRL (RNA MOD ABSENT) SORTED BAM FILE
output_dir=$5
log2fc=$6
odds=$7
padj=$8
m6A_status=$9 
fraction_diff=${10} 
visualization=${11}
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

if [ -d "$output_dir"/complete_analysis ]; then
    echo "Directory /path/to/dir exists."
else
    mkdir "$output_dir"/complete_analysis
fi


## ONE LINER TO DETERMINE SEQUENCE LENGTH
length=$(awk '/^>/ {if (seqlen){print seqlen}; seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' "$genome_file")

#### UNIQUE MAP FILTER METHOD NOTE: REQUIRES -alt.sorted.bam
samtools view -b "$test_file" "$name":1-"$length" -o "$output_dir"/map/"$name".TEST.bam
samtools sort -o "$output_dir"/map/"$name".TEST.sorted.bam "$output_dir"/map/"$name".TEST.bam
samtools index "$output_dir"/map/"$name".TEST.sorted.bam

samtools view -b "$control_file" "$name":1-"$length" -o "$output_dir"/map/"$name".CTRL.bam
samtools sort -o "$output_dir"/map/"$name".CTRL.sorted.bam "$output_dir"/map/"$name".CTRL.bam
samtools index "$output_dir"/map/"$name".CTRL.sorted.bam

"$DIR"/../modules/bam-readcount -f "$genome_file" "$output_dir"/map/"$name".TEST.sorted.bam "$name" > "$output_dir"/bam_readcount/"$name".TEST.bamreadcount.txt
"$DIR"/../modules/bam-readcount -f "$genome_file" "$output_dir"/map/"$name".CTRL.sorted.bam "$name" > "$output_dir"/bam_readcount/"$name".CTRL.bamreadcount.txt


#### JONATHAN TO CHECK BELOW

input_bamreadcounts="$output_dir"/bam_readcount/$name
#echo $transcript_name
#echo $control_file
#echo $test_file

#Creates a directory called filter and subdir of transcript name

if [ -f "$input_bamreadcounts.TEST.bamreadcount.txt" -a -f "$input_bamreadcounts.CTRL.bamreadcount.txt" ]
then
python3 "$DIR"/../modules/readcount_filter.py -i $input_bamreadcounts -o $output_dir
else 
    echo "Failed at bam_readcount"
fi

#merged_transcripts=merged/$transcript_name.*
merged_transcripts="$output_dir"/merged/$name.*
#echo $transcript_name
#echo $merged_transcripts
if [ -f $merged_transcripts ] 
then
python3 "$DIR"/../modules/odds_ratio.py -i $merged_transcripts -o $output_dir
else 
    echo "Failed at odds ratio"
fi

#motif_information / make new directory with motif information added to individual reads
#python3 motif_information.py -i odds_ratio/AdPol.merged.filtered.odds_ratio.txt

motif_transcripts="$output_dir"/odds_ratio/$name.*
#motif_transcripts=odds_ratio/$transcript_name.*
#echo $motif_transcripts
if [ -f $motif_transcripts ]; then
python3 "$DIR"/../modules/motif_information.py -i $motif_transcripts -o $output_dir -m $m6A_status
else 
    echo "Failed at motif information"
fi

#Gtest/ make new directory of gtest values added to individual reads
gtest_transcripts="$output_dir"/motif_information/$name.*
#gtest_transcripts=motif_information/$transcript_name.*
#python3 create_output_file.py -i _ -o gTest
#echo $gtest_transcripts
# Rscript ../modules/Gtest.R $gtest_transcripts "$output_dir"/gTest/$transcript_name.merged.odds_ratio.motif_information.gtest.csv

if [ -f $gtest_transcripts ]; then
python3 "$DIR"/../modules/Gtest.py -i $gtest_transcripts -o $output_dir
else 
    echo "Failed at gtest"
fi

candidate_transcripts="$output_dir"/gTest/$name.*
# python3 "$DIR"/../modules/find_candidates.py -i $candidate_transcripts -r $odds -p $padj -d $fraction_diff -o $output_dir/$name.complete.txt
if [ -f $candidate_transcripts ]; then
python3 "$DIR"/../modules/find_candidates.py -i $candidate_transcripts -r $odds -p $padj -d $fraction_diff -o $output_dir/complete_analysis/$name.complete.txt
else 
    echo "Failed at find_candidates"
fi
#python3 "$DIR"/../modules/genomic_locations.py -i $output_dir/$id.complete.txt -t $list -o $output_dir/$id.complete.txt

if [ "$(ls -A $output_dir/complete_analysis/)" ] 
then
python3 "$DIR"/../extras/summary.py -i $output_dir/complete_analysis/ -o $output_dir/summary.txt -m $m6A_status
fi


if [ -f $output_dir/complete_analysis/$name.complete.txt ]; then
if [ $visualization = "True" ]
then
python3 "$DIR"/../modules/candidate_visualization.py -i $output_dir/complete_analysis/$name.complete.txt -o $output_dir -m $m6A_status -t exome
fi 
fi


if [ -f "$output_dir/summary.txt" ]; then
if [ $m6A_status = "True" ]
then
python3 "$DIR"/../extras/m6a_summary_plot.py -d $output_dir/complete_analysis/ -i $output_dir/summary.txt -o $output_dir/summary_visualization_m6A.pdf
fi 
fi

rm -r "$output_dir"/bam_readcount "$output_dir"/filtered "$output_dir"/gTest "$output_dir"/map "$output_dir"/merged "$output_dir"/motif_information "$output_dir"/odds_ratio 





