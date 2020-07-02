#!/bin/bash 


#### EXOME LEVEL ANALYSIS

genome_file=$1 ### Input genome FASTA file
name=$2 ### working ID - must match FASTA header
test_file=$3 ### INPUT CTRL (RNA MOD PRESENT) SORTED BAM FILE
control_file=$4  ### INPUT CTRL (RNA MOD ABSENT) SORTED BAM FILE
output_dir=$5
log2fc=$6
odds=$7
padj=$8

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

if [ -d "$output_dir"/gTest ]; then
    echo "Directory /path/to/dir exists."
else
    mkdir "$output_dir"/gTest
fi

if [ -d "$output_dir"/odds_ratio ]; then
    echo "Directory /path/to/dir exists."
else
    mkdir "$output_dir"/odds_ratio
fi

if [ -d "$output_dir"/motif_information ]; then
    echo "Directory /path/to/dir exists."
else
    mkdir "$output_dir"/motif_information
fi

if [ -d "$output_dir"/merged ]; then
    echo "Directory /path/to/dir exists."
else
    mkdir "$output_dir"/merged
fi

### ONE LINER TO DETERMINE SEQUENCE LENGTH
length=$(awk '/^>/ {if (seqlen){print seqlen}; seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' "$genome_file")

#### UNIQUE MAP FILTER METHOD NOTE: REQUIRES -alt.sorted.bam
#samtools view -b "$test_file" "$name":1-"$length" -o "$output_dir"/map/"$name".UNMOD.bam
#samtools sort -o "$output_dir"/map/"$name".UNMOD.sorted.bam "$output_dir"/map/"$name".UNMOD.bam
#samtools index "$output_dir"/map/"$name".UNMOD.sorted.bam

#samtools view -b "$control_file" "$name":1-"$length" -o "$output_dir"/map/"$name".MOD.bam
#samtools sort -o "$output_dir"/map/"$name".MOD.sorted.bam "$output_dir"/map/"$name".MOD.bam
#samtools index "$output_dir"/map/"$name".MOD.sorted.bam

#samtools index "$output_dir"/map/"$name".CTRL.sorted.filtered.bam
#samtools index "$output_dir"/map/"$name".TEST.sorted.filtered.bam

#bam-readcount -f "$genome_file" "$output_dir"/map/"$name".TEST.sorted.filtered.bam "$name" > "$output_dir"/bam_readcount/"$name".TEST.bamreadcount.txt
#bam-readcount -f "$genome_file" "$output_dir"/map/"$name".CTRL.sorted.filtered.bam "$name" > "$output_dir"/bam_readcount/"$name".CTRL.bamreadcount.txt


#### JONATHAN TO CHECK BELOW

input_bamreadcounts=$output_dir/bam_readcount/$name
###echo $transcript_name
###echo $control_file
###echo $test_file

#Creates a directory called filter and subdir of transcript name
python3 ../modules/readcount_filter.py -i $input_bamreadcounts

#paste -d "\t" "$input_bamreadcounts".TEST.filtered.txt "$input_bamreadcounts".TEST.filtered.txt > "$input_bamreadcounts".merged.filtered.txt

#odds_ratio / make new directory with odds ratio added to individual reads
#merged_transcripts=merged/$transcript_name.*
merged_transcripts=$output_dir/merged/$name.*
#echo $transcript_name
#echo $merged_transcripts
python3 ../modules/odds_ratio.py -i $merged_transcripts

#motif_information / make new directory with motif information added to individual reads
#python3 motif_information.py -i odds_ratio/AdPol.merged.filtered.odds_ratio.txt

motif_transcripts=$output_dir/odds_ratio/$name.*
#motif_transcripts=odds_ratio/$transcript_name.*
#echo $motif_transcripts
python3 ../modules/motif_information.py -i $motif_transcripts

#Gtest/ make new directory of gtest values added to individual reads
gtest_transcripts=$output_dir/motif_information/$transcript_name.*
#gtest_transcripts=motif_information/$transcript_name.*
#python3 create_output_file.py -i _ -o gTest
#echo $gtest_transcripts
# Rscript ../modules/Gtest.R $gtest_transcripts "$output_dir"/gTest/$transcript_name.merged.odds_ratio.motif_information.gtest.csv
python3 ../modules/Gtest.py -i $gtest_transcripts

candidate_transcripts=$output_dir/gTest/$transcript_name.*
python3 ../modules/find_candidates.py -i $candidate_transcripts -r $odds -l $log2fc -p $padj

