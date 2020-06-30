#!/bin/bash 


#### TRANSCRIPT SPECIFIC ANALYSIS

test_file=$1
control_file=$2
transcript_name=$3


### NEED FUNCTION TO SPLITS PROVIDED FASTA FILE 



### ONE LINER TO DETERMINE SEQUENCE LENGTH
length=$(awk '/^>/ {if (seqlen){print seqlen}; seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' transcripts/"$transcript_name".fa)

mkdir -p map
#### UNIQUE MAP FILTER METHOD NOTE: REQUIRES -alt.sorted.bam
samtools view -b "$test_file" "$transcript_name":1-"$length" -o map/"$transcript_name".M3P1.bam
samtools sort -o map/"$transcript_name".M3P1.sorted.bam map/"$transcript_name".M3P1.bam
samtools index map/"$transcript_name".M3P1.sorted.bam

samtools view -b "$control_file" "$transcript_name":1-"$length" -o map/"$transcript_name".M3KO1.bam
samtools sort -o map/"$transcript_name".M3KO1.sorted.bam map/"$transcript_name".M3KO1.bam
samtools index map/"$transcript_name".M3KO1.sorted.bam

mkdir -p bam_readcount
$bam_readcount -f transcripts/"$transcript_name".fa map/"$transcript_name".M3P1.sorted.bam "$transcript_name" > bam_readcount/"$transcript_name".M3P1.bamreadcount.txt
$bam_readcount -f transcripts/"$transcript_name".fa map/"$transcript_name".M3KO1.sorted.bam "$transcript_name" > bam_readcount/"$transcript_name".M3KO1.bamreadcount.txt





input_bamreadcounts=bam_readcount/$transcript_name
#echo $transcript_name
#echo $control_file
#echo $test_file

#Creates a directory called filter and subdir of transcript name
python3 readcount_filter.py -i $input_bamreadcounts

#paste -d "\t" "$input_bamreadcounts".M3P1.filtered.txt "$input_bamreadcounts".M3KO1.filtered.txt > "$input_bamreadcounts".merged.filtered.txt

#odds_ratio / make new directory with odds ratio added to individual reads
#merged_transcripts=merged/$transcript_name.*
merged_transcripts=merged/$transcript_name.*
#echo $transcript_name
#echo $merged_transcripts
python3 odds_ratio.py -i $merged_transcripts

#motif_information / make new directory with motif information added to individual reads
#python3 motif_information.py -i odds_ratio/AdPol.merged.filtered.odds_ratio.txt

motif_transcripts=odds_ratio/$transcript_name.*
#motif_transcripts=odds_ratio/$transcript_name.*
#echo $motif_transcripts
python3 motif_information.py -i $motif_transcripts



#Gtest/ make new directory of gtest values added to individual reads
gtest_transcripts=motif_information/$transcript_name.*
#gtest_transcripts=motif_information/$transcript_name.*
#python3 create_output_file.py -i _ -o gTest
#echo $gtest_transcripts
mkdir -p gTest
Rscript Gtest.R $gtest_transcripts gTest/$transcript_name.merged.odds_ratio.motif_information.gtest.csv
