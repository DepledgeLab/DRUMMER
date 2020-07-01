# DRUMMER
Detection of RNA modifications mediated by error rates (DRUMMER) is a software packaged designed to identify RNA modifications at both transcript and nucleotide-level resolution through the comparative analysis of Nanopore direct RNA sequencing (DRS) datasets. 

DRUMMER was designed and implemented by Jonathan Abebe and [Daniel Depledge](https://med.nyu.edu/faculty/daniel-p-depledge)


## Introduction
Experimental design defines the context of analysis and the modifications likely to be identified. For instance, comparison of DRS datasets from a parental cell line and METTL3-knockout cell line allows for the detection of m6A.


## Installation 

### Pre-requisites
DRUMMER requires the following packages to be installed and available in the users path: 
[SAMTools](http://www.htslib.org/) v1.3 or higher
[bam-readcount](https://github.com/genome/bam-readcount) v0.8 or higher

Python3 and modules: pandas, numpy

### Setting up DRUMMER
git clone https://github.com/DepledgeLab/DRUMMER

COMING SOON...

## Data preparation

### Alignment and filtering

DRUMMER requires sorted.bam files containing transcriptome- or genome-level read alignments for each of the two experimental conditions being compared. Note that this is **_the_** most critical consideration when running DRUMMER (or comparable tools). Working to the adage that what you give is what you get, it is vital that **_you_** are confident that the following conditions are satisified.

- each read has a single (primary) alignment. All secondary/supplemental alignments have been filtered out
- each read is aligned against the correct transcript sequence (transcriptome analysis only)

The choice of alignment & filtering parameters for generated the input sorted.bam files will always be context dependent. For instance, our work on [mapping the locations of m6A modifications on the adenovirus Ad5 transcriptome](https://www.biorxiv.org/content/10.1101/865485v1) required careful filtering of sequence alignments due to the preence of numerous overlapping transcripts that shared the same 5' and/or 3' ends. Following a detailed [characterization of the Adenovirus AD5 transcriptome](https://www.biorxiv.org/content/10.1101/2019.12.13.876037v1), we subequently aligned our nanopore DRS datasets against the transcriptome as follows

```
minimap2 -t 8 -ax map-ont -p 0.99 Ad5.transcriptome.fasta dataset1.reads.fq > dataset1.aligned.sam

samtools view -h dataset1.aligned.sam | awk ' ( $5 > 0 || $1 ~ /SQ/ ) ' | awk ' ( $2 == 0 || $1 ~ /SQ/ ) ' | awk ' ( $6 !~ /H/ || $1 ~ /SQ/ ) ' | samtools view -b - > dataset1.aligned.bam

samtools sort -o dataset1.aligned.sorted.bam dataset1.aligned.bam

samtools index dataset1.aligned.sorted.bam
```

### Setting up a transcript list file

COMING SOON...



## Running DRUMMER
DRUMMER can be run in either exome or isoform mode. Exome mode (-m exome) uses DRS read alignments against the genome of a given organism to identify putatively modified bases while isoform mode (-m isoform) uses DRS read alignments against the transcriptome of a given organism to provide a high resolution mapping. While isoform mode is superior, it is also slower. 

Usage:
```
Usage: drummer.sh -r [FASTA] -l|-n [TARGETS] -c [CONTROL] -t [TEST] -o [OUTPUT] -m [RUNMODE] (OPTIONS)
```
Required flags
```
-r              fasta format reference genome (exome) or transcriptome (isoforms)

-u              list of transcripts (isoform) to be examined (single column or XX-column format)
OR
-n              name of genome (exome) - must match fasta file header

-c              sorted.bam file - control (RNA modification(s) present)
-t              sorted.bam file - treatment (RNA modification(s) absent)
-o              output directory
-m              runmode (exome|isoform)

```
Optional flags
```
-x              specify log2fc required (default >= 0.5)
-y              specify odds ratio requirement (default >= 1)
-z              specify adjusted p_value (G-test) requirement (default<= 0.05)
```






## Interpreting output

COMING SOON...


## Troubleshooting

COMING SOON...




Detim imim finyish du wa ting, im ye fo sèmpere







