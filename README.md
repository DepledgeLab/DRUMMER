# DRUMMER
DRUMMER is designed to identify RNA modifications at nucleotide-level resolution on distinct transcript isoforms through the comparative analysis of basecall errors in Nanopore direct RNA sequencing (DRS) datasets. 

DRUMMER was designed and implemented by Jonathan S. Abebe and [Daniel P. Depledge](https://med.nyu.edu/faculty/daniel-p-depledge)

## Updates
DRUMMER v0.2 will shortly be released with the following improvements
- a single output summary file containing all candidate sites
- improved visualizations for both transcript and genome level analyses
- improved speed through parallelization
- improved tutorials for both human and viral datasets

## Table of contents
- [Introduction](#introduction)
- [Installation](#installation)
  * [Pre-requisites](#pre-requisites)
  * [Installation](#installation-1)
- [Running DRUMMER](#running-drummer)
- [Output](#output)
- [Running DRUMMER with the test datasets](#running-drummer-with-the-test-datasets)
- [Data preparation](#data-preparation)
  * [Alignment and filtering](#alignment-and-filtering)
  * [Setting up a transcript list file (isoform mode only)](#setting-up-a-transcript-list-file--isoform-mode-only-)
- [Troubleshooting](#troubleshooting)
- [Wisdom](#wisdom)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>



## Introduction
Experimental design defines the context of analysis and the modifications likely to be identified. For instance, comparison of DRS datasets from a parental cell line and METTL3-knockout cell line allows for the detection of m6A.


## Installation 

### Pre-requisites
DRUMMER requires the following packages to be installed and available in the users path: 
[SAMTools](http://www.htslib.org/) v1.3 or higher
[BEDTools](https://bedtools.readthedocs.io/en/latest/) v2.26 or higher

BASH v4.2 or higher

Python3 and modules: pandas, numpy, scipy

### Installation
git clone https://github.com/DepledgeLab/DRUMMER

## Running DRUMMER
DRUMMER requires two co-ordinate sorted and indexed BAM files as input. These should contain read alignments for the test (RNA modification absent) and control (RNA modification present) datasets (see Data Preparation section below). DRUMMER can be run in either exome or isoform mode. Exome mode (-m exome) uses DRS read alignments against the genome of a given organism to identify putatively modified bases while isoform mode (-m isoform) uses DRS read alignments against the transcriptome of a given organism to provide a high resolution mapping. While isoform mode is superior, it is also slower. 

Usage:
```
Usage: drummer.sh -r [FASTA] -u|-n [TARGETS] -c [CONTROL] -t [TEST] -o [OUTPUT] -m [RUNMODE] (OPTIONS)
```
Required flags
```
-r              fasta format reference genome (exome) or transcriptome (isoforms)

-u              list of transcripts (isoform) to be examined (single column or six-column format)
OR
-n              name of genome (exome) - must match fasta file header

-c              sorted.bam file - control (RNA modification(s) present)
-t              sorted.bam file - treatment (RNA modification(s) absent)
-o              output directory
-m              runmode (exome|isoform)

```
Optional flags
```
-y              specify odds ratio requirement (default >= 1.5)
-z              specify adjusted p_value requirement for both G-test and Odds Ratio (default<= 0.05)
-a              m6A mode (default = True), set to False to ignore m6A information
-d              reference fraction difference between unmodified and modified (default = 0.01)
-f              Candidate site visualization (default = False), set to True to visualize candidate calls for each individual transcript
```

## Output

When run to completion, DRUMMER generates a single tab-seperated text file (summary.txt) containing all predicted candidate RNA modification sites with contextual information (genome position, isoform position, sequence motif, etc). When run in m6A mode (-a TRUE), a distribution plot is also generated in an accompanying .pdf file (summary_visualization_m6A.pdf). The output directory 'complete_analysis' contains individual data files for each reference sequence provided. A second (optionall) directory 'visualization' contains individual plots of G-test scores versus position for each individual reference sequence. 

A detailed description of column headers in the summary.txt file is shown below. For the individual outputs, please see the accompanying file 'individual_output_headers.txt' for a full description of headers. 
```
transcript_id: name of transcript (isoform mode only)
chromosome: name of chromosome
reference_base: reference nucleotide at this position
pos_mod: position of nucleotide on transcript (isoform mode) or genome (exome mode)
depth_mod: read depth at this position (RNA modification present dataset)
ref_fraction_mod: fraction of reads with reference base at this position (RNA modification present dataset)
depth_unmod: read depth at this position (RNA modification absent dataset)
ref_fraction_unmod: fraction of reads with reference base at this position (RNA modification absent dataset)
frac_diff: difference between ref_fraction_unmod and ref_fraction_mod
odds_ratio: odds ratio
OR_padj: odds ratio adjusted p-value (bonferroni)
eleven_bp_motif: sequence (11-mer) centered on current position
G_test: result of 2x5 G-test
G_padj: G-test adjusted p-value (bonferroni)
candidate_site: values limited to candidate, [candidate masked], or empty based on cutoffs chosen
nearest_ac: (m6A only) distance (nt) to nearest AC motif (-ve indicates upstream, +ve indicates downstream)
nearest_ac_motif: (m6A only) sequence (5-mer) of nearest AC motif (centered on A)
genomic_position: position of nucleotide on genome (isoform mode)
```

## Running DRUMMER with the test datasets

A simple test dataset is included in the DRUMMER repository and can be used to verify DRUMMER is working correctly in your environment. Both exome-mode and isoform-mode analyses should complete in a matter of minutes.

exome mode ( runtime < 20 mins on a 'standard' desktop )
```
./drummer.sh -r TESTDATA/Adenovirus-Ad5.fasta -n Ad5 -c TESTDATA/exome.Ad5.MOD.bam -t TESTDATA/exome.Ad5.UNMOD.bam -o OUTPUTDIR_exome -m exome
```

isoform mode ( runtime < 20 mins on a 'standard' desktop )
```
./drummer.sh -r TESTDATA/Ad5_v9.1_complete.fasta -u TESTDATA/list.txt -c TESTDATA/isoform.sample.MOD.bam -t TESTDATA/isoform.sample.UNMOD.bam -o OUTPUTDIR_isoform -m isoform
```
 
## Data preparation

### Alignment and filtering

DRUMMER requires sorted.bam files containing transcriptome- or genome-level read alignments for each of the two experimental conditions being compared. Note that this is **_the_** most critical consideration when running DRUMMER (or comparable tools). Working to the adage that what you give is what you get, it is vital that **_you_** are confident that the following conditions are satisified.

- each read has a single (primary) alignment. All secondary/supplemental alignments have been filtered out
- each read is aligned against the correct transcript sequence (transcriptome analysis only)

The choice of alignment & filtering parameters for generated the input sorted.bam files will always be context dependent. For instance, our work on [mapping the locations of m6A modifications on the adenovirus Ad5 transcriptome](https://www.biorxiv.org/content/10.1101/865485v1) required careful filtering of sequence alignments due to the preence of numerous overlapping transcripts that shared the same 5' and/or 3' ends. Following a detailed [characterization of the Adenovirus AD5 transcriptome](https://www.biorxiv.org/content/10.1101/2019.12.13.876037v1), we subsequently used [minimap2](https://github.com/lh3/minimap2) to align our nanopore DRS datasets against the transcriptome as follows:

```
minimap2 -t 8 -ax map-ont -p 0.99 Ad5.transcriptome.fasta dataset1.reads.fq > dataset1.aligned.sam

samtools view -h dataset1.aligned.sam | awk ' ( $5 > 0 || $1 ~ /SQ/ ) ' | awk ' ( $2 == 0 || $1 ~ /SQ/ ) ' | awk ' ( $6 !~ /H/ || $1 ~ /SQ/ ) ' | samtools view -b - > dataset1.aligned.bam

samtools sort -o dataset1.aligned.sorted.bam dataset1.aligned.bam

samtools index dataset1.aligned.sorted.bam

```

### Setting up a transcript list file (isoform mode only)
When running DRUMMER in isoform mode, a list of transcript IDs should be supplied with the -u flag. Note these IDs must match exactly the headers present in the transcriptome database. Each ID should be supplied on seperate lines (single-column) e.g.
```
E1A-s
E1A-l
E1A-9s
```
While isoform mode identifies modified bases in a given transcript and reports the position of the modification within the transcript (e.g. position 142), additional genomic context is often useful, particularly where transcripts overlap (i.e. two transcripts may have modification at different positions within each transcript but at the same genomic location). DRUMMER can calculate the genomic position of all modifications if provided with a six-column file derived from a BED12 file. This can be easily generated by aligning the transcriptome multi-fasta file against the relevant genome fasta file, parsing to BED12 format, and then extracting relevant columns e.g.
```
### Align transcript sequences (e.g. TESTDATA/Ad5_v9.1_complete.fasta) against the genome (e.g. TESTDATA/Adenovirus-Ad5.fasta). Note that -uf flag should be replaced with -un when working with coronaviruses or other viruses with non-canonical splicing
minimap2 -ax splice -k14 -uf --secondary=no Adenovirus-Ad5.fasta Ad5_v9.1_complete.fasta > t2g.sam

### Convert to BAM (retaining only primary alignments)
samtools view -F 2304 -b -o t2g.bam t2g.sam
samtools sort -o t2g.sorted.bam t2g.bam
samtools index t2g.sorted.bam

### Convert to BED
bamToBed -bed12 -i t2g.sorted.bam > t2g.sorted.bed

### Extract relevant columns to transcripts.txt input file
cut -f1,4,6,7,10,11,12 t2g.sorted.bed > transcripts.txt
```


## Troubleshooting

COMING SOON...



## Wisdom
Detim imim finyish du wa ting, im ye fo sèmpere





