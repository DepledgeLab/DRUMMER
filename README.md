# DRUMMER
DRUMMER is a software packaged designed to identify RNA modifications at both transcript and nucleotide-level resolution through the comparative analysis of basecall errors in Nanopore direct RNA sequencing (DRS) datasets. 

DRUMMER was designed and implemented by Jonathan Abebe and [Daniel Depledge](https://med.nyu.edu/faculty/daniel-p-depledge)


## Introduction
Experimental design defines the context of analysis and the modifications likely to be identified. For instance, comparison of DRS datasets from a parental cell line and METTL3-knockout cell line allows for the detection of m6A.


## Installation 

### Pre-requisites
DRUMMER requires the following packages to be installed and available in the users path: 
[SAMTools](http://www.htslib.org/) v1.3 or higher
[BEDTools](https://bedtools.readthedocs.io/en/latest/) v2.26 or higher

BASH v4.2 or higher

Python3 and modules: pandas, numpy, scipy

### Setting up DRUMMER
git clone https://github.com/DepledgeLab/DRUMMER

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
cut -f4,6,7,10,11,12 t2g.sorted.bed > transcripts.txt

```

## Running DRUMMER
DRUMMER can be run in either exome or isoform mode. Exome mode (-m exome) uses DRS read alignments against the genome of a given organism to identify putatively modified bases while isoform mode (-m isoform) uses DRS read alignments against the transcriptome of a given organism to provide a high resolution mapping. While isoform mode is superior, it is also slower. 

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
-z              specify adjusted p_value (G-test) requirement (default<= 0.05)
-a              m6A mode (default = True), set to False to ignore m6A information
-d              reference fraction difference between unmodified and modified (default = 0.01)
-f              Candidate site visualization (default = False), set to True to visualize candidate calls
```

## Output

When run to completion, DRUMMER outputs a single text file containing the results. A detailed description of the column headers follow.
```
chr_unmod: name of chromosome or transcript
pos_unmod: position of nucleotide on chromosome/transcript
ref_unmod: nucleotide at this position
depth_unmod: read depth at this position
A_unmod: number of reads supporting an Adenine at this position
C_unmod: number of reads supporting an Cytosine at this position
G_unmod: number of reads supporting an Guanine at this position
T_unmod: number of reads supporting an Thymine at this position
N_unmod: number of reads supporting an indel (N) at this position
ref_fraction_unmod: fraction of reads matching the reference nucleotide
chr_mod: name of chromosome or transcript
pos_mod: position of nucleotide on chromosome/transcript
ref_mod: nucleotide at this position
depth_mod: read depth at this position
A_mod: number of reads supporting an Adenine at this position
C_mod: number of reads supporting an Cytosine at this position
G_mod: number of reads supporting an Guanine at this position
T_mod: number of reads supporting an Thymine at this position
N_mod: number of reads supporting an indel (N) at this position
ref_fraction_mod: fraction of reads matching the reference nucleotide
ratio_unmod: 
ratio_mod: 
fold_change: fold change difference between ratio_unmod and ratio_mod
log2_fc: log2 fold change
odds_ratio: odds ratio
p_values_OR: pvalue calculated by odds_ratio
nearest_ac: (m6A only) distance (nt) to nearest AC motif (-ve indicates upstream, +ve indicates downstream)
nearest_ac_motif: (m6A only) sequence (5-mer) of nearest AC motif (centered on A)
five_bp_motif: sequence (5-mer) centered on current position
eleven_bp_motif: sequence (11-mer) centered on current position
G_test: Result of 2x5 G-test
p_val: p-value of G-test
padj: bonferroni-corrected p-value 
p_values_OR_adj: bonferroni-corrected p-value (odds ratio)
candidate_site: Nucleotide predicted to be modified based on supplied cutoffs for padj, odds_ratio, and log2_fc
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
 


## Troubleshooting

COMING SOON...



### Wisdom
Detim imim finyish du wa ting, im ye fo sèmpere







