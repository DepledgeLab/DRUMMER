usage(){
echo "
DRUMMER (2020 - present)

Description:  Fast and accurate splice-aware read aligner.
Please read bbmap/docs/guides/BBMapGuide.txt for more information.

To index:     bbmap.sh ref=<reference fasta>
To map:       bbmap.sh in=<reads> out=<output sam>
To map without writing an index:
    bbmap.sh ref=<reference fasta> in=<reads> out=<output sam> nodisk

in=stdin will accept reads from standard in, and out=stdout will write to
standard out, but file extensions are still needed to specify the format of the
input and output files e.g. in=stdin.fa.gz will read gzipped fasta from
standard in; out=stdout.sam.gz will write gzipped sam.

Indexing Parameters (required when building the index):
nodisk=f                Set to true to build index in memory and write nothing
                        to disk except output.
ref=<file>              Specify the reference sequence.  Only do this ONCE,
                        when building the index (unless using 'nodisk').
build=1                 If multiple references are indexed in the same directory,
                        each needs a unique numeric ID (unless using 'nodi
Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter
any problems, or post at: http://seqanswers.com/forums/showthread.php?t=41057




Usage: [-r REFERENCE] [-c CONTROL] [-t TEST]

Required flags:
-r 		fasta format reference genome (exome) or transcriptome (isoforms)

-l 		list of transcripts to be examined (single column or XX-column format)
OR
-n		name of genome (exome) - must match fasta file header

-c 		sorted.bam file - control (RNA modification(s) present)
-t 		sorted.bam file - treatment (RNA modification(s) absent)
-o 		output directory
-m 		runmode (exome|isoform)\n

Tuning flags:
log2fc		specify log2fc required (default >= 0.5)
odds            specify odds ratio requirement (default >= 1)
padj            specify adjusted p_value (G-test) requirement (default<= 0.05)

DRUMMER was written by Jonathan S. Abebe & Daniel P. Depledge.
If you encounter any problems with DRUMMER, please log them as an issue at the GitHub page https://github.com/DepledgeLab/DRUMMER
"
}


### Function for showing usage
if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
        usage
        exit
fi


### GETOPTS
while getopts ":rctomnlxyz:" opt; do
   case "$opt" in
      r ) reference_file="$OPTARG" ;;
      n ) name="$OPTARG" ;;
      c ) control_file="$OPTARG" ;;
      t ) test_file="$OPTARG" ;;
      l ) transcripts="$OPTARG" ;;
      o ) output_dir="$OPTARG" ;;
      x ) log2fc="$OPTARG" ;;
      y ) odds="$OPTARG" ;;
      z ) padj="$OPTARG" ;;
      m ) runmode="$OPTARG" ;;
      ? ) usage ;; # Print usage in case parameter is non-existent
   esac
done

### CHECK THAT MANDATORY FLAGS ARE PROVIDED
if [ "x" != "x$reference_file" ]; then
  echo "-r [flag] is required"
  exit
fi

if [ "x" != "x$test_file" ]; then
  echo "-t [flag] is required"
  exit
fi

if [ "x" != "x$control_file" ]; then
  echo "-l [flag] is required"
  exit
fi

if [ "x" != "x$output_dir" ]; then
  echo "-o [flag] is required"
  exit
fi


### APPLY VALUES TO OPTIONAL PARAMETERS OR REVERT TO DEFAULTS

if [ "x" != "x$odds" ]; then
  odds=1
  exit
fi

if [ "x" != "x$padj" ]; then
  padj="0.05"
  exit
fi

if [ "x" != "x$log2fc" ]; then
  log2fc="0.5"
  exit
fi



### DETERMINE RUN PATH (EXOME VS ISOFORM)
if [[ $runmode == "exome" ]]; then
  if [ "x" != "x$name" ]; then
    echo "-m exome requires -n [flag]"
    exit
  else
    echo "Hallelujah"
    #xargs -P 8 -n 1 ./core/drummer-core-exome.sh $reference_file $name $test_file $control_file $output_dir $log2fc $odds $padj
  fi
elif [[ $runmode == "isoform" ]]; then
  if [ "x" != "x$transcripts" ]; then
    echo "-m exome requires -x [flag]"
    exit
  else
    echo "Hallelujah"
    #xargs -P 8 -n 1 ./core/drummer-core-isoform.sh $reference_file $transcripts $test_file $control_file $output_dir $log2fc $odds $padj
  fi
else
echo "ERROR:: User must specify either -runmode exome or runmode isoform"
echo ""
usage
fi

### EXOME TEST 
#./wrapper.sh -r /gpfs/data/depledgelab/reference-genomes/Adenovirus-Ad5.fasta -n Ad5 -c /gpfs/scratch/depled01/m6A/Ad5_topstrand.M3P1.sorted.bam -t /gpfs/scratch/depled01/m6A/Ad5_topstrand.M3KO1.sorted.bam -o /gpfs/scratch/depled01/DRUMTEST2 -m exome

### ISOFORM TEST







