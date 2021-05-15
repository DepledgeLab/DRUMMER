import argparse
import modules.test_isoform as test_isoform
import os
import modules.support 
import modules.test_exome as test_exome


ap = argparse.ArgumentParser(description = 'Takes in the output from the pipeline and determines candidate sites \
using log2fc, odds_ratio and padj')

requiredGrp = ap.add_argument_group('required arguments')

requiredGrp = ap.add_argument_group('required arguments')
requiredGrp.add_argument("-r",'--reference_fasta', required=True, help="genome or transcriptome fasta file")
requiredGrp.add_argument("-l",'--list', required=False, help="list of chromosome/transcript IDs for analysis [ must match IDs of supplied fasta file ]")
requiredGrp.add_argument("-t",'--treatment_bam', required=True, help="sorted.bam file(s) for treatment [ e.g. knockdown, knockout, overexpression ]",action="store",nargs = "*")
requiredGrp.add_argument("-c",'--control_bam', required=True, help="sorted.bam file(s) for control",action="store",nargs = "*")
requiredGrp.add_argument("-o",'--output_dir', required=True, help="output directory (will be made / overwritten)")
# requiredGrp.add_argument("-x",'--log2fc', required=False, help="Log2FC value (not used)",default=1,type = float)
requiredGrp.add_argument("-z",'--odds', required=False, help="odds ratio cutoff (default = +/- 1.5)",default=1.5,type = float)
requiredGrp.add_argument("-p",'--padj', required=False, help="padj cutoff (default = < 0.05 for both OR and G-test)",default=0.05,type = float)
requiredGrp.add_argument("-m",'--m6A', required=False, help="run in m6A mode (default = false)")
requiredGrp.add_argument("-f",'--fraction_diff', required=False, help="Reference fraction difference (default = 0.01)",default=.01,type = float)
requiredGrp.add_argument("-v",'--visualization', required=False, help="produce visualizations for individual transcripts (default = false)")
# requiredGrp.add_argument("-i",'--additional_information', required=False, help="Additional information to display for summary file",action="store",nargs = "*")
requiredGrp.add_argument("-a",'--analysis_mode', required=True, help="isoform or exome")
requiredGrp.add_argument("-n",'--name', required=False, help="Chromosome name, exome mode")
requiredGrp.add_argument("-i",'--indel_filter', required=False, help="Filter for indels or retain")
args = vars(ap.parse_args())



modules.support.check_samtools()
#transcriptome_file,test_file,control_file,log2fc,odds,padj,m6A_status,fraction_diff,output_dir,additional_columns,mode

if args['analysis_mode'].lower() == 'exome':
    print('IN EXOME')
    transcriptome_file = args['reference_fasta']
    treatment_bam = args['treatment_bam']
    control_bam = args['control_bam']
#     log2fc = args['log2fc']
    odds = args['odds']
    padj = args['padj']
    m6A = modules.support.handle_booleans(args['m6A'],False)
    deletion = modules.support.handle_booleans(args['indel_filter'],False)
    fraction_diff = args['fraction_diff']
    output_dir = args['output_dir']
#     additional_columns = args['additional_information']
    name = args['name']
    if name is None:
        raise modules.support.MyException('Need -n argument when running in exome mode')
    modules.support.print_logo('exome\n')
    test_exome.main(transcriptome_file,treatment_bam,name,control_bam,odds,padj,m6A,fraction_diff,output_dir,'exome',deletion)
elif args['analysis_mode'].lower() == 'isoform':
    print('isoform')
    transcriptome_file = args['reference_fasta']
    treatment_bam = args['treatment_bam']
    path_transcripts = args['list']
    control_bam = args['control_bam']
    deletion = modules.support.handle_booleans(args['indel_filter'],False)
#     log2fc = args['log2fc']
    odds = args['odds']
    padj = args['padj']
    m6A = modules.support.handle_booleans(args['m6A'],False)
    fraction_diff = args['fraction_diff']
    visualization_input = modules.support.handle_booleans(args['visualization'],False)
    output_dir = args['output_dir']
#     additional_columns = args['additional_information']
    if path_transcripts is None:
        raise modules.support.MyException('Need -u argument when running in isoform mode')
    modules.support.print_logo('isoform\n')
    test_isoform.main(transcriptome_file,treatment_bam,path_transcripts,control_bam,odds,padj,m6A,fraction_diff,visualization_input,output_dir,'isoform',deletion)
    
    
    
    
    
    
    
    
    
    
    
    