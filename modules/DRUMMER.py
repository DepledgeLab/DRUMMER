from __future__ import division
import pandas as pd
import numpy as np
import argparse
import math
import re
import argparse
import os
import pandas as pd
import subprocess
import pandas as pd
import math
import numpy as np
import scipy.stats as stats
import concurrent.futures
import timeit
import re
import math
import genomic
import gtest
import merge
import motif
import candidates
import odds
import sys
import tqdm
#from colorama import Fore
from termcolor import colored
#python3 DRUMMER.py -u /gpfs/home/ja3539/depledge/Jonathan/DRUMMER/TESTDATA/Ad5.sample.transcripts.txt -r /gpfs/home/ja3539/depledge/Jonathan/DRUMMER/TESTDATA/Ad5_v9.1_complete.fasta -o test -c /gpfs/home/ja3539/depledge/Jonathan/DRUMMER/TESTDATA/isoform.Ad5.MOD.bam -t /gpfs/home/ja3539/depledge/Jonathan/DRUMMER/TESTDATA/isoform.Ad5.UNMOD.bam
# import modules.readcount_filter
# import modules.create_output
ap = argparse.ArgumentParser(description = 'Takes in the output from the pipeline and determines candidate sites \
using log2fc, odds_ratio and padj')

requiredGrp = ap.add_argument_group('required arguments')
requiredGrp.add_argument("-r",'--transcriptome_file', required=True, help="input file location")
requiredGrp.add_argument("-u",'--list', required=True, help="guide rna sequence")
requiredGrp.add_argument("-t",'--test_file', required=True, help="reference file",action="store",nargs = "*")
requiredGrp.add_argument("-c",'--control_file', required=True, help="output file location",action="store",nargs = "*")
requiredGrp.add_argument("-o",'--output_dir', required=True, help="output file location")
requiredGrp.add_argument("-x",'--log2fc', required=False, help="reference file")
requiredGrp.add_argument("-y",'--odds', required=False, help="output file location")
requiredGrp.add_argument("-z",'--padj', required=False, help="guide rna sequence")
requiredGrp.add_argument("-a",'--m6A_status', required=False, help="reference file")
requiredGrp.add_argument("-d",'--fraction_diff', required=False, help="output file location")
requiredGrp.add_argument("-f",'--visualization', required=False, help="reference file")
requiredGrp.add_argument("-m",'--mode_of_analysis', required=True, help="Isoform or Exome")
requiredGrp.add_argument("-i",'--additional_information', required=False, help="output file location",action="store",nargs = "*")
args = vars(ap.parse_args())



def str2bool(v):
#Taken from https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1','True'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0','False'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

transcriptome_file = args['transcriptome_file']
path_transcripts = args['list']
test_file = args['test_file']
control_file = args['control_file']

log2fc = args['log2fc']
odds = args['odds']
padj = args['padj']
m6A_status = args['m6A_status']
fraction_diff = args['fraction_diff']
visualization_input = args['visualization']
output_dir = args['output_dir']
mode_of_analysis = args['mode_of_analysis']
additional_columns = args['additional_information']

if log2fc == None:
    log2fc = 0.5
    
if odds == None:
    odds = 1.5
    
if padj == None:
    padj = 0.05

if fraction_diff == None:
    fraction_diff = 0.01

if m6A_status == None:
    m6A_status = True
else:
	print('m6A_status',m6A_status)
	m6A_status = str2bool(m6A_status)
	
if visualization_input == None:
    visualization_input = False
else: 
	print('Visualization_input',visualization_input)
	visualization_input = str2bool(visualization_input)


all_permutations = [(x,y) for x in control_file for y in test_file]
all_permutations_dict = {'rep_'+str(i+1):all_permutations[i] for i in range(len(all_permutations))}




def which(program):
#Taken from https://github.com/pinellolab/CRISPResso2/blob/master/CRISPResso2/CRISPRessoPooledCORE.py
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def check_samtools():
#Taken from https://github.com/pinellolab/CRISPResso2/blob/master/CRISPResso2/CRISPRessoPooledCORE.py
    cmd_path=which('samtools')
    if cmd_path:
        #print('SAMTOOLS IN PATH')
        return True
    else:
        sys.stdout.write('\nDRUMMER requires samtools')
        sys.stdout.write('\n\nPlease install samtools and add it to your path following the instructions at: http://www.htslib.org/download/')
        return False
check_samtools()



devnull = open(os.devnull, 'w')

transcript_directory = output_dir + '/transcripts/'

def isoform_mode(transcriptome_file):
	with open(transcriptome_file, 'r') as reader:
		all_transcripts = [l.rstrip("\n") for l in reader]
	transcripts_dict = dict([all_transcripts[i:i+2] for i in range(0,len(all_transcripts),2)])

	length_dictionary = {}
	for transcript,sequence in transcripts_dict.items():
		file_name = transcript.strip('>')
		length_dictionary[file_name] = len(sequence)
		os.makedirs(transcript_directory, exist_ok=True) 
		f = open(transcript_directory + file_name + '.fa', "w")
		write_file = [transcript+'\n',sequence]
		f.writelines(write_file)
		f.close()

	list_of_transcripts_df = pd.read_csv(path_transcripts,sep = '\t',header = None)
	shape_lst = list_of_transcripts_df.shape[-1]
	if shape_lst == 1:
		list_of_transcripts_df.columns = ['name']
	else:
		list_of_transcripts_df.columns = ['chrom','name','strand','thickStart','blockCount','blockSizes','blockStarts']
	return list_of_transcripts_df,length_dictionary,shape_lst

def exome_mode():
	pass

if mode_of_analysis == 'Isoform':
	iterate,length_dictionary,shape_lst = isoform_mode(transcriptome_file)
else: 
	iterate = all_permutations_dict.keys()
	
	
bam_readcount = "/Users/mac/Desktop/DRUMMER/modules/bam-readcount"
from odds import run_odds
from motif import run_motif

pooled_string = colored('''
 ________________________________
| __   __                 __  __ |
||  \ |_/ |  | |\/| |\/| |_  |_/ |
||__/ | \ |__| |  | |  | |__ | \ |
|________________________________|
''','green')
print(pooled_string)

def do_work(i,rep,comp,comp2):

	current_file = transcript_directory + i + '.fa'
	length = length_dictionary[i]
# 	print('printing X',rep)
# 	print('printing Y1',comp)
# 	print('printing Y2',comp2)
	map_directory = output_dir+'/' +rep+ '/map/'
	os.makedirs(map_directory, exist_ok=True)
# 	print(i)
	 
	view_test_format = "samtools view -b {} {}:1-{} -o {}/{}/map/{}.UNMOD.bam".format(comp,i,length,output_dir,rep,i).split(' ')
	sort_test_format = "samtools sort -o {}/{}/map/{}.UNMOD.sorted.bam {}/{}/map/{}.UNMOD.bam".format(output_dir,rep,i,output_dir,rep,i).split(' ')
	index_test_format = "samtools index {}/{}/map/{}.UNMOD.sorted.bam".format(output_dir,rep,i).split(' ')
# 	print(view_test_format)
	subprocess.call(view_test_format)
	subprocess.call(sort_test_format)
	subprocess.call(index_test_format)
	
	view_ctrl_format = "samtools view -b {} {}:1-{} -o {}/{}/map/{}.MOD.bam".format(comp2,i,length,output_dir,rep,i).split(' ')
	sort_ctrl_format = "samtools sort -o {}/{}/map/{}.MOD.sorted.bam {}/{}/map/{}.MOD.bam".format(output_dir,rep,i,output_dir,rep,i).split(' ')
	index_ctrl_format = "samtools index {}/{}/map/{}.MOD.sorted.bam".format(output_dir,rep,i).split(' ')

	subprocess.call(view_ctrl_format)
	subprocess.call(sort_ctrl_format)
	subprocess.call(index_ctrl_format)

	bam_readcount_dir = output_dir +'/'+rep+ '/bam_readcount/'

	os.makedirs(bam_readcount_dir,exist_ok = True)
	
	bam_readcount_unmod =  ".././modules/bam-readcount -f {}/transcripts/{}.fa {}/{}/map/{}.UNMOD.sorted.bam -w 1".format(output_dir,i,output_dir,rep,i).split(' ')
	bam_readcount_mod =  ".././modules/bam-readcount -f {}/transcripts/{}.fa {}/{}/map/{}.MOD.sorted.bam -w 1".format(output_dir,i,output_dir,rep,i).split(' ')
	
	with open(bam_readcount_dir+i+'.UNMOD.bamreadcount.txt','w') as fout: #context manager is OK since `call` blocks :)
		subprocess.call(bam_readcount_unmod,stdout=fout, stderr=subprocess.DEVNULL)

	with open(bam_readcount_dir+i+'.MOD.bamreadcount.txt','w') as fout: #context manager is OK since `call` blocks :)
		subprocess.call(bam_readcount_mod,stdout=fout, stderr=subprocess.DEVNULL)

	merged_df = merge.merged_dataframes(bam_readcount_dir+i+'.MOD.bamreadcount.txt',bam_readcount_dir+i+'.UNMOD.bamreadcount.txt',i)
    
	os.makedirs(output_dir+'/'+rep+'/odds/',exist_ok = True)
	odds_df = run_odds(merged_df)
	odds_df.to_csv(output_dir+'/'+rep+'/odds/' +i+'.txt',sep = '\t',index =False)
	from motif import run_motif
	os.makedirs(output_dir+'/'+rep+'/motif/',exist_ok = True)
	motif_df = motif.run_motif(odds_df,m6A_status)
	motif_df.to_csv(output_dir+'/'+rep+'/motif/' +i+'.txt',sep = '\t',index =False)

	from gtest import run_gtest
	gtest_df = gtest.run_gtest(motif_df)
	os.makedirs(output_dir+'/'+rep+'/gTest/',exist_ok = True)
	gtest_df.to_csv(output_dir+'/'+rep+'/gTest/' +i+'.txt',sep = '\t',index =False)
	
	from candidates import run_find_candidates
	candidates_df = candidates.run_find_candidates(gtest_df,odds,padj,fraction_diff)
	os.makedirs(output_dir+'/'+rep+'/complete_analysis/',exist_ok = True)
	candidates_df.to_csv(output_dir+'/'+rep+'/complete_analysis/' +i+'.complete.txt',sep = '\t',index =False)


	#os.makedirs(output_dir+'/'+rep+'/complete_analysis/',exist_ok = True) shape_lst
	if shape_lst != 1:
		from genomic import run_genomics
		genomics_df = genomic.run_genomics(output_dir+'/'+rep+'/complete_analysis/' +i+'.complete.txt',path_transcripts)
		genomics_df.to_csv(output_dir+'/'+rep + '/complete_analysis/'+i+'.complete.txt',sep = '\t',index = False)

# 	print(genomics_df)
	from visualization import run_visualization
# 	print('visualization_input',visualization_input,type(visualization_input))
	if visualization_input == True:
		os.makedirs(output_dir+'/'+rep+'/visualization/', exist_ok = True)
		run_visualization(genomics_df,m6A_status,output_dir+'/'+rep + '/visualization/'+i+'.pdf',mode_of_analysis,i)
	#print('Completed {} sample {}'.format(rep,i))
from itertools import repeat
from summary import run_summary
from m6a_summary import run_plotting


def progressbar(it, prefix="", size=60, file=sys.stdout):
    count = len(it)
    def show(j):
        x = int(size*j/count)
        file.write("%s[%s%s] %i/%i\r" % (prefix, "#"*x, "."*(size-x), j, count))
        file.flush()        
    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)
    file.write("\n")
    file.flush()
    

from tqdm.contrib.concurrent import process_map
from p_tqdm import p_map
for replicate,comparisons in all_permutations_dict.items():
    try:
        l = process_map(do_work, iterate.name,repeat(replicate),repeat(comparisons[0]),repeat(comparisons[1]), max_workers=16,file=sys.stdout,miniters=1,desc = replicate,colour = 'red')

        #print('Working on {} --> {} -- {} \n...'.format(replicate,comparisons[0],comparisons[1]))
        #with concurrent.futures.ProcessPoolExecutor() as executor:
        #    for _ in tqdm.tqdm(executor.map(do_work, iterate.name,repeat(replicate),repeat(comparisons[0]),repeat(comparisons[1])),total = len(iterate.name),file=sys.stdout,miniters=1,desc = replicate):
        #        pass
    ##Replace True with m6A status and additional_columns
    except UnboundLocalError:
        continue
for repl in all_permutations_dict.keys():
    summary_df = run_summary(output_dir+'/'+repl + '/complete_analysis/',m6A_status,additional_columns)
    summary_df.to_csv(output_dir+'/'+repl + '/summary.txt',sep = '\t')
    if m6A_status == True:
        run_plotting(output_dir+'/'+repl + '/complete_analysis/',summary_df,output_dir+'/'+repl + '/m6A_plot.pdf')
        


# for replicate,comparisons in all_permutations_dict.items():
#     with concurrent.futures.ProcessPoolExecutor() as executor:
#         executor.map(do_work, iterate.name,repeat(replicate),repeat(comparisons[0]),repeat(comparisons[1]))
# ##Replace True with m6A status and additional_columns
#     summary_df = run_summary(output_dir+'/'+replicate + '/complete_analysis/',m6A_status,additional_columns)
#     summary_df.to_csv(output_dir+'/'+replicate + '/summary.txt',sep = '\t')
#     if m6A_status == True:
#         run_plotting(output_dir+'/'+replicate + '/complete_analysis/',summary_df,output_dir+'/'+replicate + '/m6A_plot.pdf')











