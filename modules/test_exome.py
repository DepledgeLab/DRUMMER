from __future__ import division
import argparse
import modules.support 
from modules.odds import run_odds
from modules.motif import run_motif
import modules.merge
from modules.gtest import run_gtest
from modules.candidates import run_find_candidates
from modules.genomic import run_genomics
from modules.visualization import run_visualization
from itertools import repeat
from modules.summary import run_summary
import modules.candidates
from modules.m6a_summary import run_plotting
#from tqdm.contrib.concurrent import process_map
#from p_tqdm import p_map
import sys
import math
import numpy as np
import os 
import pandas as pd
import subprocess
import modules.genomic
import concurrent.futures
import modules.multiple_combine
from modules.plot_multi import run_main
from Bio import SeqIO
from modules.parsing_check import exome_mode_parsing
devnull = open(os.devnull, 'w')
def exome_mode(transcriptome_file,transcript_directory):
    length_dictionary = {}
    f_open = open(transcriptome_file,'rU')
    for rec in SeqIO.parse(f_open,'fasta'):
        id = rec.id
        seq = rec.seq
        out_file = transcript_directory + id + '.fa'
        id_file = open(out_file,'w')
        length_dictionary[id] = len(seq)
        id_file.write(">"+str(id)+"\n"+str(seq))
        id_file.close()
    return length_dictionary
#     with open(transcriptome_file, 'r') as reader:
#         all_transcripts = [l.rstrip("\n") for l in reader]
#     file_name = all_transcripts[0].strip('>')
#     sequence = ''.join(all_transcripts[1:])
#     output = '\n'.join(all_transcripts)
# #     print('SEQ',sequence)
#     length_dictionary = len(sequence)
# #     print('all_transcripts',all_transcripts)
#     f = open(transcript_directory + file_name + '.fa','w') 
#     f.writelines(output)
#     f.close()
    # return length_dictionary

def run_samtools(comp,comp2,transcript_id,length,output_dir,rep):
	view_test_format = "samtools view -b {} {}:1-{} -o {}/{}/map/{}.UNMOD.bam".format(comp,transcript_id,length,output_dir,rep,transcript_id).split(' ')
	bedtools_unmod = "bedtools genomecov -d -split -ibam {}".format(comp).split(' ')
	sort_test_format = "samtools sort -o {}/{}/map/{}.UNMOD.sorted.bam {}/{}/map/{}.UNMOD.bam".format(output_dir,rep,transcript_id,output_dir,rep,transcript_id).split(' ')
	index_test_format = "samtools index {}/{}/map/{}.UNMOD.sorted.bam".format(output_dir,rep,transcript_id).split(' ')

	subprocess.call(view_test_format)
	subprocess.call(sort_test_format)
	subprocess.call(index_test_format)
	
	view_ctrl_format = "samtools view -b {} {}:1-{} -o {}/{}/map/{}.MOD.bam".format(comp2,transcript_id,length,output_dir,rep,transcript_id).split(' ')
	bedtools_mod = "bedtools genomecov -d -split -ibam {}".format(comp2).split(' ')
	sort_ctrl_format = "samtools sort -o {}/{}/map/{}.MOD.sorted.bam {}/{}/map/{}.MOD.bam".format(output_dir,rep,transcript_id,output_dir,rep,transcript_id).split(' ')
	index_ctrl_format = "samtools index {}/{}/map/{}.MOD.sorted.bam".format(output_dir,rep,transcript_id).split(' ')
	subprocess.call(view_ctrl_format)
	subprocess.call(sort_ctrl_format)
	subprocess.call(index_ctrl_format)
# 	print(comp)
	print(output_dir+'/'+rep +'/map/' + transcript_id + '.unmod.genomecov.txt')
	with open(output_dir+'/'+rep +'/map/' + transcript_id + '.unmod.genomecov.txt','w') as fout:
	    subprocess.call(bedtools_unmod,stdout=fout)
	with open(output_dir+'/'+rep +'/map/' + transcript_id + '.mod.genomecov.txt','w') as fout:
	    subprocess.call(bedtools_mod,stdout=fout)
	#subprocess.call(bedtools_mod)
	#subprocess.call(bedtools_unmod)

    
def run_bamreadcounts(output_dir,transcript_id,rep,length):
	bam_readcount_unmod =  "modules/bam-readcount -q 0 -b 0 -d 1000000 -w 1 -f {}/transcript/{}.fa {}/{}/map/{}.UNMOD.sorted.bam {}:1-{}".format(output_dir,transcript_id,output_dir,rep,transcript_id,transcript_id,length).split(' ')
	bam_readcount_mod =  "modules/bam-readcount -q 0 -b 0 -d 1000000 -w 1 -f {}/transcript/{}.fa {}/{}/map/{}.MOD.sorted.bam {}:1-{}".format(output_dir,transcript_id,output_dir,rep,transcript_id,transcript_id,length).split(' ')
# 	print('OUTPUT_DIR',output_dir)
# 	print('transcript_id',transcript_id)
# 	print('rep',rep)
# 	print('bam_read_countDIR',bam_readcount_dir)
	with open(bam_readcount_dir+transcript_id+'.UNMOD.bamreadcount.txt','w') as fout: #context manager is OK since `call` blocks :)
		subprocess.call(bam_readcount_unmod,stdout=fout)
	with open(bam_readcount_dir+transcript_id+'.MOD.bamreadcount.txt','w') as fout: #context manager is OK since `call` blocks :)
		subprocess.call(bam_readcount_mod,stdout=fout, stderr=subprocess.DEVNULL)
		
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
    
def split_rep(lst):
    new_lst = []
    for i in lst:
        fixed_names = []
        for j in i:
            if '/' in j:
                file_name = j.split('/')[-1]
                fixed_names.append(file_name.strip('.sorted.bam'))
            else:
                return_name = j.strip('.sorted.bam')
                fixed_names.append(return_name)
        new_lst.append('-'.join(fixed_names))
    return new_lst

def do_work(all_permuations_w_replicates,i,transcript_directory,output_dir,m6A_status,odds,padj,fraction_diff,mode_of_analysis,deletion_filter):
	rep = all_permuations_w_replicates[0]
	comp = all_permuations_w_replicates[-1][1] # control
	comp2 = all_permuations_w_replicates[-1][0] # treatment
# 	print('AALLL PERMUAT',all_permuations_w_replicates)
# 	print('AALLL PERMUAT rep',all_permuations_w_replicates[0])
# 	print('AALLL PERMUAT comp',all_permuations_w_replicates[-1][0])
# 	print('AALLL PERMUAT comp2',all_permuations_w_replicates[-1][1])
	#os.makedirs(output_dir+'/'+rep+'/map/',exist_ok = True)

	
	current_file = transcript_directory + i + '.fa'
	print('Current_file',current_file)
	os.makedirs(output_dir+'/'+rep+'/map/',exist_ok = True)
	os.makedirs(output_dir+'/'+rep+'/complete_analysis/',exist_ok = True)
	os.makedirs(output_dir+'/'+rep+'/gTEST/',exist_ok = True)
	global bam_readcount_dir
	bam_readcount_dir = output_dir +'/' +rep+ '/bam_readcount/'
	os.makedirs(bam_readcount_dir,exist_ok = True)
	run_samtools(comp,comp2,i,length,output_dir,rep)
	run_bamreadcounts(output_dir,i,rep,length)
	os.makedirs(output_dir+'/'+rep+'/MERGED/',exist_ok = True)
	os.makedirs(output_dir+'/'+rep+'/ODDS/',exist_ok = True)
	os.makedirs(output_dir+'/'+rep+'/MOTIF/',exist_ok = True)
	#os.makedirs(output_dir+'/'+rep+'/gTEST/',exist_ok = True)
	#os.makedirs(output_dir+'/'+rep+'/complete_analysis/',exist_ok = True)
	merged_df = modules.merge.merged_dataframes(bam_readcount_dir+i+'.MOD.bamreadcount.txt',bam_readcount_dir+i+'.UNMOD.bamreadcount.txt',i,deletion_filter)
	merged_df.to_csv(output_dir+'/'+rep+'/MERGED/' +i+'.txt',sep = '\t',index =False)
	print('MERGED')
	odds_df = run_odds(merged_df)
	odds_df.to_csv(output_dir+'/'+rep+'/ODDS/' +i+'.txt',sep = '\t',index =False)
	print('ODDS')
	motif_df = run_motif(odds_df,m6A_status)
	motif_df.to_csv(output_dir+'/'+rep+'/MOTIF/' +i+'.txt',sep = '\t',index =False)
	print('MOTIF')
	#gtest_df = run_gtest(motif_df)
	#gtest_df.to_csv(output_dir+'/'+rep+'/gTEST/' +i+'.txt',sep = '\t',index =False)
	gtest_df = run_gtest(motif_df)
	gtest_df.to_csv(output_dir+'/'+rep+'/gTEST/' +i+'.txt',sep = '\t',index =False)
	candidates_df = modules.candidates.run_find_candidates(output_dir+'/'+rep+'/gTEST/' +i+'.txt',odds,padj,fraction_diff)
	#print(candidates_df)
	#os.makedirs(output_dir+'/'+rep+'/complete_analysis/',exist_ok = True)
	candidates_df.to_csv(output_dir+'/'+rep+'/complete_analysis/' +i+'.complete.txt',sep = '\t',index =False)
	os.makedirs(output_dir+'/'+rep+'/complete_analysis/',exist_ok = True)
	candidates_df.to_csv(output_dir+'/'+rep+'/complete_analysis/' +i+'.complete.txt',sep = '\t',index =False)
	
	os.makedirs(output_dir+'/'+rep+'/visualization/', exist_ok = True)
# 	run_visualization(candidates_df,m6A_status,output_dir+'/'+rep + '/visualization/'+i+'.pdf',mode_of_analysis,i)



def main(transcriptome_file,test_file,name,control_file,odds,padj,m6A_status,fraction_diff,output_dir,mode,deletion_filter):
    transcript_directory = output_dir + '/transcript/'
    global length
    global shape_lst
    os.makedirs(output_dir+'/'+'/transcript/',exist_ok = True)
    all_permutations = [(x,y) for x in control_file for y in test_file]
    replicate_names = split_rep(all_permutations)
    #replicate_names = ['rep1','rep2','rep3','rep4']
    all_permutations_w_rep = list(zip(replicate_names,all_permutations))
    length = exome_mode(transcriptome_file,transcript_directory)[name]
    print('all_permutations_w_rep',all_permutations_w_rep)
    print('****LENGTH****',length)
    print(all_permutations)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(do_work, all_permutations_w_rep,repeat(name),repeat(transcript_directory),repeat(output_dir),repeat(m6A_status),repeat(odds),repeat(padj),repeat(fraction_diff),repeat('exome'),repeat(deletion_filter))
    #l = process_map(do_work, all_permutations_w_rep,repeat(name),repeat(transcript_directory),repeat(output_dir),repeat(m6A_status),repeat(odds),repeat(padj),repeat(fraction_diff),repeat('exome'))
#    l = process_map(do_work, iterate.name,repeat(replicate),repeat(comparisons[0]),repeat(comparisons[1]),repeat(transcript_directory),repeat(output_dir),repeat(m6A_status),
#    repeat(odds),repeat(padj),repeat(fraction_diff),repeat(path_transcripts),repeat(visualization_input),repeat(mode), max_workers=16,file=sys.stdout,miniters=1,desc = replicate)

    for repl in replicate_names:
        summary_df = run_summary(output_dir+'/'+repl + '/complete_analysis/',m6A_status)
        summary_df.to_csv(output_dir+'/'+repl + '/summary.txt',sep = '\t',index = False)
        exome_mode_parsing(output_dir+'/'+repl,name)
        if m6A_status == True:
            run_plotting(output_dir+'/'+repl + '/complete_analysis/',summary_df,output_dir+'/'+repl + '/m6A_plot.pdf')
    if len(replicate_names) > 1:
        final_df = modules.multiple_combine.main(output_dir,replicate_names,m6A_status)
        final_df.to_csv(output_dir +'/'+'multiple_comp.txt',sep = '\t',na_rep='NaN',index= False)
        run_main(final_df,output_dir)
    
    #python3 FAST-DRUMMER.py -r /gpfs/data/depledgelab/Jonathan/Evaluation/transcripts/Adenovirus-Ad5.fasta -n Ad5 -o exome-test -c /gpfs/home/ja3539/depledge/Jonathan/DRUMMER/TESTDATA/exome.Ad5.MOD.bam  /gpfs/home/ja3539/depledge/Jonathan/DRUMMER/TESTDATA/exome2.Ad5.MOD.bam -t /gpfs/home/ja3539/depledge/Jonathan/DRUMMER/TESTDATA/exome.Ad5.UNMOD.bam /gpfs/home/ja3539/depledge/Jonathan/DRUMMER/TESTDATA/exome2.Ad5.UNMOD.bam  -a exome -m True
    #python3 FAST-DRUMMER.py -r /gpfs/data/depledgelab/Jonathan/Evaluation/transcripts/Ad5_v9.1_complete.fasta -l /gpfs/home/ja3539/depledge/Jonathan/DRUMMER/TESTDATA/Ad5.sample.transcripts.txt -o isoform-test -c /gpfs/home/ja3539/depledge/Jonathan/DRUMMER/TESTDATA/isoform.Ad5.MOD.bam /gpfs/home/ja3539/depledge/Jonathan/DRUMMER/TESTDATA/isoform2.Ad5.MOD.bam -t /gpfs/home/ja3539/depledge/Jonathan/DRUMMER/TESTDATA/isoform.Ad5.UNMOD.bam /gpfs/home/ja3539/depledge/Jonathan/DRUMMER/TESTDATA/isoform2.Ad5.UNMOD.bam  -a isoform -m True
if __name__ == "__main__":
    main('/Users/mac/Desktop/DRUMMER/TESTDATA/Adenovirus-Ad5.fasta',['/Users/mac/Desktop/DRUMMER/TESTDATA/exome.Ad5.MOD.bam'],
    'Ad5',['/Users/mac/Desktop/DRUMMER/TESTDATA/exome.Ad5.UNMOD.bam'],1,1,True,1,'test-exome-runs','exome')