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
from os import listdir
from os.path import isfile, join
import pandas as pd
import subprocess
import modules.genomic
import concurrent.futures
import modules.multiple_combine
from modules.plot_multi import run_main
from Bio import SeqIO
from modules.parsing_check import isoform_mode_parsing
from functools import reduce

def isoform_mode(transcriptome_file,transcript_directory,path_transcripts):
    f_open = open(transcriptome_file, "rU")
    length_dictionary = {}
    for rec in SeqIO.parse(f_open, "fasta"):
        id = rec.id
        seq = rec.seq
        length_dictionary[id] = len(seq)

        #print(id,len(seq))
        #        f = open(transcript_directory + file_name + '.fa', "w")
        out_file = transcript_directory + id + '.fa'
        id_file = open(out_file, "w")
        id_file.write(">"+str(id)+"\n"+str(seq)) 
        id_file.close()
    list_of_transcripts_df = pd.read_csv(path_transcripts,sep = '\t',header = None)
    shape_lst = list_of_transcripts_df.shape[-1]
    if shape_lst == 1:
        list_of_transcripts_df.columns = ['name']
    else:
        list_of_transcripts_df.columns = ['chrom','name','strand','thickStart','blockCount','blockSizes','blockStarts']
    return list_of_transcripts_df,length_dictionary,shape_lst

def run_samtools(comp,comp2,transcript_id,length,output_dir,rep):
	view_test_format = "samtools view -b {} {}:1-{} -o {}/{}/map/{}.UNMOD.bam".format(comp,transcript_id,length,output_dir,rep,transcript_id).split(' ')
	#bedtools_unmod = "bedtools genomecov -d -split -ibam {}".format(comp).split(' ')
	sort_test_format = "samtools sort -o {}/{}/map/{}.UNMOD.sorted.bam {}/{}/map/{}.UNMOD.bam".format(output_dir,rep,transcript_id,output_dir,rep,transcript_id).split(' ')
	index_test_format = "samtools index {}/{}/map/{}.UNMOD.sorted.bam".format(output_dir,rep,transcript_id).split(' ')

	subprocess.call(view_test_format)
	subprocess.call(sort_test_format)
	subprocess.call(index_test_format)
	
	view_ctrl_format = "samtools view -b {} {}:1-{} -o {}/{}/map/{}.MOD.bam".format(comp2,transcript_id,length,output_dir,rep,transcript_id).split(' ')
	#bedtools_mod = "bedtools genomecov -d -split -ibam {}".format(comp2).split(' ')
	sort_ctrl_format = "samtools sort -o {}/{}/map/{}.MOD.sorted.bam {}/{}/map/{}.MOD.bam".format(output_dir,rep,transcript_id,output_dir,rep,transcript_id).split(' ')
	index_ctrl_format = "samtools index {}/{}/map/{}.MOD.sorted.bam".format(output_dir,rep,transcript_id).split(' ')
	subprocess.call(view_ctrl_format)
	subprocess.call(sort_ctrl_format)
	subprocess.call(index_ctrl_format)

def run_bedtools(comp,comp2,output_dir,rep):
	bedtools_unmod = "bedtools genomecov -d -split -ibam {}".format(comp).split(' ')
	bedtools_mod = "bedtools genomecov -d -split -ibam {}".format(comp2).split(' ')
	with open(output_dir+'/'+rep +'/map/'  + 'FULL.unmod.genomecov.txt','w') as fout:
	    subprocess.call(bedtools_unmod,stdout=fout)
	with open(output_dir+'/'+rep +'/map/' + 'FULL.mod.genomecov.txt','w') as fout:
	    subprocess.call(bedtools_mod,stdout=fout)   
    
    
    
def run_bamreadcounts(output_dir,transcript_id,rep,length):
	bam_readcount_unmod =  "modules/bam-readcount -q 0 -b 0 -d 1000000 -w 1 -f {}/transcripts/{}.fa {}/{}/map/{}.UNMOD.sorted.bam {}:1-{}".format(output_dir,transcript_id,output_dir,rep,transcript_id,transcript_id,length).split(' ')
	bam_readcount_mod =  "modules/bam-readcount -q 0 -b 0 -d 1000000 -w 1 -f {}/transcripts/{}.fa {}/{}/map/{}.MOD.sorted.bam {}:1-{}".format(output_dir,transcript_id,output_dir,rep,transcript_id,transcript_id,length).split(' ')
	print(bam_readcount_unmod)
	print(bam_readcount_mod)
	with open(bam_readcount_dir+transcript_id+'.UNMOD.bamreadcount.txt','w') as fout: #context manager is OK since `call` blocks :)
		subprocess.call(bam_readcount_unmod,stdout=fout, stderr=subprocess.DEVNULL)
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

def do_work(i,rep,comp2,comp,transcript_directory,output_dir,m6A_status,odds,padj,fraction_diff,path_transcripts,visualization_input,mode_of_analysis,deletion_filter):
	#print('started..{}'.format(i))
	current_file = transcript_directory + i + '.fa'
	length = length_dictionary[i]
	global bam_readcount_dir
	bam_readcount_dir = output_dir +'/'+rep+ '/bam_readcount/'
	os.makedirs(bam_readcount_dir,exist_ok = True)
	os.makedirs(output_dir+'/'+rep+'/MERGED/',exist_ok = True)
	os.makedirs(output_dir+'/'+rep+'/ODDS/',exist_ok = True)
	os.makedirs(output_dir+'/'+rep+'/MOTIF/',exist_ok = True)
	run_samtools(comp,comp2,i,length,output_dir,rep)
	run_bamreadcounts(output_dir,i,rep,length)
	merged_df = modules.merge.merged_dataframes(bam_readcount_dir+i+'.MOD.bamreadcount.txt',bam_readcount_dir+i+'.UNMOD.bamreadcount.txt',i,deletion_filter)
# 	print(merged_df.head(5))
	merged_df.to_csv(output_dir+'/'+rep+'/MERGED/' +i+'.txt',sep = '\t',index =False)
	odds_df = run_odds(merged_df)
	motif_df = run_motif(odds_df,m6A_status)
	gtest_df = run_gtest(motif_df)
	os.makedirs(output_dir+'/'+rep+'/gTEST/',exist_ok = True)
	gtest_df.to_csv(output_dir+'/'+rep+'/gTEST/' +i+'.txt',sep = '\t',index =False)
	candidates_df = modules.candidates.run_find_candidates(output_dir+'/'+rep+'/gTEST/' +i+'.txt',odds,padj,fraction_diff)
	#print(candidates_df)
	os.makedirs(output_dir+'/'+rep+'/complete_analysis/',exist_ok = True)
	candidates_df.to_csv(output_dir+'/'+rep+'/complete_analysis/' +i+'.complete.txt',sep = '\t',index =False)
	#print('shape_lst',shape_lst)
	if shape_lst != 1:
		#print('IN GENOMICS')
		genomics_df = genomic.run_genomics(output_dir+'/'+rep+'/complete_analysis/' +i+'.complete.txt',path_transcripts)
		genomics_df.to_csv(output_dir+'/'+rep + '/complete_analysis/'+i+'.complete.txt',sep = '\t',index = False)
	#print('visualization_input',visualization_input)

	if visualization_input == True:
		#print('IN visualizaation')
		os.makedirs(output_dir+'/'+rep+'/visualization/', exist_ok = True)
		run_visualization(candidates_df,m6A_status,output_dir+'/'+rep + '/visualization/'+i+'.pdf',mode_of_analysis,i)
	print('completed..{}'.format(i))

def return_top3_transcripts(path):
    transcript_max = {}
    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    for i in onlyfiles:
        current_df = pd.read_csv(path+i,sep = '\t',usecols = ['depth_ctrl','depth_treat'])
        current_df['mean'] = current_df.mean(axis = 1)
        max_mean = max(current_df['mean'])
        transcript_max[i] = max_mean
    top_3_trans = sorted(transcript_max, key=transcript_max.get, reverse=True)[:3]
    print(top_3_trans)
    top_3_trans = list(map(lambda x:x.replace('.complete.txt',''),top_3_trans))
    return top_3_trans

def main(transcriptome_file,test_file,path_transcripts,control_file,odds,padj,m6A_status,fraction_diff,visualization_input,output_dir,mode,deletion_filter):
    transcript_directory = output_dir + '/transcripts/'
    global length_dictionary
    global shape_lst
    os.makedirs(output_dir+'/'+'/transcripts/',exist_ok = True)
    #os.makedirs(output_dir+'/'+rep+'/map/',exist_ok = True)
    all_permutations = [(x,y) for x in control_file for y in test_file]
    replicate_names = split_rep(all_permutations)
    all_permutations_w_rep = list(zip(replicate_names,all_permutations))
    #all_permutations_dict = {'rep_'+str(i+1):all_permutations[i] for i in range(len(all_permutations))}
    iterate,length_dictionary,shape_lst = isoform_mode(transcriptome_file,transcript_directory,path_transcripts)
    #print('iterate',iterate)
    print(all_permutations_w_rep)
    [os.makedirs(output_dir+'/'+output_loc+'/map/',exist_ok = True) for output_loc,compare in all_permutations_w_rep]
    [run_bedtools(compare[0],compare[1],output_dir,output_loc) for output_loc,compare in all_permutations_w_rep]
    for replicate,comparisons in all_permutations_w_rep:
        try:
            with concurrent.futures.ProcessPoolExecutor() as executor:
                executor.map(do_work, iterate.name,repeat(replicate),repeat(comparisons[0]),repeat(comparisons[1]),repeat(transcript_directory),repeat(output_dir),repeat(m6A_status),
                repeat(odds),repeat(padj),repeat(fraction_diff),repeat(path_transcripts),repeat(visualization_input),repeat(mode),repeat(deletion_filter))
        except UnboundLocalError:
            continue

    for repl in replicate_names:
        summary_df = run_summary(output_dir+'/'+repl + '/complete_analysis/',m6A_status)
        summary_df.to_csv(output_dir+'/'+repl + '/summary.txt',sep = '\t',index = False)
        top_3 = return_top3_transcripts(output_dir+'/'+repl + '/complete_analysis/')
        isoform_mode_parsing(output_dir+'/'+repl,top_3)
        if m6A_status == True:
            run_plotting(output_dir+'/'+repl + '/complete_analysis/',summary_df,output_dir+'/'+repl + '/m6A_plot.pdf')
    if len(replicate_names) > 1:
        final_df = modules.multiple_combine.main(output_dir,replicate_names,m6A_status)
        if shape_lst != 1:
            mapping = dict(zip(iterate['name'],iterate['strand']))
            final_df.insert(3, 'strand', final_df["transcript_id"].map(mapping))
        final_df.to_csv(output_dir +'/'+'multiple_comp.txt',sep = '\t',na_rep='NaN',index= False)
        run_main(final_df,output_dir)
    
if __name__ == "__main__":
    main('/gpfs/home/ja3539/depledge/Jonathan/DRUMMER/TESTDATA/Ad5_v9.1_complete.fasta',['/gpfs/data/depledgelab/Jonathan/Adeno-DRUMMER/ALIGNMENTS/Ad5_transcript.M3P1-40.sorted.bam'],
    '/gpfs/home/ja3539/depledge/Jonathan/DRUMMER/TESTDATA/Ad5.sample.transcripts.txt',['/gpfs/data/depledgelab/Jonathan/Adeno-DRUMMER/ALIGNMENTS/Ad5_transcript.M3KO1-40.sorted.bam'],
    1,1,1,True,1,True,'test-runs-full',None,'isoform')