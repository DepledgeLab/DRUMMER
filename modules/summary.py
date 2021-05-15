import pandas as pd
import argparse
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import math
import os
from os import listdir
from os.path import isdir, join, isfile
import warnings
warnings.filterwarnings("ignore")


def run_summary(target_dir,mode):
    onlydir = ['/'+f for f in listdir(target_dir) if isfile(join(target_dir, f))]
    #df[((df["Name"]=="Tom") & (df["Age"]<=42)) | (df["Age"]<=34)]
    all_candidates = pd.DataFrame()
    for i in onlydir:
        full_path = target_dir + i
        df = pd.read_csv(full_path,sep = '\t')
        #df = df[df['candidate_site'] == 'candidate']
        #df = df[(df['candidate_site'] == 'candidate') | (df['candidate_site'] == '[candidate_masked]')]
        df = df[(df['accumulation'] == 'accumulation') | (df['depletion'] == 'depletion')]
        all_candidates = pd.concat([all_candidates,df])
    all_candidates = all_candidates.reset_index(drop = True)
    if mode == True:
        keep_columns = ['chr_ctrl','Chromosome','ref_ctrl','pos_ctrl','genomic_position','depth_ctrl','depth_treat','ref_fraction_ctrl','ref_fraction_treat','frac_diff','odds_ratio','log2_(OR)','p_values_OR_adj','eleven_bp_motif','G_test','padj','accumulation','depletion','nearest_ac','nearest_ac_motif','homopolymer','is_SNP']
        
    else:
        keep_columns = ['chr_ctrl','Chromosome','ref_ctrl','pos_ctrl','genomic_position','depth_ctrl','depth_treat','ref_fraction_ctrl','ref_fraction_treat','frac_diff','odds_ratio','log2_(OR)','p_values_OR_adj','eleven_bp_motif','G_test','padj','accumulation','depletion','homopolymer','is_SNP']

    if 'Chromosome' not in all_candidates.columns:
        keep_columns.remove('Chromosome')
    if 'genomic_position' not in all_candidates.columns:
        keep_columns.remove('genomic_position')

#     if additional_columns != None:
#         keep_columns += additional_columns


    final_candidates = all_candidates[keep_columns]

    final_candidates['ref_fraction_treat'] = final_candidates['ref_fraction_treat'].round(3)
    final_candidates['ref_fraction_ctrl'] = final_candidates['ref_fraction_ctrl'].round(3)
    final_candidates['odds_ratio'] = final_candidates['odds_ratio'].round(3)
    final_candidates['p_values_OR_adj'] = final_candidates['p_values_OR_adj'].astype(float).map('{:0.2e}'.format)
    #final_candidates['ref_fraction_unmod'] = final_candidates['ref_fraction_unmod'].round(3)
    final_candidates['padj'] = final_candidates['padj'].astype(float).map('{:0.2e}'.format)
    final_candidates['frac_diff'] = final_candidates['frac_diff'].round(3)
    final_candidates=final_candidates.rename(columns = {'p_values_OR_adj':'OR_padj','padj':'G_padj','chr_ctrl':"transcript_id",'ref_ctrl':'reference_base','pos_ctrl':'transcript_pos'})
    return final_candidates
	#final_candidates = final_candidates.sort_values('chr_mod',ascending = True)
	#final_candidates.to_csv(output_dir, sep ='\t',index = False)
	
if __name__ == "__main__":
	ap = argparse.ArgumentParser(description = 'Takes in the output from the pipeline and determines candidate sites using log2fc, odds_ratio and padj')

	requiredGrp = ap.add_argument_group('required arguments')
	requiredGrp.add_argument("-i",'--Nanocompore_input', required=True, help="Input nanocompore output")
	requiredGrp.add_argument("-o",'--output', required=True, help="output location (.txt extension)")
	args = vars(ap.parse_args())

	input_path = args['Nanocompore_input'] 
	output_dir = args['output']  
	df = run_summary(input_path,True)
	df.to_csv(output_dir,sep = '\t',index = None)