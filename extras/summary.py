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


my_parser = argparse.ArgumentParser()
my_parser.add_argument('-i','--input', action='store',help = "Directory containing DRUMMER outputs")
my_parser.add_argument("-c", "--columns", action="store",nargs = "*",help="Additional columns to add to summary file")
my_parser.add_argument('-o','--output', action='store', help="Output file ")
my_parser.add_argument("-m", "--mode", action="store",help="Mode of analyis")

args = my_parser.parse_args()

target_dir = args.input
additional_columns = args.columns
output_dir = args.output
mode = args.mode

cwd = os.getcwd()
onlydir = ['/'+f for f in listdir(target_dir) if isfile(join(target_dir, f))]

all_candidates = pd.DataFrame()
for i in onlydir:
    full_path = target_dir + i
    df = pd.read_csv(full_path,sep = '\t')
    df = df[df['candidate_site'] == 'candidate']
    all_candidates = pd.concat([all_candidates,df])
all_candidates = all_candidates.reset_index(drop = True)
if mode == 'True':
	keep_columns = ['chr_mod','Chromosome','ref_mod','pos_mod','depth_mod','ref_fraction_mod','depth_unmod','ref_fraction_unmod','frac_diff','odds_ratio','p_values_OR_adj','eleven_bp_motif','G_test','padj','candidate_site','nearest_ac','nearest_ac_motif','genomic_position']
else:
	keep_columns = ['chr_mod','Chromosome','ref_mod','pos_mod','depth_mod','ref_fraction_mod','depth_unmod','ref_fraction_unmod','frac_diff','odds_ratio','p_values_OR_adj','eleven_bp_motif','G_test','padj','candidate_site','genomic_position']

if 'Chromosome' not in all_candidates.columns:
    	keep_columns.remove('Chromosome')
if 'genomic_position' not in all_candidates.columns:
    	keep_columns.remove('genomic_position')
    	
if additional_columns != None:
	keep_columns += additional_columns
	

final_candidates = all_candidates[keep_columns]

final_candidates['ref_fraction_mod'] = final_candidates['ref_fraction_mod'].round(3)
final_candidates['ref_fraction_unmod'] = final_candidates['ref_fraction_unmod'].round(3)
final_candidates['odds_ratio'] = final_candidates['odds_ratio'].round(3)
final_candidates['p_values_OR_adj'] = final_candidates['p_values_OR_adj'].astype(float).map('{:0.2e}'.format)
final_candidates['ref_fraction_unmod'] = final_candidates['ref_fraction_unmod'].round(3)
final_candidates['padj'] = final_candidates['padj'].astype(float).map('{:0.2e}'.format)
final_candidates['frac_diff'] = final_candidates['frac_diff'].round(3)

final_candidates=final_candidates.rename(columns = {'p_values_OR_adj':'OR_padj','padj':'G_padj','chr_mod':"transcript_id",'ref_mod':'reference_base'})

#final_candidates = final_candidates.sort_values('chr_mod',ascending = True)
final_candidates.to_csv(output_dir, sep ='\t',index = False)
