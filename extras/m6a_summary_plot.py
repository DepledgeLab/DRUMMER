import pandas as pd
import argparse
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import math
import os
from os import listdir
from os.path import isdir, join, isfile
import random
import statistics
import pandas as pd
from scipy import stats
from statsmodels.stats import weightstats as stests
import statistics
from scipy.stats import bartlett
from scipy.stats import levene
import warnings
warnings.filterwarnings("ignore")

ap = argparse.ArgumentParser(description = 'Takes in the output from the pipeline and determines candidate sites \
using log2fc, odds_ratio and padj')
requiredGrp = ap.add_argument_group('required arguments')
requiredGrp.add_argument("-d",'--directory', required=True, help="directory containing all files")
requiredGrp.add_argument("-i",'--input', required=True, help="input file location")
requiredGrp.add_argument("-o",'--output', required=True, help="output file location")

args = vars(ap.parse_args())
only_candidates_path = args['input']
target_dir = args['directory']
output = args['output']

random.seed(30)
#target_dir = '/Users/mac/Desktop/DRUMMER_Figures/transcript_level_sites/P1-KO1-36_update/'
onlydir = ['/'+f for f in listdir(target_dir) if isfile(join(target_dir, f))]

candidates_df = pd.read_csv(only_candidates_path,sep = '\t')

if len(candidates_df) > 50:
	all_sites = pd.DataFrame()
	for i in onlydir:
		full_path = target_dir + i
		df = pd.read_csv(full_path,sep = '\t')
	#     df = df[df['candidate_site'] == 'candidate']
		all_sites = pd.concat([all_sites,df])
	all_sites = all_sites.reset_index(drop = True)

	#only_candidates_path = '/Users/mac/Desktop/DRUMMER_Figures/transcript_level_sites/ac_distance_text.txt'


	get_shape_candidates = candidates_df.shape[0]
	get_shape_random = all_sites.shape[0]

	randomlist = random.sample(range(0, get_shape_random), get_shape_candidates)
	random_df = all_sites.iloc[randomlist]

	candiates_distance = candidates_df['nearest_ac'].value_counts()
	random_distance = random_df['nearest_ac'].value_counts()

	stat, p = levene(candidates_df['nearest_ac'], random_df['nearest_ac'])
	_,pval = stats.ks_2samp(candidates_df['nearest_ac'], random_df['nearest_ac'])

	pval = '{:0.3e}'.format(pval)
	p = '{:0.3e}'.format(p)

	lst_counts = [candiates_distance,random_distance]
	colors = ['darkred','slategrey']
	type_line = ['solid','solid']
	labels = ['Candidates','Random']
	alpha = [1,.5]
	edge = ['black','black']
	total_sums = [' ({})'.format(sum(candiates_distance)),'']


	fig, ax = plt.subplots(figsize=(20, 10))
	for i in range(len(lst_counts)):
		current_df = lst_counts[i].sort_index()
		x = current_df.index
		y = list(current_df)
	#     new_x,new_y = filtered_for_10(x,y)
		ax.bar(x, y, linestyle=type_line[i], color=colors[i], lw=2.5,
			alpha=alpha[i], label=labels[i]+total_sums[i],edgecolor = edge[i])
		plt.xlim(-25, 25) 
	ax.legend(fontsize = 15)
	ax.set_xlabel('AC distance from Candidate Site',fontsize = 20)
	ax.set_ylabel('Count of each Candidate Site calls',fontsize = 20)
	ax.set_title('Distance of each called candidate site',fontsize = 36)
	ax.tick_params(labelsize=15)
	top_limit = plt.ylim()[-1]

	ax.text(-25, max(candiates_distance*.95), ('  Kolmogorov–Smirnov test p-value: {}').format(pval), fontsize=14,fontweight='bold')
	ax.text(-25, max(candiates_distance*.85), ('  Levene test p-value: {}').format(p), fontsize=14,fontweight='bold')

	plt.savefig(output)













