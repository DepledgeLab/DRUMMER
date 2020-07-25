import pandas as pd
import numpy as np
import argparse
import math
import re
from create_output_file import create_output

ap = argparse.ArgumentParser(description = 'Takes in the output from the pipeline and determines candidate sites \
using log2fc, odds_ratio and padj')
requiredGrp = ap.add_argument_group('required arguments')
requiredGrp.add_argument("-i",'--input', required=True, help="input file location")
requiredGrp.add_argument("-r",'--odds_ratio', required=True, help="input file location")
#requiredGrp.add_argument("-l",'--log2_fc', required=False, help="input file location")
requiredGrp.add_argument("-p",'--padj', required=True, help="input file location")
requiredGrp.add_argument("-o",'--output', required=True, help="output file location")



args = vars(ap.parse_args())
input = args['input']
odds_ratio = float(args['odds_ratio'])
#log2fc = float(args['log2_fc'])
padj = float(args['padj'])
output = args['output']

# print(log2fc, odds_ratio, padj)
def is_candidate(df,odds_ratio,padj):
    '''Takes in a dataframe looks at log2fc, odds and padj to determine if the site is a candidate'''
#     print('log2fc',df[df['log2_fc']<log2fc])
    idx = ((df['odds_ratio']>odds_ratio) &(df['padj']< padj) & (df['ref_fraction_mod'] < df['ref_fraction_unmod']))
	#idx = (df['log2_fc']>log2fc) & (df['odds_ratio']>odds_ratio) &(df['padj']< padj)
    df['candidate_site'] = ['candidate' if i == True else '' for i in idx ]
    return df

def Diff(li1, li2): 
    return (list(set(li1) - set(li2))) 
include_candidate_df = pd.read_csv(input,sep = '\t')
include_candidate_df = is_candidate(include_candidate_df,odds_ratio,padj)

index_candidates = list(include_candidate_df[include_candidate_df['candidate_site'] == 'candidate'].index)

if len(index_candidates) > 2:
	final_list = [[index_candidates.pop(0)]]

	for ind in index_candidates:
		if ind - final_list[-1][0] < 5:
			final_list[-1].append(ind)
		else:
			final_list.append([ind])
		
	index_of_highest = []
	for i in final_list:
		k = []
		for j in i:
			k.append(include_candidate_df.iloc[j]['G_test'])
		index_of_highest.append(i[np.argmax(k)])
	
	include_candidate_df['candidate_site'] = ''

	include_candidate_df.loc[index_of_highest,'candidate_site'] = 'candidate'

	flattened_list = [indx for lsts in final_list for indx in lsts]

	candidate_masks = Diff(flattened_list,index_of_highest)

	include_candidate_df.loc[candidate_masks,'candidate_site'] = '[candidate_masked]'
else:
	print('Found {} candidate sites'.format(len(index_candidates)))
	
#include_candidate_df['candidate_site'].value_counts()

print('Candidate sites:\n',include_candidate_df['candidate_site'].value_counts())
# include_candidate_df.to_csv(output,index=False)

# output = create_output(input,input,'candidates')
print(input)
include_candidate_df.to_csv(output,sep = '\t', index = False)
