import pandas as pd
import numpy as np
import argparse
import math
import re
from create_output_file import create_output
warnings.filterwarnings("ignore")

ap = argparse.ArgumentParser(description = 'Takes in the output from the pipeline and determines candidate sites \
using log2fc, odds_ratio and padj')
requiredGrp = ap.add_argument_group('required arguments')
requiredGrp.add_argument("-i",'--input', required=True, help="input file location")
requiredGrp.add_argument("-r",'--odds_ratio', required=True, help="odds ratio cutoff")
#requiredGrp.add_argument("-l",'--log2_fc', required=False, help="input file location")
requiredGrp.add_argument("-p",'--padj', required=True, help="padj cutoff")
requiredGrp.add_argument("-d",'--fraction_diff', required=True, help="fraction difference")
requiredGrp.add_argument("-o",'--output', required=True, help="output file location")



args = vars(ap.parse_args())
input = args['input']
odds_ratio = float(args['odds_ratio'])
#log2fc = float(args['log2_fc'])
padj = float(args['padj'])
output = args['output']
fraction_diff = args['fraction_diff']
#print("FRACTION_DIFF",fraction_diff)
# print(log2fc, odds_ratio, padj)
def is_candidate(df,odds_ratio,padj):
    '''Takes in a dataframe looks at log2fc, odds and padj to determine if the site is a candidate'''
#     print('log2fc',df[df['log2_fc']<log2fc])
    idx = ((df['odds_ratio']>odds_ratio) &(df['padj']< padj) & ( (df['ref_fraction_unmod'] - df['ref_fraction_mod'] ) > float(fraction_diff)) & (df['p_values_OR_adj']<padj))
	#idx = (df['log2_fc']>log2fc) & (df['odds_ratio']>odds_ratio) &(df['padj']< padj)
    df['candidate_site'] = ['candidate' if i == True else '' for i in idx ]
    return df

def Diff(li1, li2): 
    return (list(set(li1) - set(li2))) 
include_candidate_df = pd.read_csv(input,sep = '\t')
include_candidate_df = is_candidate(include_candidate_df,odds_ratio,padj)

index_candidates = list(include_candidate_df[include_candidate_df['candidate_site'] == 'candidate'].index)

if len(index_candidates) > 1:
    #print('Found {} candidate sites'.format(len(index_candidates)))
    full_list = []
    lst = []
    for k,value in enumerate(index_candidates):
        if k > 0:
            if index_candidates[k] - index_candidates[k-1] < 5:
                if index_candidates[k] not in lst:
                    lst.append(index_candidates[k])
            else:
                full_list.append(lst)
                lst = []
                if index_candidates[k] not in lst:
                    lst.append(index_candidates[k])
        else:
            lst.append(index_candidates[k])
    full_list.append(lst)
    index_of_highest = []
    for i in full_list:
        k = []
        for j in i:
            k.append(include_candidate_df.iloc[j]['G_test'])
        index_of_highest.append(i[np.argmax(k)])
    
    include_candidate_df['candidate_site'] = ''

    include_candidate_df.loc[index_of_highest,'candidate_site'] = 'candidate'

    flattened_list = [indx for lsts in full_list for indx in lsts]

    candidate_masks = Diff(flattened_list,index_of_highest)

    include_candidate_df.loc[candidate_masks,'candidate_site'] = '[candidate_masked]'



def check_homopolymer(string):
    return_val = []
    for i in ['TTT','AAA','GGG','CCC']:
        if i in string:
            return True
include_candidate_df['homopolymer'] = [check_homopolymer(i) for i in include_candidate_df['eleven_bp_motif']]
#include_candidate_df['candidate_site'].value_counts()
#ide_candidate_df['frac_diff'] = include_candidate_df['ref_fraction_unmod'] - include_candidate_df['ref_fraction_mod']
include_candidate_df.insert(20, 'frac_diff', include_candidate_df['ref_fraction_mod'] - include_candidate_df['ref_fraction_unmod'])
# print('Candidate sites:\n',include_candidate_df['candidate_site'].value_counts())
# include_candidate_df.to_csv(output,index=False)

#output = create_output(input,input,'candidates')
# print(input)
# make_dir = output_location = output +'/all_transcripts_complete/'
# os.makedirs(make_dir, exist_ok = True)
# output_location = output +'/visualization/'+ sample + '.pdf'
print(output)
include_candidate_df.to_csv(output,sep = '\t', index = False)
