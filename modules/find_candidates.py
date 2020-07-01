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
ap.add_argument("-r","--odds_ratio", required=False,nargs='?',const=1.3,type=int,default=1.3, help="odds ratio cutoff")
ap.add_argument("-l","--log2_fc", required=False,nargs='?', const=-0.5,type=int,default=-0.5,help="Log2 fold change cutoff")
ap.add_argument("-p","--padj", required=False,nargs='?',const=0.05,type=int,default=0.05, help="padj cutoff")




args = vars(ap.parse_args())
input = args['input']
odds_ratio = args['odds_ratio']
log2fc = args['log2_fc']
padj = args['padj']

# print(log2fc, odds_ratio, padj)
def is_candidate(df,log2fc,odds_ratio,padj):
    '''Takes in a dataframe looks at log2fc, odds and padj to determine if the site is a candidate'''
#     print('log2fc',df[df['log2_fc']<log2fc])
    idx = (df['log2_fc']>log2fc) & (df['odds_ratio']>odds_ratio) &(df['padj']< padj)
    df['candidate_site'] = ['Candidate' if i == True else '' for i in idx ]
    return df

include_candidate_df = pd.read_csv(input,sep = '\t')
include_candidate_df = is_candidate(include_candidate_df,log2fc,odds_ratio,padj)
# print(include_candidate_df.head())

print('Candidate sites:\n',include_candidate_df['candidate_site'].value_counts())
# include_candidate_df.to_csv(output,index=False)

output = create_output(input,'candidates')
include_candidate_df.to_csv(output,sep = '\t', index = False)