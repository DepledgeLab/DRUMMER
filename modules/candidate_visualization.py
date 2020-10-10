import pandas as pd
import os 
import argparse
from os import listdir
from os.path import isdir, join, isfile
import glob
import itertools
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import math

ap = argparse.ArgumentParser(description = 'Takes in the output from the pipeline and determines candidate sites \
using log2fc, odds_ratio and padj')
requiredGrp = ap.add_argument_group('required arguments')
requiredGrp.add_argument("-i",'--input', required=True, help="input file location")
requiredGrp.add_argument("-o",'--output', required=True, help="output file location")

args = vars(ap.parse_args())
input = args['input']
output = args['output']

sample = input.split('/')[-1].strip('.complete.txt')

df = pd.read_csv(input,sep = '\t')
#df = df[['chr_mod','pos_mod','candidate_site','G_test','padj']]
candidate_count = df['candidate_site'].value_counts().to_dict()
# print(df.head())

def return_shape_size(df,color):
    legend_dict = defaultdict(int)
    col_g = []
    size = []
    count = 0
    for index,row in df.iterrows():
        if row['candidate_site'] == 'candidate' and row['homopolymer'] == True:
            count +=0
            col_g.append('blue')
            size.append(50)
            legend_dict['homopolymer'] += 1
        elif row['candidate_site'] == 'candidate' and row['homopolymer'] != True:
            count += 0
            col_g.append('red')
            size.append(50)
            legend_dict['no_homopolymer'] += 1
        elif row['candidate_site'] == '[candidate_masked]':
            col_g.append('red')
            size.append(10)
            legend_dict['masked'] +=1
        else:
            col_g.append('black')
            size.append(1)
            legend_dict['none'] += 1
    return [col_g,size,legend_dict,count]
    
    
fig, ax1 = plt.subplots(figsize=(20, 10))

ax1.set_title('Visualization of Candidate and Masked sites ({})'.format(sample),fontsize = 36)
ax1.set_xlabel('Position',fontsize = 20)
ax1.set_ylabel('gTEST',fontsize = 20)
plt.xticks(np.arange(0, max(df['pos_mod'])+50, 100),rotation = 45)
col_g,size,legend_dict,count = return_shape_size(df,'red')
# print(legend_dict)
scatter1 = ax1.scatter(list(df['pos_mod']), list(df['G_test']), c=col_g,s = size)

# ax1.legend((scatter1,scatter1,scatter1,scatter1), ('Candidate & Significant = {}'.format(legend_dict['candidate_sign'])
#                                            , 'Only significant = {}'.format(legend_dict['sign']),
#                                            "None = {}".format(legend_dict['none']),
#                                         'masked = {}'.format(legend_dict['masked'])),
#                                            loc='upper left', bbox_to_anchor=(0, 1))

ax1.legend((scatter1,scatter1,scatter1,scatter1), ('Candidate & homopolymer = {}'.format(legend_dict['homopolymer']),
                                                   'Candidate & no homopolymer = {}'.format(legend_dict['no_homopolymer']),
                                           "Candidate Masked = {}".format(legend_dict['masked']),
                                        'none = {}'.format(legend_dict['none'])),
                                           loc='upper left', bbox_to_anchor=(0, 1))
leg = ax1.get_legend()
leg.legendHandles[0]._sizes = [50]
leg.legendHandles[1]._sizes = [50]
leg.legendHandles[2]._sizes = [10]
leg.legendHandles[3]._sizes = [1]
# leg.legendHandles[3]._sizes = [1]
leg.legendHandles[0].set_color('blue')
leg.legendHandles[1].set_color('red')
leg.legendHandles[2].set_color('red')
leg.legendHandles[3].set_color('black')
# leg.legendHandles[3].set_color('grey')

scaling_x = (int(math.ceil(df.shape[0] / 10.0)) * 10)/100
scaling_y = math.ceil(max(df.G_test)/100)
for index,value in df.iterrows():
    if value['candidate_site'] == 'candidate':
        ax1.annotate(value.pos_mod,(value.pos_mod-scaling_x,value.G_test+scaling_y),size = 7.5)
#         ax1.annotate(value.five_bp_motif,(value.pos_mod-25,value.G_test+50),size = 7.5)
# plt.savefig('/Users/mac/Desktop/DRUMMER_Figures/Fix_masked_issues/homopolymer-US-1-candidate-site_visualization.pdf')\
make_dir = output_location = output +'/visualization/'
os.makedirs(make_dir, exist_ok = True)
output_location = output +'/visualization/'+ sample + '.pdf'
plt.savefig(output_location)