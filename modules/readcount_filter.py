import os
from os import listdir
from os.path import isfile, join

import pandas as pd
import numpy as np
import argparse
from create_output_file import create_output
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

columns_names = ['chr','pos','ref','depth','A','C','G','T','N']

ap = argparse.ArgumentParser(description = 'Takes in the output of bam-readcount \
and returns a text file containing the count of each nucleotide at the position and the \
fraction of the reference nucleotide among all reads.')
requiredGrp = ap.add_argument_group('required arguments')
requiredGrp.add_argument("-i","--input", required=True, help="input file location")
requiredGrp.add_argument("-o","--output", required=True, help="output directory location")

args = vars(ap.parse_args())
input = args['input']
output = args['output']

# Deletions subtract from the N column.
# Every Negative N column value is equal to zero.
# A flag option to disregard homopolymer indels > 3. (default is to include)


# print('output',output)
def proper_filter(read_metrics:list):
    """Takes in a list returns a dictionary with counts of each nucleotide
    """
    nucleotides = ['A','C','G','T','N']
    keep_track = {k:0 for k in nucleotides} #Initialize each nucleotide count with 0
    for indiv_metrics in read_metrics:
        if indiv_metrics[0] in nucleotides: 
            splitted_metrics = indiv_metrics.split(':')
            keep_track[splitted_metrics[0]] += int(splitted_metrics[1])
        elif indiv_metrics[0] not in nucleotides and indiv_metrics[0] != '=':
            splitted_metrics = indiv_metrics.split(':')
#             print('splitted_metrics',splitted_metrics)
#             print('This is the sign', splitted_metrics[0][0])
            if splitted_metrics[0][0] == '+':
            	keep_track['N'] += int(splitted_metrics[1])
            else:
            	keep_track['N'] -= int(splitted_metrics[1])
#             	print('splitted_metrics',splitted_metrics)
#             	print(splitted_metrics[1])
#             if len(splitted_metrics[0]) <= 3: #Keep in count the -2 - +2  indels
#             	keep_track['N'] += int(splitted_metrics[1]) 
#     print('keep track N', keep_track)
    if keep_track['N'] < 0:
    	keep_track['N'] = 0
#     print(keep_track)
    return keep_track
    
def new_depth(df):
    new_depths = []
    for index,column in df.iterrows():
        new_depths.append(sum([column['A'],column['C'],column['G'],column['T'],column['N']]))
    return new_depths
    
def do_math(df:'DataFrame',index:int):
    """Takes in dataframe and row index, returns that rows reference fraction
    """
    reference_nucleotide = df.loc[index,'ref']
    depth = df.loc[index,'depth']
    if int(depth) == 0:
        return 0
    else:
        numerator = df.loc[index,reference_nucleotide]
        return int(numerator)/int(depth)
    
#argument = 'bamreadcount/E3.RIDb'.split('/')
input = input.split('/')
# print('input',input)
# print(input)
file = input.pop(-1)
mypath = '/'.join(input)
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
# print('File',file)
#Test file is always second
k = [mypath+'/'+i for i in onlyfiles if i.startswith(file)]
k = sorted(k, key=lambda x:('TEST' in x, x))

print('k',k)
try:
	lst = []
	for i in k:
		list_of_rows = []
		with open(i) as FileObj:
			for lines in FileObj:
				initial_split = lines.split('\t')
				first_4_dict = dict(zip([1,2,3,4],initial_split[:4])) #First four columns into dict with random keynames
				nucleotide_count = proper_filter(initial_split[4:])
				first_4_dict.update(nucleotide_count)
				list_of_rows.append(pd.DataFrame([first_4_dict]))
		filtered_df = pd.concat(list_of_rows).reset_index(drop=True)
		filtered_df.columns = columns_names
		filtered_df['depth'] = new_depth(filtered_df)
#		filtered_df = filtered_df[filtered_df['depth'] != 0].reset_index(drop=True)
		filtered_df['ref_fraction'] = [do_math(filtered_df,i) for i in range(len(filtered_df))]
		lst.append(filtered_df)
	#     print(i)
	#     print('file',file)
		output_path = create_output(output,i,'filtered',file)
	#     print('output',output)
		without_sub = create_output(output,i,'filtered')
	#     print('without_sub',without_sub)
	#     print(# output)
		filtered_df.to_csv(output_path,sep = '\t', index = False)

	cols = ['chr', 'pos', 'ref', 'depth', 'A', 'C', 'G', 'T', 'N', 'ref_fraction',
		   'chr.1', 'pos.1', 'ref.1', 'depth.1', 'A.1', 'C.1', 'G.1', 'T.1', 'N.1',
		   'ref_fraction.1']
	   
	# merged_df = pd.concat([lst[0],lst[1]],axis =1)
	# merged_df.columns = cols
	merged_df = lst[0].merge(lst[1],on = 'pos',suffixes = ('','.1'))
	merged_dir = output + '/' + 'merged/'
	# print('merged_dir',merged_dir)
	os.makedirs(merged_dir, exist_ok = True)

	output_dir = merged_dir + file + '.merged.txt'
	# print('output_dir',output_dir)
	#merged_dir = without_sub.split('.')
	#merged_dir.insert(-1,'merged')
	#merged_dir = '.'.join(merged_dir)

	merged_df.to_csv(output_dir,sep = '\t', index = False)
except ValueError:
	print("Problem merging:",file)









