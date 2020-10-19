import pandas as pd
import argparse
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import math

my_parser = argparse.ArgumentParser()
my_parser.add_argument('--inputs', action='store', nargs="*")
my_parser.add_argument('--output', action='store', help="output file (.pdf))")
my_parser.add_argument("-p", "--positional", action="store", default=True,
                    help="positional labeling")
my_parser.add_argument("-hp", "--homopolymer", action="store", default=True,
                    help="homopolymer based coloring")
my_parser.add_argument("-a", "--m6A_mode", action="store", default=True,
                    help="arrow pointing to AC motifs")
my_parser.add_argument("-l", "--sample_labels", action="store",nargs = "*",
                    help="sample labels must equal length of inputs")
my_parser.add_argument('-r','--range', action='store',help = "start and end range for plotting", nargs=2)
my_parser.add_argument("-s", "--line_space", action="store", default=100,
                    help="linespace for plotting")
args = my_parser.parse_args()

print('linespace',args.line_space)
# print(args.inputs)
# print(args.output)
# print(args.positional)
# print(args.homopolymer)
# print(args.sample_labels)
linespace= int(args.line_space)
if args.sample_labels == None:
	names = [i for i in range(len(args.inputs))]
else:
	names = args.sample_labels
# print('names',names)
lst_of_paths = args.inputs
homopolymer = args.homopolymer
pos_label = args.positional
m6A = args.m6A_mode

def roundup(x):
    return int(math.ceil(x / 10.0)) * 10
def filter_df(df,start,end):
    df = df[df['pos_mod'] > start]
    df = df[df['pos_mod'] < end]
    return df

all_max = []
if args.range != None and len(args.range) ==2:
	for i in lst_of_paths:
		df = filter_df(pd.read_csv(i,sep = '\t'),int(args.range[0]),int(args.range[-1]))
		current_max = roundup(max(df.G_test))
		all_max.append(current_max)
	start = int(args.range[0])
	end = int(args.range[-1])
	max_shape_arrow = max(all_max)
	max_shape = max_shape_arrow + 10
else:
	max_shape_arrow = roundup(max([max(pd.read_csv(i,sep = '\t').G_test) for i in lst_of_paths]))
	max_shape = max_shape_arrow + 10
	start = 0

def filter_df(df,start,end):
    df = df[df['pos_mod'] > start]
    df = df[df['pos_mod'] < end]
    return df

def return_shape_size(df,color,homopolymer):
    legend_dict = defaultdict(int)
    col_g = []
    size = []
    count = 0
    for index,row in df.iterrows():
        if row['candidate_site'] == 'candidate' and row['homopolymer'] == True:
            count +=0
            size.append(100)
            if homopolymer == True:
                col_g.append('blue')
            else:
                col_g.append('red')
        elif row['candidate_site'] == 'candidate' and row['homopolymer'] != True:
            count += 0
            col_g.append('red')
            size.append(100)
            legend_dict['no_homopolymer'] += 1
        elif row['candidate_site'] == '[candidate_masked]' and row['homopolymer'] == True: 
            size.append(10)
            legend_dict['masked'] +=1
            if homopolymer == True:
                col_g.append('blue')
            else:
                col_g.append('red')
        elif row['candidate_site'] == '[candidate_masked]' and row['homopolymer'] != True: 
            col_g.append('red')
            size.append(10)
            legend_dict['masked'] +=1
        else:
            col_g.append('black')
            size.append(10)
            legend_dict['none'] += 1
    return [col_g,size,legend_dict,count]
    
fig, ax1 = plt.subplots(figsize=(50, 20),sharey= True)
row = len(lst_of_paths)
col = 1
cnt = 1
max_value = []
for i in range(1,row+1):
    df = pd.read_csv(lst_of_paths[i-1],sep = '\t')
    if args.range != None and len(args.range) == 2:
    	df = filter_df(df,start,end)
    plt.subplot(row,col,cnt)
#     plt.xticks([], [])
    plt.ylabel('{} \n G-Test score'.format(names[i-1]),size = 25)
    col_g,size,legend_dict,count = return_shape_size(df,'red',homopolymer)
    plt.scatter(list(df['pos_mod']), list(df['G_test']), c=col_g,s = size)
    plt.yticks(np.arange(0, max_shape, step=10),size = 20)
    if pos_label == True:
        for index,value in df.iterrows():
            if value['candidate_site'] == 'candidate':
                plt.annotate(value.pos_mod,(value.pos_mod+5,value.G_test-5),size = 12)
    fig.tight_layout()
    cnt += 1
    if cnt - 1 == row:
    	plt.xticks(np.arange(start, max(df['pos_mod']), linespace),size=40)
    else:
    	plt.xticks([], [])
plt.subplot(row,col,1)
if m6A == True:
    for i,j in df.iterrows():
        if j['five_bp_motif'][2:4] == 'AC':
            plt.annotate('', xy=(j['pos_mod'], max_shape_arrow),xycoords='data',xytext=(j['pos_mod'], max_shape),
					textcoords='data',
					arrowprops=dict(arrowstyle= '-|>',color='slategrey',lw=2.5,ls='--'))
plt.subplot(row,col,row)
plt.xlabel('Transcript Position',size = 40)
fig.tight_layout()
plt.savefig(args.output)