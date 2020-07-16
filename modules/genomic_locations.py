import pandas as pd
import numpy as np
import argparse
import math
import re

ap = argparse.ArgumentParser(description = 'Takes in the output from the pipeline and determines candidate sites \
using log2fc, odds_ratio and padj')
requiredGrp = ap.add_argument_group('required arguments')
requiredGrp.add_argument("-t",'--transcript_file', required=True, help="transcript file containing genomic information")
requiredGrp.add_argument("-i",'--input', required=True, help="input csv file containing candidate sites")
requiredGrp.add_argument("-o",'--output', required=True, help="input csv file containing candidate sites")

args = vars(ap.parse_args())
transcripts = args['transcript_file']
path = args['input']
output = args['output']
include_candidate_df = pd.read_csv(path)
transcripts = pd.read_csv(transcripts,sep = '\t',header= None)
candidate_positions = include_candidate_df[include_candidate_df['candidate_site'] == 'Candidate']['pos_unmod'].tolist()

# print('transcripts',transcripts)
# transcript_of_interest = transcripts[transcripts['transcript'] == chromo_name]



def sum_exons(df,strand):
    '''Get the sum of all exons in transcript file'''
    number_exons = df['num_exons']
    get_exons = df['length_exon'].tolist()[0].split(',')
    get_exons = list(map(int, get_exons))
    if strand == '-':
        get_exons = list(reversed(get_exons))
    exon_length = [0]
    running_sum = 0
    for i in get_exons:
        running_sum += i
        exon_length.append(running_sum)
    return exon_length
    
def get_exon_spans(exon_list:list):
    '''Returns a list of spans indicating the start and end of each exon'''
    exon_dic = {}
    exon_dic = {exon:span for exon,span in enumerate(zip(exon_list,exon_list[1:]))}
#     print('get_exon_spans',exon_dic)
#         if index > 0 and index < len(lst_of_exons):
#             lst_of_exon_coordinates = list(range(lst_of_exons[index-1]+1,lst_of_exons[index]+1))
# #             print(lst_of_exon_coordinates)
#             exon_dic[index] = [lst_of_exon_coordinates[0],lst_of_exon_coordinates[-1]]
    return exon_dic
    
def find_exon(exon_dic:dict,candidate_sites:list):
    '''Finds which exon the candidate site is in
    '''
    exon_number = {}
    for key,values in exon_dic.items():
        exon_candidate = []
        for candidate in candidate_sites:
            if values[0] <= candidate <= values[-1]:
                exon_candidate.append(candidate)
        exon_number[key] = exon_candidate
    return exon_number
    
def summing_values(df,exon_num,location,strand,candidate):
    '''Sums all relevant exons to get final location
    '''
    genomic_start = df.iloc[0]['genome_start']
    exon_starts = list(map(int,df.iloc[0]['start_exon'].split(',')))
    exon_length = list(map(int,df.iloc[0]['length_exon'].split(',')))
#     print('candidate',candidate)
#     print('genomic start',genomic_start)
#     print('exon_starts',exon_starts)
#     print('exon_num',exon_num)
#     print('location',location)
#     print('length',exon_length)
    if strand == '+':
        #print('summing',sum(exon_length[:exon_num-1]))
#         print('exon_starts[exon_num]', exon_num)
        summing_values = location - sum(exon_length[:exon_num]) + exon_starts[exon_num] + genomic_start
#         print(summing_values)
#         summing_values = genomic_start + sum(exon_length[:exon_num-1]) + sum(exon_starts[:exon_num]) + location
    else:
        [exon_starts,exon_ends] = adjust_negative(genomic_start,exon_starts,exon_length)
        reverse_lengths = list(reversed(exon_length))
#         print('specific exon',exon_ends[exon_num])
#         print('location',location)
#         print('reverse_lengths',reverse_lengths)
#         print('test_output',exon_ends[exon_num] - location - sum(reverse_lengths[:exon_num]))
        summing_values = exon_ends[exon_num] -  (location-sum(reverse_lengths[:exon_num]) + exon_num)
#         exon_num = len(exon_length) - exon_num
#         summing_values = genomic_start + sum(exon_length[:exon_num]) + sum(exon_starts[:exon_num+1]) + location 
    return summing_values
    
def adjust_negative(g_starts,e_starts,e_length):
    e_starts = list(reversed(e_starts))
    e_length = list(reversed(e_length))
    new_starts = []
    new_ends = []
    for i in range(len(e_starts)):
        new_ends.append(g_starts + e_starts[i] + e_length[i])
        new_starts.append(g_starts + e_starts[i]+1)
    return [new_starts,new_ends]


args = vars(ap.parse_args())
transcripts = args['transcript_file']
path = args['input']
output = args['output']
include_candidate_df = pd.read_csv(path)
transcripts = pd.read_csv(transcripts,sep = '\t',header= None)
candidate_positions = include_candidate_df[include_candidate_df['candidate_site'] == 'Candidate']['pos_unmod'].tolist()

if len(transcripts.columns) > 1:
	col_names = ['transcript','strand','genome_start','num_exons','length_exon','start_exon']
	transcripts.columns = col_names

	#Find relevant transcript using chromosome column of dataframe
	if len(set(include_candidate_df['chr_unmod'])) == 1:
		chromo_name = include_candidate_df['chr_unmod'][0]
	else:
		print('multiple chromosome names in csv file')

	print('chromo_name',chromo_name)
	transcript_of_interest = transcripts[transcripts['transcript'] == chromo_name].reset_index(drop=True)
	print(transcript_of_interest)
	strand = transcript_of_interest.loc[0]['strand']

	lst_of_exons = sum_exons(transcript_of_interest,strand)
	# candidate_positions = [724,754]
	lst_of_exons = sorted(list(set(lst_of_exons)))
	# print('list of exons',lst_of_exons)
	# print(lst_of_exons)
	exon_spans = get_exon_spans(lst_of_exons)
	# print(exon_spans)

	finding_exons = find_exon(exon_spans,candidate_positions)

	final_coordinates = {}
	updated_coordinates = []
	for key,values in finding_exons.items():
		summed_vals = [summing_values(transcript_of_interest,key,location,strand,candidate_positions) for location in values]
		final_coordinates[key] = summed_vals
		updated_coordinates.append(summed_vals)

	updated_coordinates = [location for locations in updated_coordinates for location in locations]

	include_candidate_df['genomic_position'] = ' '
	candidate_and_updated = list(zip(candidate_positions,updated_coordinates))
	for cand,updt in candidate_and_updated:
		include_candidate_df['genomic_position'].loc[include_candidate_df['pos_unmod'] == cand] = updt
	include_candidate_df.to_csv(output,sep = '\t')

    
    
    
    
