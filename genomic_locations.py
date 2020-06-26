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

args = vars(ap.parse_args())
transcripts = args['transcript_file']
path = args['input']

include_candidate_df = pd.read_csv(path)
transcripts = pd.read_csv(transcripts,sep = '\t',header= None)
candidate_positions = include_candidate_df[include_candidate_df['candidate_site'] == 'Candidate']['pos_unmod'].tolist()

# print('transcripts',transcripts)
# transcript_of_interest = transcripts[transcripts['transcript'] == chromo_name]

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


def sum_exons(df):
    '''Get the sum of all exons in transcript file'''
    number_exons = df['num_exons']
    get_exons = df['length_exon'].tolist()[0].split(',')
    get_exons = list(map(int, get_exons)) 
    exon_length = [0]
    running_sum = 0
    for i in get_exons:
        running_sum += i
        exon_length.append(running_sum)
    return exon_length
    
def get_exon_spans(exon_list:list):
    '''Returns a list of spans indicating the start and end of each exon'''
    exon_dic = {}
    for index in range(len(lst_of_exons)+1):
        if index > 0 and index < len(lst_of_exons):
            lst_of_exon_coordinates = list(range(lst_of_exons[index-1]+1,lst_of_exons[index]+1))
#             print(lst_of_exon_coordinates)
            exon_dic[index] = [lst_of_exon_coordinates[0],lst_of_exon_coordinates[-1]]
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
    
    
def new_locations(df,exon_number,strand):
    '''Determines the new location of candidate site within corresponding exon'''
    new_exon_locations = {}
    exon_lengths = exon_length = list(map(int,df.iloc[0]['length_exon'].split(',')))
    for exon,locations in exon_number.items():
        if strand == '+':
            adding_prev_exons = sum(exon_lengths[:exon-1])
        elif strand == '-':
            adding_prev_exons = sum(exon_length[exon-1:])
#             adding_prev_exons = 
        new_locations = [location-adding_prev_exons for location in locations]
        new_exon_locations[exon] = new_locations
    return new_exon_locations
    
def summing_values(df,exon_num,location,strand):
    '''Sums all relevant exons to get final location
    '''
    genomic_start = df.iloc[0]['genome_start']
    exon_starts = list(map(int,df.iloc[0]['start_exon'].split(',')))
    exon_length = list(map(int,df.iloc[0]['length_exon'].split(',')))
    if strand == '+':
        summing_values = genomic_start + sum(exon_length[:exon_num-1]) + sum(exon_starts[:exon_num]) + location
    else:
        exon_num = len(exon_length) - exon_num
        summing_values = genomic_start + sum(exon_length[:exon_num]) + sum(exon_starts[:exon_num+1]) + location
    return summing_values
    
    
lst_of_exons = sum_exons(transcript_of_interest)
lst_of_exons = list(set(lst_of_exons))
# print('list of exons',lst_of_exons)
exon_spans = get_exon_spans(lst_of_exons)
finding_exons = find_exon(exon_spans,candidate_positions)
new_exon_locations = new_locations(transcript_of_interest,finding_exons,strand)
# print(new_exon_locations)

final_coordinates = {}
updated_coordinates = []
for key,values in new_exon_locations.items():
    summed_vals = [summing_values(transcript_of_interest,key,location,strand) for location in values]
    final_coordinates[key] = summed_vals
    updated_coordinates.append(summed_vals)

# print(len(candidate_positions))
# print(updated_coordinates)
updated_coordinates = [location for locations in updated_coordinates for location in locations]
# print(len(updated_coordinates))
# print(final_coordinates.values())
include_candidate_df['genomic_position'] = ' '
candidate_and_updated = list(zip(candidate_positions,updated_coordinates))
for cand,updt in candidate_and_updated:
    include_candidate_df['genomic_position'].loc[include_candidate_df['pos_unmod'] == cand] = updt
include_candidate_df.to_csv('/Users/mac/Desktop/week_June22/E1A-s_with_genomic.csv')
    
    
    
    
    
    
    
    