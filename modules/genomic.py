import pandas as pd
import numpy as np
import argparse
import math
import re
import warnings
warnings.filterwarnings("ignore")



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
    if strand == '+':
        summing_values = location - sum(exon_length[:exon_num]) + exon_starts[exon_num] + genomic_start
    else:
        [exon_starts,exon_ends] = adjust_negative(genomic_start,exon_starts,exon_length)
        reverse_lengths = list(reversed(exon_length))
        summing_values = exon_ends[exon_num] -  (location-sum(reverse_lengths[:exon_num]) + exon_num)
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

def return_homology(df):
    index_of_interest = []
    for index,rows in df.iterrows():
#         print(rows['eleven_bp_motif'])
        return_expression = re.findall(r'((\w)\2{2,})', rows['eleven_bp_motif'])
        if len(return_expression) > 0 and return_expression[0][-1] != 'N' and rows['candidate_site'] == 'candidate':
            index_of_interest.append(index)
    return index_of_interest

#args = vars(ap.parse_args())
#transcripts = args['transcript_file']
#path = args['input']
#output = args['output']
#include_candidate_df = pd.read_csv(path,sep = '\t')
#transcripts = pd.read_csv(transcripts,sep = '\t',header= None)

#include_candidate_df['candidate_site']=include_candidate_df['candidate_site'].fillna(' ')
#candidate_positions = include_candidate_df['pos_mod'].tolist()


def run_genomics(path,transcripts):
#	print('RUN GENOMICS')
#	print(transcripts)
	include_candidate_df = pd.read_csv(path,sep = '\t')
	transcripts = pd.read_csv(transcripts,sep = '\t',header= None)
	#include_candidate_df['candidate_site']=include_candidate_df['candidate_site'].fillna(' ')
	include_candidate_df['accumulation'] = include_candidate_df['accumulation'].fillna(' ')
	include_candidate_df['depletion'] = include_candidate_df['depletion'].fillna(' ')
	candidate_positions = include_candidate_df['pos_ctrl'].tolist()
	if len(transcripts.columns) > 1 and len(candidate_positions) > 0:
#		print('\n###Genomic locations (ENABLED)###\n')
#		print(include_candidate_df)
		if len(set(include_candidate_df['chr_ctrl'])) == 1:
#			print(include_candidate_df)
			chromo_name = include_candidate_df['chr_ctrl'][0]
		else:  
			print('multiple chromosome names in csv file')
#		print(include_candidate_df)
		if len(transcripts.columns) == 7:
			col_names = ['chromosome','transcript','strand','genome_start','num_exons','length_exon','start_exon']
			transcripts.columns = col_names
			final_chomo = transcripts[transcripts['transcript'] == chromo_name]['chromosome'].values[0]
			include_candidate_df.insert(0, 'Chromosome', final_chomo)
		else:
			col_names = ['transcript','strand','genome_start','num_exons','length_exon','start_exon']
#		print(include_candidate_df)
		transcripts.columns = col_names


#		print('chromo_name',chromo_name)
		transcript_of_interest = transcripts[transcripts['transcript'] == chromo_name].reset_index(drop=True)
#		print(transcript_of_interest)
		strand = transcript_of_interest.loc[0]['strand']
		#print(transcript_of_interest)
		num_exon = transcript_of_interest['num_exons']
		#print('number of exons',num_exon.loc[0])
		lst_of_exons = sum_exons(transcript_of_interest,strand)

		lst_of_exons = sorted(list(set(lst_of_exons)))

		exon_spans = get_exon_spans(lst_of_exons)

		finding_exons = find_exon(exon_spans,candidate_positions)
		#print(finding_exons)

		final_coordinates = {}
		updated_coordinates = []
		for key,values in finding_exons.items():
			summed_vals = [summing_values(transcript_of_interest,key,location,strand,candidate_positions) for location in values]
			final_coordinates[key] = summed_vals
			updated_coordinates.append(summed_vals)

		updated_coordinates = [location for locations in updated_coordinates for location in locations]
		include_candidate_df['genomic_position'] = ' '
		candidate_and_updated = list(zip(candidate_positions,updated_coordinates))
		for i,j in candidate_and_updated:
			if i in [36,65,122,213,317,1071]:
				print(i,j+1)
		for cand,updt in candidate_and_updated:
			if strand == '+':
				add_val = [k for k,val in finding_exons.items() if cand in val]
				include_candidate_df['genomic_position'].loc[include_candidate_df['pos_ctrl'] == int(cand)] = updt + add_val[0]
			else:
				include_candidate_df['genomic_position'].loc[include_candidate_df['pos_ctrl'] == int(cand)] = updt + 1

		#include_candidate_df.to_csv('/Users/mac/Desktop/genomic_test/L5-XYZ-Fiber.original.genome.txt',sep = '\t',index = None)
#		print(include_candidate_df)
		return include_candidate_df
		#include_candidate_df.to_csv(output,sep = '\t',index=False)
if __name__ == "__main__":
    path = '/Users/mac/Desktop/Fast-DRUMMER/test-jan23-isoform-full/Ad5_transcript.M3P2-422-Ad5_transcript.M3KO2-422/complete_analysis/L5-XY-Fiber.complete.txt'
    transcripts = '/Users/mac/Desktop/DRUMMER/TESTDATA/Ad5.complete.transcript.chr.txt'
    run_genomics(path,transcripts)
#df = pd.read_csv('E3.10K.txt',sep = '\t')
#run_genomics(df,'/gpfs/home/ja3539/depledge/Jonathan/DRUMMER/TESTDATA/Ad5.sample.transcripts.txt')
