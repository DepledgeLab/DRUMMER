import pandas as pd
import numpy as np
import argparse
import math
import re
from create_output_file import create_output

ap = argparse.ArgumentParser(description = 'Takes in the output of bam-readcount \
and returns a text file containing the count of each nucleotide at the position and the \
fraction of the reference nucleotide among all reads.')
requiredGrp = ap.add_argument_group('required arguments')
requiredGrp.add_argument("-i","--input", required=True, help="input file location")
requiredGrp.add_argument("-o","--output", required=True, help="output file directory")

args = vars(ap.parse_args())
input = args['input']
output = args['output']

def closest_ac(df):
    seq = ''.join(df['ref_unmod'])
    ac_location = [m.start() for m in re.finditer('AC', seq)]
    nearest_ac = []
    five_base_motif = []
    #print(ac_location)
    for i in range(len(seq)):
        absolute_difference_function = lambda list_value : abs(list_value - i)
        closest_value = min(ac_location, key=absolute_difference_function)
        zero_index_correction = i + 1
        five_base_motif.append(seq[closest_value-2:closest_value+3])
        #print('closest_value',closest_value,'index',zero_index_correction)
        nearest_ac.append(zero_index_correction - closest_value + 1)
    return nearest_ac,five_base_motif
    
def get_kmers(df,num_kmers1,num_kmers2):
    kmers1_pad,kmers2_pad = list(map(lambda x: math.floor(x/2),[num_kmers1,num_kmers2]))
    longest_value = max(kmers1_pad,kmers2_pad)
    #seq = longest_value * 'N' + seq + longest_value * 'N'
    motif_length = [kmers1_pad,kmers2_pad]
    all_motifs = []
    for num_motif in motif_length:
        current_kmer_motif = []
        seq = ''.join(df['ref_unmod'])
        seq = num_motif * 'N' + seq + num_motif * 'N'
        for index,string in enumerate(seq[num_motif:len(seq)-num_motif]):
            correct_index = index + num_motif
            motif = seq[correct_index-num_motif:correct_index+num_motif+1]
            current_kmer_motif.append(motif)
        all_motifs.append(current_kmer_motif)
    return all_motifs


df = pd.read_csv(input,sep = '\t')
# print(df.columns)
df = df.dropna()
near_ac,five_bp_motif = closest_ac(df)
five,eleven = get_kmers(df,5,11)

df['nearest_ac'] = near_ac
df['nearest_ac_motif'] = five_bp_motif
df['five_bp_motif'] = five
df['eleven_bp_motif'] = eleven

output = create_output(output,input,'motif_information')
df.to_csv(output,sep = '\t', index = False)
