def closest_ac(df):
    seq = ''.join(df['ref_treat'])
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
        nearest_ac.append(zero_index_correction - closest_value - 1)
    return nearest_ac,five_base_motif
    
def get_kmers(df,num_kmers1,num_kmers2):
    kmers1_pad,kmers2_pad = list(map(lambda x: math.floor(x/2),[num_kmers1,num_kmers2]))
    longest_value = max(kmers1_pad,kmers2_pad)
    #seq = longest_value * 'N' + seq + longest_value * 'N'
    motif_length = [kmers1_pad,kmers2_pad]
    all_motifs = []
    for num_motif in motif_length:
        current_kmer_motif = []
        seq = ''.join(df['ref_treat'])
        seq = num_motif * 'N' + seq + num_motif * 'N'
        for index,string in enumerate(seq[num_motif:len(seq)-num_motif]):
            correct_index = index + num_motif
            motif = seq[correct_index-num_motif:correct_index+num_motif+1]
            current_kmer_motif.append(motif)
        all_motifs.append(current_kmer_motif)
    return all_motifs
    
import pandas as pd
import re 
import math
#m6A_status = "True"
#path = '/Users/mac/Desktop/DRUMMER/virus_7_vis_m6a_True/odds_ratio/E3.10K.merged.odds_ratio.txt'
#df = pd.read_csv(path,sep = '\t')
def run_motif(df,m6A_status):
    df = df.dropna()
    #print(df.head())
    near_ac,five_bp_motif = closest_ac(df)
    five,eleven = get_kmers(df,5,11)
    if m6A_status == True:
#         print('M6A STATUS',m6A_status)
        df['nearest_ac'] = near_ac
        df['nearest_ac_motif'] = five_bp_motif
    df['five_bp_motif'] = five
    df['eleven_bp_motif'] = eleven
    #df = df[(df['depth_treat'] > 100) & (df['depth_ctrl'] > 100)]
    #print(df.head())
    return df

if __name__ == "__main__":
    df = pd.read_csv('/Users/mac/Desktop/DRUMMER_Figures/march14th/Ad5.complete.filter.txt',sep = '\t',dtype={"depletion": "string", "accumulation": "string"})
    new_df = run_motif(df,True)
    #print(new_df.head())
    new_df.to_csv('/Users/mac/Desktop/DRUMMER_Figures/march14th/Ad5.complete.filter.m6a.txt',sep = '\t',index = None)