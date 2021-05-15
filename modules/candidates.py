import pandas as pd
import numpy as np


def is_candidate(df,odds_ratio,padj,fraction_diff):
    '''Takes in a dataframe looks at log2fc, odds and padj to determine if the site is a candidate'''
    accum_cutoff = 1/odds_ratio
    idx_snp = ((df['ref_fraction_treat'] < 0.6) | (df['ref_fraction_ctrl'] < 0.6))
    df['is_SNP'] = ['snp' if i == True else "" for i in idx_snp]
    
    idx_depletion = ((df['odds_ratio']>odds_ratio) &(df['padj']< padj) & ( (df['ref_fraction_treat'] - df['ref_fraction_ctrl'] ) > float(fraction_diff)) & (df['p_values_OR_adj']<padj))
    idx_accum = ((df['odds_ratio']<accum_cutoff) &(df['padj']< padj) & ( (df['ref_fraction_treat'] - df['ref_fraction_ctrl'] ) < float(fraction_diff)) & (df['p_values_OR_adj']<padj))
    df['accumulation'] = ['accumulation' if i == True else '' for i in idx_accum ]
    df['depletion'] = ['depletion' if i == True else '' for i in idx_depletion ]
    return df

def Diff(li1, li2): 
    return (list(set(li1) - set(li2))) 
    
def get_masked(include_candidate_df,candidate):
    """Take in dataframe and either accumulation or depletion, returns masked versions of sites
    """

#     if include_candidate_df[candidate].empty == False:
    index_candidates = list(include_candidate_df[include_candidate_df[candidate] == candidate].index)
    #print('include_candidate_df.shape',include_candidate_df.shape)
    #print(index_candidates)
#     print(index_candidates)
#     else:
#         index_candidates = []
    #index_deplet = list(include_candidate_df[include_candidate_df['Depletion'] == 'depletion'].index)
    #index_candidates = index_accum + index_deplet
    if len(index_candidates) > 1:
        #print('IN IF STATEMENT 0')
        #print('Found {} candidate sites'.format(len(index_candidates)))
        full_list = []
        lst = []
        #print('IN IF STATEMENT 0.01')
        for k,value in enumerate(index_candidates):
            #print('IN IF STATEMENT1')
            if k > 0:
                if index_candidates[k] - index_candidates[k-1] < 5:
                    #print('IN IF STATEMENT2')
                    if index_candidates[k] not in lst:
                        lst.append(index_candidates[k])
                else:
                    #print('IN IF STATEMENT3')
                    full_list.append(lst)
                    lst = []
                    if index_candidates[k] not in lst:
                        lst.append(index_candidates[k])
            else:
                #print('IN IF STATEMENT5')
                lst.append(index_candidates[k])
        full_list.append(lst)
        index_of_highest = []
        for i in full_list:
            #print('IN IF STATEMENT6')
            k = []
            for j in i:
                #print('IN IF STATEMENT 6.1')
                k.append(include_candidate_df.iloc[j]['G_test'])
            #print('IN IF STATEMENT 6.3')
            index_of_highest.append(i[np.argmax(k)])
        include_candidate_df[candidate] = ''
        include_candidate_df.loc[index_of_highest,candidate] = candidate
        flattened_list = [indx for lsts in full_list for indx in lsts]
        candidate_masks = Diff(flattened_list,index_of_highest)
        masked_rep = '[{}_masked]'.format(candidate)
        include_candidate_df.loc[candidate_masks,candidate] = masked_rep
        #print(include_candidate_df.head())
    return include_candidate_df

def check_homopolymer(string):
    return_val = []
    for i in ['TTT','AAA','GGG','CCC']:
        if i in string:
            return True
            
def run_find_candidates(include_candidate_df,odds_ratio,padj,fraction_diff):
    #include_candidate_df = include_candidate_df.reset_index(drop = True)
    #include_candidate_df = is_candidate(include_candidate_df,odds_ratio,padj,fraction_diff)
    if type(include_candidate_df) == str:
#        print('INSIDE TYPE')
        include_candidate_df = pd.read_csv(include_candidate_df,sep = '\t')
        include_candidate_df = is_candidate(include_candidate_df,odds_ratio,padj,fraction_diff)
    else:
        include_candidate_df = is_candidate(include_candidate_df,odds_ratio,padj,fraction_diff)
    include_candidate_df = get_masked(include_candidate_df,'accumulation')
#     print('PAST ACCUMULATION',include_candidate_df.head())
    include_candidate_df = get_masked(include_candidate_df,'depletion')
    #print(include_candidate_df.head())
#     print('ACCUMULATION and DEPLETION done',include_candidate_df.head())
#     print('ACCUMULATION and DEPLETION done',include_candidate_df.head())
    #index_candidates = list(include_candidate_df[include_candidate_df['candidate_site'] == 'candidate'].index)
    include_candidate_df['homopolymer'] = [check_homopolymer(i) for i in include_candidate_df['eleven_bp_motif']]
    include_candidate_df.insert(20, 'frac_diff', include_candidate_df['ref_fraction_treat'] - include_candidate_df['ref_fraction_ctrl'])

    return include_candidate_df
    
    
if __name__ == "__main__":
#     path = '/gpfs/home/ja3539/depledge/Jonathan/Final-DRUMMER/test-accumulation-depletion/Ad5-M3KO1.422.raw.AD5.rev.gen-Ad5-M3P1.422.raw.AD5.rev.gen/gTEST/Ad5-RevComp.txt'
    path = '/Users/mac/Desktop/Fast-DRUMMER/Ad5.txt'
    df = pd.read_csv(path,sep = '\t')
    run_find_candidates(df,1.5,0.05,0.01)