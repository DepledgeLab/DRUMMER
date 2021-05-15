import pandas as pd

def proper_filter(read_metrics:list,filter_large_indels):
    """Takes in a list returns a dictionary with counts of each nucleotide
    """
    nucleotides = ['A','C','G','T','N','adjN']
    #print(type(filter_large_indels))
    keep_track = {k:0 for k in nucleotides} #Initialize each nucleotide count with 0
    for indiv_metrics in read_metrics:
        if indiv_metrics[0] in nucleotides: 
            splitted_metrics = indiv_metrics.split(':')
            keep_track[splitted_metrics[0]] += int(splitted_metrics[1])
        elif indiv_metrics[0] not in nucleotides and indiv_metrics[0] != '=':
            if filter_large_indels == False and str(filter_large_indels).isnumeric() == False:
                adjust_N, actual_N = filter_nofilter(indiv_metrics)
                keep_track['adjN'] += adjust_N
                keep_track['N'] += actual_N
            else:
                #print('here')
                adjust_N, actual_N = filter_large(indiv_metrics,filter_large_indels)
                keep_track['adjN'] += adjust_N
                keep_track['N'] += actual_N
    return keep_track
def filter_large(current_alignment,value):
    adjust_N = 0
    actual_N = 0
    nucleotide_info = current_alignment.split(':')
    #print('value',value)
    #print(nucleotide_info)
    if nucleotide_info[0][0] == '+' and len(nucleotide_info[0]) < int(value)+2:
        #print(nucleotide_info)
        adjust_N += int(nucleotide_info[1])
        actual_N += int(nucleotide_info[1])
    elif nucleotide_info[0][0] == '-' and len(nucleotide_info[0]) < int(value)+2:
        #print(nucleotide_info)
        adjust_N += int(nucleotide_info[1])
        actual_N += int(nucleotide_info[1])
    return adjust_N, actual_N

def filter_nofilter(current_alignment):
    adjust_N = 0
    actual_N = 0
    nucleotide_info = current_alignment.split(':')
    #print(nucleotide_info)
    if nucleotide_info[0][0] == '+':
        #print(nucleotide_info)
        adjust_N += int(nucleotide_info[1])
        actual_N += int(nucleotide_info[1])
    elif nucleotide_info[0][0] == '-':
        #print(nucleotide_info)
        adjust_N += int(nucleotide_info[1])
        actual_N += int(nucleotide_info[1])
    return adjust_N, actual_N  
    
def do_math(df:'DataFrame',index:int):
    """Takes in dataframe and row index, returns that rows reference fraction
    """
    current_df = df[df.index == index]
#     print('current DF',current_df)
    reference_nucleotide = list(current_df['ref'])[0]
    depth = int(list(current_df['depth'])[0])
    numerator = int(list(current_df[reference_nucleotide])[0])
    if int(depth) == 0:
        return 0
    else:
        return int(numerator)/int(depth)

def merged_dataframes(mod,unmod,file,filter_bool):
    '''Takes in the mod and unmod dataframes, returns them parsed, and merged.
    '''
    try:
        columns_names = ['chr','pos','ref','depth','A','C','G','T','N','adjN']
        seperate_dfs = []
        for i in [mod,unmod]:
            list_of_rows = []
            with open(i) as FileObj:
                for lines in FileObj:
                    initial_split = lines.split('\t')
#                     print('here@@@')
                    first_4_dict = dict(zip([1,2,3,4],initial_split[:4])) #First four columns into dict with random keynames
                    nucleotide_count = proper_filter(initial_split[4:],filter_bool)
                    first_4_dict.update(nucleotide_count)
                    list_of_rows.append(pd.DataFrame([first_4_dict]))
                filtered_df = pd.concat(list_of_rows).reset_index(drop=True)
#                 print('here1')
                filtered_df.columns = columns_names
                filtered_df['depth'] = filtered_df[['A','C','T','G','adjN']].astype(int).sum(axis=1)
                filtered_df = filtered_df.drop(['adjN'],axis =1)
                upper_ref = [i.upper() for i in filtered_df['ref']]
                filtered_df['ref']  = upper_ref
#                 print('here2')
                filtered_df['ref_fraction'] = [do_math(filtered_df,indx) for indx in filtered_df.index]
                filtered_df = filtered_df[filtered_df['depth'] != 0].reset_index(drop=True)
                seperate_dfs.append(filtered_df)
        merged_df = seperate_dfs[0].merge(seperate_dfs[1],on = 'pos',suffixes = ('','.1'))
    except ValueError:
        print("Problem merging:",file,'\n')
        print(mod,unmod)
    return merged_df

if __name__ == "__main__":
    merged_dataframes('ExomeFwd-filter-march13/Ad5-M3P1.422.raw.AD5.fwd.gen-Ad5-M3KO1.422.raw.AD5.fwd.gen/bam_readcount/Ad5.MOD.bamreadcount.txt','ExomeFwd-filter-march13/Ad5-M3P1.422.raw.AD5.fwd.gen-Ad5-M3KO1.422.raw.AD5.fwd.gen/bam_readcount/Ad5.UNMOD.bamreadcount.txt','Ad5',False)












