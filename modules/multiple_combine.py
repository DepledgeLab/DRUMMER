import pandas as pd
from collections import defaultdict
import math

def merging_dataframes(lst_of_paths):
    suffix_lst = map(str, list(range(1,len(lst_of_paths)+1)))
    suffix_lst = ['.' + i for i in suffix_lst]
    original_path = lst_of_paths[0]
    original = pd.read_csv(original_path,sep = '\t')
    original['name'] = original['transcript_id'] +'_'+ original['transcript_pos'].astype(str)
    new_lst_of_paths = lst_of_paths[1:]
    for i in range(len(new_lst_of_paths)):
        suff_ind = i + 1 
        current_df = pd.read_csv(new_lst_of_paths[i],sep = '\t')
        current_df['name'] = current_df['transcript_id'] +'_'+ current_df['transcript_pos'].astype(str)
        original = original.merge(current_df,on='name',how = 'outer',suffixes=((suffix_lst[suff_ind-1],suffix_lst[suff_ind])))
    return original,suffix_lst


def get_highest_gtest(combine_df,suff_lst,m6A):
    mapping_dict = dict(zip(list(range(0,len(suff_lst))),suff_lst))
    final_dict = {}
    index_replicate = {}
    only_gtest = [i for i in combine_df.columns if 'G_test' in i]
    #print(combine_df.columns)
    for index,values in combine_df.iterrows():
        update_dic = {}
        values = values.fillna(0)
        #my_list = [values['G_test_x'],values['G_test_y'],values['G_test_a'],values['G_test_b']]
        my_list = [values[col] for col in only_gtest]
        none_zero = [i for i, gtest in enumerate(my_list) if gtest != 0]
        max_value = max(my_list)
        max_index = my_list.index(max_value)
        which_letter = mapping_dict[max_index]
        update_dic['eleven_bp_motif'] = values['eleven_bp_motif'+which_letter]
        update_dic['depletion'] = values['depletion'+which_letter]
        update_dic['accumulation'] = values['accumulation'+which_letter]
        if 'genomic_position.1' in combine_df.columns:
#             print('get_highest_gtest--multiplecobine')
            update_dic['genomic_position'] = values['genomic_position'+which_letter]
        if m6A == True:
            update_dic['nearest_ac_motif'] = values['nearest_ac_motif'+which_letter]
            update_dic['nearest_ac'] = values['nearest_ac'+which_letter]
        update_dic['G_test'] = values['G_test'+which_letter]
        update_dic['G_padj'] = values['G_padj'+which_letter]
        update_dic['odds_ratio'] = values['odds_ratio'+which_letter]
        update_dic['OR_padj'] = values['OR_padj'+which_letter]
        update_dic['frac_diff'] = values['frac_diff'+which_letter]
        update_dic['is_SNP'] = values['is_SNP'+which_letter]
        update_dic['count'] = len(none_zero)
        final_dict[index] = update_dic
        index_replicate[index] = none_zero
    combine_df_f = pd.DataFrame.from_dict(final_dict, orient='index')
    combine_df_f['check'] = combine_df_f.index
    combine_df_f[['transcript_id','position']]=combine_df_f.check.str.split('_',expand=True)
    if m6A == True:
        combine_df_f = combine_df_f[['transcript_id','position','eleven_bp_motif','nearest_ac_motif','nearest_ac','G_test','G_padj','odds_ratio','OR_padj','check','count','frac_diff','accumulation','depletion','is_SNP']]
    else:
        combine_df_f = combine_df_f[['transcript_id','position','eleven_bp_motif','G_test','G_padj','odds_ratio','OR_padj','check','count','accumulation','depletion','frac_diff','is_SNP']]
    return combine_df_f,index_replicate,final_dict
    
def grouping_positions(index_candidates):
    full_list = []
    lst = []
    for k,value in enumerate(index_candidates):
        if k > 0:
            if index_candidates[k] - index_candidates[k-1] < 5:
                if index_candidates[k] not in lst:
                    lst.append(index_candidates[k])
            else:
                full_list.append(lst)
                lst = []
                if index_candidates[k] not in lst:
                    lst.append(index_candidates[k])
        else:
            lst.append(index_candidates[k])
    full_list.append(lst)
    return [i for i in full_list if len(i) > 1]
    
def return_final_value(test,transcript,lst_candidates):
    lst_candidates = list(map(int,lst_candidates))
    current_ = test[test['transcript_id'] == transcript]
    current_ = current_[current_['position'].astype(int).isin(lst_candidates)]
    current_.index = current_['position'].reset_index(drop = True)
    max_gtest_idx = current_[['G_test']].idxmax()[0]
    lower_gtest_idxs = [i for i in lst_candidates if i!= int(max_gtest_idx)]
    total_count = sum(current_['count'])
    return max_gtest_idx,lower_gtest_idxs, total_count

def Diff(li1, li2):
    return (list(list(set(li1)-set(li2)) + list(set(li2)-set(li1))))


def collate_information(compare_name,df_path):
    df = pd.read_csv(df_path,sep = '\t')
    out_dict = {}
    df[['transcript_pos','G_test','G_padj','odds_ratio','OR_padj','frac_diff']] = df[['transcript_pos','G_test','G_padj','odds_ratio','OR_padj','frac_diff']].astype(str)
    for index,values in df.iterrows():
        new_lst = [compare_name,values['transcript_pos'],values['G_test'],values['G_padj'],values['odds_ratio'],values['OR_padj'],values['frac_diff']]
        name = values['transcript_id'] +'_'+ values['transcript_pos']
        combined_info = ':'.join(new_lst)
        out_dict[name] = combined_info
    return out_dict
def return_odds_adj(df,new_columns):
    max_odds = []
    max_odds_padj = []
    for key,val in df[new_columns].iterrows():
        val = val.dropna()
        #print(val)
        get_OR = [i.split(':')[-3] for i in val]
        get_OR_padj = [i.split(':')[-2] for i in val]
        convert_float_OR = list(map(float,get_OR))
        convert_float_OR_padj = list(map(float,get_OR_padj))
        max_odds.append(max(convert_float_OR))
        max_odds_padj.append((max(convert_float_OR_padj)))
    return max_odds,max_odds_padj

def get_comparisons(lst_of_dicts,keep_drop_sites,columns,lst_reps):
    specific_site, all_sites = keep_drop_sites[0],keep_drop_sites[1]
    final_dict_all = defaultdict(list)
    test_dictionary = pd.DataFrame(columns=lst_reps,index= [specific_site])
    for k in lst_of_dicts:
        for key,values in k.items():
            new_dict = {}
            if key in all_sites:
                splitted_values = values.split(':')
#                 print(splitted_values[0])
                test_dictionary.loc[specific_site,splitted_values[0]] = ':'.join(return_round(splitted_values[1:]))
#     print(test_dictionary)
    return test_dictionary

def return_round(lst):
    new_val_lst = []
    for i in lst:
        if 'e' in i:
            new_val = "%.3g" % float(i)
            new_val_lst.append(new_val)
        else:
            new_val_lst.append(str(round(float(i),4)))
    return new_val_lst

def get_comparisons_single(lst_of_dicts,single_site,columns,lst_reps):
    final_dict_all = defaultdict(list)
    test_dictionary = pd.DataFrame(columns=lst_reps,index= [single_site])
    for k in lst_of_dicts:
        for key,values in k.items():
            new_dict = {}
            if key == single_site:
#                 print(values)
                splitted_values = values.split(':')
                test_dictionary.loc[single_site,splitted_values[0]] = ':'.join(return_round(splitted_values[1:]))
    return test_dictionary


def convert_column_names(lst_cols):
    return ['{}-pos:Gtest:padj:OR:ORpadj:frac_diff'.format(i) for i in lst_cols]

def return_new_support(df,columns):
    final_count = []
    for key,vals in df.iterrows():
        all_attributes = [vals[i] for i in columns]
        rid_nans = [no_nan for no_nan in all_attributes if str(no_nan) != 'nan']
        final_count.append(len(rid_nans))
    return final_count

def main(outputdir,lst_rep_dir,m6A):
    all_paths = [outputdir + '/' + i + '/summary.txt' for i in lst_rep_dir]
    df,suffix_list = merging_dataframes(all_paths)
    df.index = df.name
    combined_df,index_replicates,final_dict = get_highest_gtest(df,suffix_list,m6A)

    lst_transcripts = set(combined_df['transcript_id'])
    group = []
    dropped_lst = []
    keep_lst = {}
    together = {}
    for i in lst_transcripts:
        current_df = combined_df[combined_df['transcript_id'] == i]
        lst_positions = sorted(current_df['position'].astype(int))
        lst_grouped_positions = grouping_positions(lst_positions)
        for j in lst_grouped_positions:
            group.append(j)
            max_gtest_idx,dropped_indices, total_counts = return_final_value(combined_df,i,j)
            dropping = [i + '_'+str(drop_i) for drop_i in dropped_indices]
            keeping = i + '_'+str(max_gtest_idx)
            dropped_lst.append(dropping)
            keep_lst[keeping] = total_counts
            dropping.append(keeping)
            together[keeping] = dropping
    dropped_lst = [item for sublist in dropped_lst for item in sublist]

    real_delete = Diff(keep_lst,dropped_lst)

    comp_path = dict(zip(all_paths,all_paths))

    updated_dict = [collate_information(k,v) for k,v in comp_path.items()]

    fianl_df = pd.DataFrame(columns = all_paths).T
    for collated_sites in together.items():
        current_df = get_comparisons(updated_dict,collated_sites,all_paths,lst_rep_dir)
        fianl_df = fianl_df.join(current_df.T, how='outer')
    fianl_df = fianl_df.T

    singlefianl_df = pd.DataFrame(columns = all_paths).T
    all_together = [item for sublist in list(group) for item in sublist]
    single = Diff(list(combined_df.index), all_together)

    singlefianl_df = pd.DataFrame(columns = all_paths).T
    for single_site in single:
        current_df = get_comparisons_single(updated_dict,single_site,all_paths,lst_rep_dir)

        singlefianl_df = singlefianl_df.join(current_df.T,how = 'outer')
    singlefianl_df = singlefianl_df.T

    singlefianl_df = singlefianl_df.drop(dropped_lst)
    new_df = pd.concat([fianl_df,singlefianl_df])


    new_columns = convert_column_names(lst_rep_dir)

    new_cols = dict(zip(all_paths,new_columns))
    new_df = new_df.rename(columns=new_cols)
    final_out = combined_df.join(new_df, how='left')
    final_out['count'].update(pd.Series(keep_lst))
    final_out = final_out.drop(real_delete)
    final_out = final_out[~final_out.index.duplicated(keep='first')]
    final_out=final_out.rename(columns = {'G_test':'max-G_test','G_padj':'max-G_padj','count':'support'})

    other_test = pd.DataFrame.from_dict(final_dict, orient='index')
    if 'genomic_position' in other_test.columns:
#         print('other_test.columns')
        final_out = final_out.join(other_test['genomic_position'], how='left')
        cols = ['transcript_id','position','genomic_position', 'eleven_bp_motif','nearest_ac_motif','nearest_ac','max-G_test','max-G_padj','support','accumulation','depletion','frac_diff','is_SNP']
    else:
        cols = ['transcript_id','position', 'eleven_bp_motif','nearest_ac_motif','nearest_ac','max-G_test','max-G_padj','support','accumulation','depletion','frac_diff','is_SNP']
    final_out = final_out.sort_values(['check'], ascending=[True])
    final_cols = cols + new_columns
    if m6A == False:
        final_cols.remove('nearest_ac_motif')
        final_cols.remove('nearest_ac')
        cols.remove('nearest_ac_motif')
        cols.remove('nearest_ac')
    final_out = final_out[final_cols]

    if 'genomic_position' in other_test.columns:
        final_out['genomic_position'] =  final_out['genomic_position'].astype(int)
    columns = Diff(cols,final_out.columns)
    final_out['support'] = return_new_support(final_out,columns)
    max_odds,max_odds_padj = return_odds_adj(final_out,new_columns)
    #print(final_out.head())
    #print(max_odds_padj)
    final_out.insert(9,'max_odds',max_odds)
    final_out.insert(10,'max_odds_padj',max_odds_padj)
    return final_out


    
if __name__ == "__main__":
    import argparse 
    ap = argparse.ArgumentParser(description = 'Takes in the output from the pipeline and determines candidate sites \
    using log2fc, odds_ratio and padj')

    requiredGrp = ap.add_argument_group('required arguments')

    requiredGrp = ap.add_argument_group('required arguments')
    requiredGrp.add_argument("-d",'--directory', required=True, help="directory of DRUMMER run")
    requiredGrp.add_argument("-c",'--comparison_names', required=True, help="names of replicates",action="store",nargs = "*")
    requiredGrp.add_argument("-o",'--output', required=True, help="output name")
    args = vars(ap.parse_args())
    directory = args['directory']
    comparison = args['comparison_names']
    output = args['output']
    print('directory',directory)
    print('comparison',comparison)
    print('output',output)
    df = main(directory,comparison,True)
    df.to_csv(output+'.txt',sep = '\t')















