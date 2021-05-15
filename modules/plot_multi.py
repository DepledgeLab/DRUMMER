import pandas as pd
from collections import Counter,defaultdict
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


consistent_columns = ['transcript_id', 'position', 'eleven_bp_motif', 'nearest_ac_motif',
       'nearest_ac', 'max-G_test', 'max-G_padj', 'support']

def Diff(li1, li2): 
    return (list(set(li1) - set(li2))) 

def return_gtest(attributes):
    print('ATTRIBUTES',attributes)
    print(attributes.split(':')[1])
    return float(attributes.split(':')[1])

def run_plot(df,output_dir):
    output_dir = output_dir + '/candidates_present.pdf'
    df_dict_no_candidates = Counter(df['support'])
    df_no_candidates = pd.DataFrame.from_dict(df_dict_no_candidates,orient = 'index',columns = ['count']).sort_index()
    fig, ax = plt.subplots()
    right_side = ax.spines["right"]
    top_side = ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    x = list(df_no_candidates.index)
    y = list(df_no_candidates['count'])
    ax.bar(x,y)
    plt.xticks(np.arange(1,5,step = 1))
    plt.ylabel('Number of candidate sites')
    plt.xlabel('Number of comparison in which candidate is present')
    plt.savefig(output_dir)
    
def run_gtest(df,columns,output_dir):
    output_dir = output_dir + '/gtest-comparison.pdf'
    y = defaultdict(list)
    for k,values in df.iterrows():
        x = [values[i] for i in columns]
#         print('XXX',x)
        comparison_count = values['support']
        cleanedList = values['max-G_test']
        #cleanedList = [return_gtest(i) for i in x if str(i) != 'nan']
        y[comparison_count].append(cleanedList)
    y = dict(y)
    new_dict = y
#     for k,v in y.items():
#         flatlist = [i for sub in v for i in sub]
#         new_dict[k] = flatlist
    df_gtests_candidates = pd.DataFrame.from_dict(new_dict,orient = 'index').sort_index().T
    fig, ax = plt.subplots(figsize=(50,10))
    ax = sns.catplot(data=df_gtests_candidates,palette=['royalblue'])
    plt.ylabel('G-test score')
    plt.xlabel('Number of comparisons')
    plt.savefig(output_dir,bbox_inches='tight')
    


def run_main(df,output_dir):
    check_columns = Diff(list(df.columns) ,consistent_columns)
    run_plot(df,output_dir)
    run_gtest(df,check_columns,output_dir)
    
if __name__ == "__main__":
    df = pd.read_csv('/Users/mac/Desktop/Fast-DRUMMER/multiple_comp.txt',sep = '\t')
    output_dir = '/Users/mac/Desktop/Fast-DRUMMER/DRUMMER-new/'
    run_main(df,output_dir)