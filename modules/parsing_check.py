import pandas as pd
import os
from os import listdir
from os.path import isfile, join
from functools import reduce
import random 
import matplotlib.pyplot as plt

directory_names = ['/map/','/bam_readcount/','/complete_analysis/']
def find_files(directory):
    dictionary = {'/map/':'genomecov.txt','/bam_readcount/':'bamreadcount.txt','/complete_analysis/':'complete.txt'}
    current_dir = dictionary[directory]
    return current_dir
    
def merging_files(root,list_of_files,endswith,transcript):
    if endswith == 'complete.txt':
        full_path = root + list_of_files[0]
        final_df = pd.read_csv(full_path,sep = '\t',usecols = ['pos_ctrl','depth_ctrl','depth_treat'])
        final_df.columns = ['position','depth_ctrl','depth_treat']
        final_df['depth-complete'] = final_df[['depth_ctrl','depth_treat']].mean(axis=1)
        final_df = final_df.drop(['depth_ctrl','depth_treat'],axis = 1)
    elif endswith == 'bamreadcount.txt':
        full_df = []
        for file in list_of_files:
            df = pd.read_csv(root+file,usecols = [1,3],sep = '\t',header = None,names = ['position','depth_bamreadcount'])
            full_df.append(df)
        final_df = full_df[0].merge(full_df[1], left_on = 'position',right_on = 'position',how = 'inner')
        final_df['depth-bamreadcount'] = final_df[['depth_bamreadcount_x','depth_bamreadcount_y']].mean(axis = 1)
        final_df = final_df.drop(['depth_bamreadcount_x','depth_bamreadcount_y'],axis = 1)
    elif endswith == 'genomecov.txt':
        full_df = []
        #print("LIST OF FILES",list_of_files)
        for file in list_of_files:
            df = pd.read_csv(root+file,usecols = [0,1,2],header = None,names = ['transcript','position','depth'],sep = '\t')
            if len(transcript) > 0:
                df = df[df['transcript']== transcript]
            df = df.drop(['transcript'],axis = 1)
            full_df.append(df)
        final_df = full_df[0].merge(full_df[1], left_on = 'position',right_on = 'position',how = 'inner')
        final_df['depth-raw'] = final_df[['depth_x','depth_y']].mean(axis = 1)
        final_df = final_df.drop(['depth_x','depth_y'],axis = 1)
    return final_df
    
def plotting(df,output_location,transcript):
    max_val = df[['depth-raw','depth-bamreadcount','depth-complete']].max().max()
    upper_limit = ylimit(max_val)
    fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(30, 10),sharex=True)

    ax1 = plt.subplot(4, 1, 1)
    ax1 = plt.plot(df['position'],df['depth-raw'],alpha = 0.25,linewidth = 5,label = 'depth-raw')
    ax1 = plt.plot(df['position'],df['depth-bamreadcount'],alpha = 0.25,linewidth = 5,label = 'bam-readcount',color = 'darkred')
    ax1 = plt.plot(df['position'],df['depth-complete'],alpha = 0.25,linewidth = 5,label = 'depth-complete',color = 'forestgreen')
    plt.title('Depth comparison ({} mean treat-ctrl)'.format(transcript),fontsize = 28)
    plt.ylim(0, upper_limit)
    plt.legend(fontsize = 15)

    ax2= plt.subplot(4, 1, 2)
    ax2 = plt.plot(df['position'],df['depth-raw'],alpha = 0.5,linewidth = 5,label = 'depth-raw')
    plt.ylim(0, upper_limit)
    plt.legend(fontsize= 15)

    ax3 = plt.subplot(4, 1, 3)
    ax3 = plt.plot(df['position'],df['depth-bamreadcount'],alpha = 0.5,linewidth = 5,label = 'bam-readcount',color = 'darkred')
    plt.ylim(0, upper_limit)
    plt.legend(fontsize = 15)

    ax4 = plt.subplot(4, 1, 4)
    ax4 = plt.plot(df['position'],df['depth-complete'],alpha = 0.5,linewidth = 5,label = 'depth-complete',color = 'forestgreen')
    plt.xlabel('Position',fontsize = 28)
    plt.ylim(0, upper_limit)
    plt.legend(fontsize = 15)
    plt.tight_layout()
    plt.savefig(output_location)
    
def ylimit(highest_val):
    converted = str(int(highest_val))
    return int(str(int(converted[0])+1) + (len(converted)-1) *'0')


def exome_mode_parsing(root_directory,name):
    current_result_lst = []
    #print(root_directory)
    for i in directory_names:
        current_directory = root_directory + i
        onlyfiles = [f for f in listdir(current_directory) if isfile(join(current_directory, f))]
        endswith_value = find_files(i)
        #print(i)
        #print('ENDSWITH_VAL',endswith_value)
        #print('ONLYFILES',onlyfiles)
        filter_ = [files for files in onlyfiles if files.endswith(endswith_value)]
        current_result = merging_files(current_directory,filter_,endswith_value,name)
        current_result_lst.append(current_result)
    df_final = reduce(lambda left,right: pd.merge(left,right,on='position'), current_result_lst)
    #print(df_final.head())
    output_dir = root_directory +'/parsing-test/'
    os.makedirs(output_dir,exist_ok = True)
    plotting(df_final,output_dir + name+'.pdf',name)
    
def isoform_mode_parsing(root_directory,list_of_transcripts):
    for transcript_names in list_of_transcripts:
        current_result_lst = []
        for i in directory_names:
            current_directory = root_directory + i
            onlyfiles = [f for f in listdir(current_directory) if isfile(join(current_directory, f))]
            endswith_value = find_files(i)
    #         print(endswith_value)
            if i == '/map/':
                current_files = [x for x in onlyfiles if x.endswith(endswith_value)]
            else:
                filter_ = [files for files in onlyfiles if files.endswith(endswith_value)]
                current_files = [x for x in filter_ if x.startswith(transcript_names)]
            current_result = merging_files(current_directory,current_files,endswith_value,transcript_names)
            current_result_lst.append(current_result)
        df_final = reduce(lambda left,right: pd.merge(left,right,on='position'), current_result_lst)
        output_dir = root_directory +'/parsing-test/'
        os.makedirs(output_dir,exist_ok = True)
        plotting(df_final,output_dir+ transcript_names+'.pdf',transcript_names)

if __name__ == "__main__":
    #exome_mode('/Users/mac/Desktop/Fast-DRUMMER/testing-parsing/xome2.Ad5.MOD-xome2.Ad5.UNMOD','Ad5')
    isoform_mode('/Users/mac/Desktop/Fast-DRUMMER/testing-parsing/isoform2.Ad5.MOD-isoform2.Ad5.UNMOD',['E3.RIDa', 'E3.10K', 'E3.12K'])