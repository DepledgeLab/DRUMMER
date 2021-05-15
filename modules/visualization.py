import pandas as pd
import os 
import argparse
from os import listdir
from os.path import isdir, join, isfile
import glob
import itertools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from collections import defaultdict
import math
import warnings
warnings.filterwarnings("ignore")

#####  JONATHAN --->> CLEAN CODE UP!!

def return_shape_size(df,color):
    legend_dict = defaultdict(int)
    col_g = []
    size = []
    count = 0
    
    for index,row in df.iterrows():
        if row['accumulation'] == 'accumulation':
            count +=0
            col_g.append('blue')
            size.append(100)
            legend_dict['accumulation'] += 1
        elif row['accumulation'] == '[accumulation_masked]':
            count += 0
            col_g.append('lightskyblue')
            size.append(10)
            legend_dict['[accumulation_masked]'] += 1
        elif row['depletion'] == 'depletion' and row['pos_ctrl'] not in top10:
            #print(row)
            col_g.append('red')
            size.append(100)
            legend_dict['depletion'] +=1
        elif row['depletion'] == 'depletion' and row['pos_ctrl'] in top10:
            count += 0
            #print('TOP!)')
            col_g.append('darkred')
            size.append(100)
            legend_dict['[depletion_masked]'] += 1
        else:
            col_g.append('black')
            size.append(1)
            legend_dict['none'] += 1
    return [col_g,size,legend_dict,count]

def return_close_ac(lst):
    return ['darkorange' if abs(i) < 6 else 'black' for i in lst]    


                

def run_visualization(df,m6A,output_location,mode,sample):
    max_gtest = max(df['G_test'])
    accumulation_data = df[(df['accumulation'] == 'accumulation') ]
    depletion_data = df[(df['depletion'] == 'depletion') ]
    accumulation_data['log2odds'] = -np.log2(accumulation_data['odds_ratio'])
    depletion_data['log2odds'] = -np.log2(depletion_data['odds_ratio'])

    dropped_inf_depletion = depletion_data['log2odds'].replace([-np.inf], np.nan).dropna(how = 'all')
    dropped_inf_accumulation = accumulation_data['log2odds'].replace([np.inf], np.nan).dropna(how = 'all')
    acc_shape = accumulation_data.shape[0]
    dep_shape = depletion_data.shape[0]
    if acc_shape>1 and dep_shape >1:
        row_vals = 3
    else:
        row_vals =2 
    combined_vals = list(dropped_inf_depletion) + list(map(abs, list(dropped_inf_accumulation)))
    if len(combined_vals) > 0 :
        #fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(30, 10),sharex=True,sharey=True,subplot_kw={'xlim': (0, 4000)})
        fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(30, 10),sharex=True,sharey=True)
        if mode == "exome":
            size = 5
            tickmarks= 10000
        else:
            size = 16
            tickmarks = 500
        cnt = 1
        if accumulation_data.shape[0] > 1:
            ax1 = plt.subplot(row_vals, 1, cnt)

            ax1.set_ylabel('gTEST Score \n(Accumulation)',fontsize = 18,color ='blue')
            col_g,size,legend_dict,count = return_shape_size(accumulation_data,'red')
            ax1.spines['right'].set_visible(False)
            ax1.spines['top'].set_visible(False)
            ax1.spines['bottom'].set_visible(False)
            ax1.spines['left'].set_color('blue')
            for t in ax1.yaxis.get_ticklines(): t.set_color('blue')
            plt.ylim(0, max_gtest+100)
            scatter1 = ax1.scatter(list(accumulation_data['pos_ctrl']), list(accumulation_data['G_test']+50), c=col_g,s = size,edgecolors = 'black',linewidth = 1)
            plt.xticks(np.arange(0, 4500, step=500),fontsize=14,weight = 'bold')
            #plt.xticks([])
            plt.yticks(np.arange(0, 2500, step=500),fontsize=14,weight = 'bold',color = 'blue')
            if m6A == True and accumulation_data.shape[0] < 500:
                for i,j in accumulation_data.iterrows():
                    if j['five_bp_motif'][2:4] == 'AC':
                        pass
                        
            cnt += 1
##################
        limits = math.ceil(max(combined_vals))
        accumulation_data['log2odds'] = -np.log2(accumulation_data['odds_ratio'])
        depletion_data['log2odds'] = -np.log2(depletion_data['odds_ratio'])
        depletion_data['log2odds'] = depletion_data['log2odds'].replace([-np.inf], limits)
        accumulation_data['log2odds'] = accumulation_data['log2odds'].replace([np.inf], -limits)
        #depletion_data = depletion_data[depletion_data['pos_ctrl'] != 160] 
        new_df = pd.concat([depletion_data,accumulation_data])
        ax2 = plt.subplot(row_vals,1,cnt)
        col_deplet,size_deplet,legend_dict,count = return_shape_size(depletion_data,'red')
        col_accum,size_accum,legend_dict,count = return_shape_size(accumulation_data,'red')
    
        color_deplet = return_close_ac(depletion_data['nearest_ac'])
        col_accum = return_close_ac(accumulation_data['nearest_ac'])
        log_depletion = depletion_data.copy()
        log_depletion = log_depletion[log_depletion['pos_ctrl'] != 160] 
        plt.ylim(-limits-2, limits+2)
        #plt.xticks([])
        col_deplet,size_deplet,legend_dict,count = return_shape_size(log_depletion,'red')
        #print('LOG',ax2.xaxis.get_majorticklabels()[0])
        if m6A == True:
            mapping_color = {True:'black',False:'blue'}
            color_deplete_copy = color_deplet = return_close_ac(log_depletion['nearest_ac'])
            scatter2 = ax2.scatter(data=log_depletion, x="pos_ctrl", y="log2odds",color = color_deplet,s = size_deplet,edgecolors = 'black',linewidth = 1)
            scatter3 = ax2.scatter(data=accumulation_data, x="pos_ctrl", y="log2odds",color = col_accum,s = size_accum,edgecolors = 'black',linewidth = 1)
            ax2.legend((scatter2,scatter2), ('Within DRACH (<5nts) ','Outside DRACH (>5nts)'),
                                                   loc='upper left', bbox_to_anchor=(0, 1),fontsize = 7.5)
            leg = ax2.get_legend()
            leg.legendHandles[0]._sizes = [50]
            leg.legendHandles[1]._sizes = [50]
            leg.legendHandles[0].set_color('darkorange')
            leg.legendHandles[1].set_color('black')
            #print(ax2.get_position())
        else: 
            scatter2 = ax2.scatter(data=depletion_data, x="pos_ctrl", y="log2odds",color = 'black',s = size_deplet,edgecolors = 'black',linewidth = 1)
            scatter3 = ax2.scatter(data=accumulation_data, x="pos_ctrl", y="log2odds",color = 'black',s = size_accum,edgecolors = 'black',linewidth = 1)
        plt.axhline(y=0, xmin=0, xmax=1,lw = .5,color = 'black')
        plt.xticks(np.arange(0, 4500, step=500),fontsize=14,weight = 'bold')
        plt.ylabel('Log2 Odds Ratio',fontsize = 18)
        ax2.spines['right'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax2.spines['bottom'].set_visible(False)
        cnt += 1
##################
        
        if depletion_data.shape[0] > 1:
            ax3 = plt.subplot(row_vals, 1, cnt)

            ax3.set_ylabel('gTEST Score \n(Depletion)',fontsize = 18,color = 'darkred')
            col_g,size,legend_dict,count = return_shape_size(depletion_data,'darkred')
            ax3.spines['right'].set_visible(False)
            ax3.spines['bottom'].set_visible(False)
            #ax3.spines['top'].set_visible(False)
            ax3.spines['left'].set_color('darkred')
            for t in ax3.yaxis.get_ticklines(): t.set_color('darkred')
            plt.ylim(0, max_gtest+100)
            scatter3 = ax3.scatter(list(depletion_data['pos_ctrl']), list(depletion_data['G_test']+50), c=col_g,s = size,edgecolors = 'black',linewidth = 1)
            ax3 = plt.gca()
            ax3.set_ylim(ax3.get_ylim()[::-1])

            plt.xticks(np.arange(0, 4500, step=500),fontsize=14,weight = 'bold')
            plt.yticks(np.arange(0, 2100,step = 500),fontsize=14,weight = 'bold',color = 'darkred')
            print(ax3.xaxis.get_majorticklabels()[0])
            if m6A == True and depletion_data.shape[0] < 500:
                for i,j in depletion_data.iterrows():
                    if j['five_bp_motif'][2:4] == 'AC':
                        pass
                        top_limit = plt.ylim()[0]
        fig.tight_layout()
        plt.savefig(output_location)

    
if __name__ == "__main__":
#     path = '/gpfs/home/ja3539/depledge/Jonathan/Final-DRUMMER/test-accumulation-depletion/Ad5-M3KO1.422.raw.AD5.rev.gen-Ad5-M3P1.422.raw.AD5.rev.gen/gTEST/Ad5-RevComp.txt'
    df = pd.read_csv('/Users/mac/Desktop/DRUMMER_Figures/Post-Updated-Figures/Figure2/L2-Penton.complete.txt',sep = '\t')
    top10 = pd.read_csv('/Users/mac/Desktop/DRUMMER_Figures/Post-Updated-Figures/Figure2/L2-Penton.top10byGtest.txt',sep = '\t')
    top10 = list(top10['pos_ctrl'])
    run_visualization(df,True,'/Users/mac/Desktop/DRUMMER_Figures/Post-Updated-Figures/Figure2/L2-Penton.pdf','isoform','L2-Penton')