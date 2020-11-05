import pandas as pd
import numpy as np
import argparse
import math
import scipy.stats as stats
from create_output_file import create_output
import warnings
warnings.filterwarnings("ignore")
#columns_names = ['chr','pos','ref','depth','A','C','G','T','N']

ap = argparse.ArgumentParser(description = 'Takes in the output of bam-readcount \
and returns a text file containing the count of each nucleotide at the position and the \
fraction of the reference nucleotide among all reads.')
requiredGrp = ap.add_argument_group('required arguments')
requiredGrp.add_argument("-i","--input", required=True, help="input file location")
requiredGrp.add_argument("-o","--output", required=True, help="output file directory")
#requiredGrp.add_argument("-i2","--input2", required=True, help="input2 file location")

args = vars(ap.parse_args())
input = args['input']
output = args['output']



#control = [[157,0,52170,8,145],[28,8,54266,20,1]]
#test = [[106, 1, 44061, 4, 84],[34,10,45661,27,1]]
def return_ratio(max_list:list,sum_list:list):
    '''Takes in a list containing max values and the sum of the counts, returns ratio.
    '''
    ratios = []
    for max_value,sum_value in zip(max_list,sum_list):
        if max_value == sum_value: #If the max value is equal to the sum value, count the ratio as 1
            ratios.append(1)
        else:
            ratios.append(max_value/(sum_value-max_value))
    return ratios

def return_odds(max_c,max_t,sum_c,sum_t):
    '''Takes in max control/test, and sum control/test returns odds.
    '''
    odds_list = []
    pvalue_list = []
    for m_c,m_t,s_c,s_t in zip(max_c,max_t,sum_c,sum_t):
        table = np.array([[m_c,s_c-m_c],[m_t,s_t-m_t]])
        oddsratio, pvalue = stats.fisher_exact(table)
        odds_list.append(oddsratio)
        pvalue_list.append(pvalue)
    return odds_list,pvalue_list

def max_sum_function(ctrl:list,test:list):
    '''Takes in a control and test list and returns odds ratio
    '''
    #Max value for each row in dataframe
    max_ctrl = [max(row) for row in ctrl]
    max_test = [max(row) for row in test]
    #print('mc',max_ctrl)
    #print('mt',max_test)
    #Sum of the values of each nucleotide count
    sum_ctrl = [sum(row) for row in ctrl]
    sum_test = [sum(row) for row in test]
    #print('sc',sum_ctrl)
    #print('st',sum_test)
    #ratioA = max_ctrl / (sum_ctrl - max_ctrl)
    ratio_ctrl = return_ratio(max_ctrl,sum_ctrl)
    ratio_test = return_ratio(max_test,sum_test)
    #Get fold change value
    fold_change = [ctrl_val/test_val for ctrl_val,test_val in zip(ratio_ctrl,ratio_test)]
    #Get log2fc value
    log2_fc = [math.log2(fc) for fc in fold_change]
    #$odds = 1 /(($maxA*($sumB-$maxB))/(($maxB*($sumA-$maxA)))); 
    odds_vals,pvalues = return_odds(max_ctrl,max_test,sum_ctrl,sum_test)
    return [ratio_ctrl,ratio_test,fold_change,log2_fc,odds_vals,pvalues]
    
def return_ratio(max_list:list,sum_list:list):
    '''Takes in a list containing max values and the sum of the counts, returns ratio.
    '''
    ratios = []
    for max_value,sum_value in zip(max_list,sum_list):
        if max_value == sum_value: #If the max value is equal to the sum value, count the ratio as 1
            ratios.append(1)
        else:
            ratios.append(max_value/(sum_value-max_value))
    return ratios
            

    
df = pd.read_csv(input,sep='\t')
df.columns = [i.replace('.1','_unmod') if '.1' in i else i+'_mod' for i in df.columns ]

control = [[row['A_unmod'],row['C_unmod'],row['G_unmod'],row['T_unmod'],row['N_unmod']] for index,row in df.iterrows()]
test = [[row['A_mod'],row['C_mod'],row['G_mod'],row['T_mod'],row['N_mod']] for index,row in df.iterrows()]


ratio_unmod,ratio_mod,fold_change, log2_fc, odds_vals,pvalues = max_sum_function(control,test)

df['ratio_unmod'] = ratio_unmod
df['ratio_mod'] = ratio_mod
df['fold_change'] = fold_change
df['log2_fc'] = log2_fc
df['odds_ratio'] = odds_vals
df['p_values_OR'] = pvalues
print('length df',len(df))
df['p_values_OR_adj'] = df['p_values_OR'] * len(df)
output = create_output(output,input,'odds_ratio')
# print(output)
df.to_csv(output,sep = '\t', index = False)










