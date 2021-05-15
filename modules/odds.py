
import pandas as pd
import math
import numpy as np
import scipy.stats as stats

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
    #max_ctrl,sum_ctrl = [(max(row),sum(row)) for row in ctrl]
    #max_test,sum_test = [(max(row),sum(row)) for row in test]
    
    max_ctrl = [max(row) for row in ctrl]
    max_test = [max(row) for row in test]
    sum_ctrl = [sum(row) for row in ctrl]
    sum_test = [sum(row) for row in test]
    
    ratio_ctrl = return_ratio(max_ctrl,sum_ctrl)
    ratio_test = return_ratio(max_test,sum_test)
    #Get fold change value
    fold_change = [ctrl_val/test_val for ctrl_val,test_val in zip(ratio_ctrl,ratio_test)]
    #Get log2fc value
    log2_fc = [-math.log2(fc) for fc in fold_change]
    odds_vals,pvalues = return_odds(max_ctrl,max_test,sum_ctrl,sum_test)
    return [ratio_ctrl,ratio_test,fold_change,log2_fc,odds_vals,pvalues]
    
def return_ratio(max_list:list,sum_list:list):
    '''Takes in a list containing max values and the sum of the counts, returns ratio.
    '''
    return [1 if max_value == sum_value else max_value/(sum_value-max_value) for max_value,sum_value in zip(max_list,sum_list)]
    


def run_odds(df):
    df.columns = [i.replace('.1','_treat') if '.1' in i else i+'_ctrl' for i in df.columns ]
    #df.columns = [i.replace('.1','_ctrl') if '.1' in i else i+'_treat' for i in df.columns ]
    #print("IN RUN ODDS")
    unmodified = [[row['A_treat'],row['C_treat'],row['G_treat'],row['T_treat'],row['N_treat']] for index,row in df.iterrows()]
    modified = [[row['A_ctrl'],row['C_ctrl'],row['G_ctrl'],row['T_ctrl'],row['N_ctrl']] for index,row in df.iterrows()]

    ratio_unmod,ratio_mod,fold_change, log2_fc, odds_vals,pvalues = max_sum_function(unmodified,modified)
    df['ratio_treat'] = ratio_unmod
    df['ratio_ctrl'] = ratio_mod
    #df['fold_change'] = fold_change
    df['log2_(OR)'] = log2_fc
    df['odds_ratio'] = odds_vals
    df['p_values_OR'] = pvalues

    return df
