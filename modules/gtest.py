import pandas as pd
from scipy.stats import chi2_contingency
import numpy as np

def run_gtest(df):
    control = [[row['A_treat'],row['C_treat'],row['G_treat'],row['T_treat'],row['N_treat']] for index,row in df.iterrows()]
    test = [[row['A_ctrl'],row['C_ctrl'],row['G_ctrl'],row['T_ctrl'],row['N_ctrl']] for index,row in df.iterrows()]
    #print('REGULAR REGULAR df',df_depth.shape[0])
    #Replace 0 with 0.0001 for G-test calculations
    control = [ list(map(lambda x: x if x != 0 else 0.0001, i)) for i in control ]
    test = [list(map(lambda x: x if x != 0 else 0.0001, i)) for i in test]
    
    #Perform G-test on each row and grab the gtest and pval
    lst = [chi2_contingency(np.array([c,t]),lambda_ = 'log-likelihood')[0:2] for c,t in zip(control,test)]
    #df_depth = df[(df['depth_treat'] > 100) & (df['depth_ctrl'] > 100)]
    
    #Save G and pval to variable and assign to column in dataframe
    df['G_test'] = [round(g) for g,p in lst]
    df['p_val'] = [p for g,p in lst] 

    #Padj is obtained by multiplying pval by length of dataframe
    #df['padj'] = df['p_val'] * df_depth.shape[0]
    #df['p_values_OR_adj'] = df['p_values_OR'] * df_depth.shape[0]
    df['padj'] = df['p_val'] * df.shape[0]
    df['p_values_OR_adj'] = df['p_values_OR'] * df.shape[0]
    #print('df_depth (post filter) shape',df_depth.shape[0])
    #print('REGULAR df',df_depth.shape[0])
    #df['p_values_OR_adj'] = df['p_values_OR'] * df_depth.shape[0]
    #print('G-test shape',df.shape[0])
    df['padj'] = [1 if i > 1 else i for i in df['padj']]
    df['p_values_OR_adj'] = [1 if i > 1 else i for i in df['p_values_OR_adj']]
    return df
if __name__ == "__main__":
    df = pd.read_csv('/Users/mac/Desktop/Fast-DRUMMER/HSV1-Kos.txt',sep = '\t')
    return_df = run_gtest(df)
    return_df.to_csv('/Users/mac/Desktop/Fast-DRUMMER/HSV1-Kos.gtest.txt',sep = '\t',index =None)