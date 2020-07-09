import pandas as pd
import numpy as np
import argparse
import math
from create_output_file import create_output
from scipy.stats import chi2_contingency

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


df = pd.read_csv(input,sep = '\t')

# print(df.columns)
#Get columns of interest
control = [[row['A_unmod'],row['C_unmod'],row['G_unmod'],row['T_unmod'],row['N_unmod']] for index,row in df.iterrows()]
test = [[row['A_mod'],row['C_mod'],row['G_mod'],row['T_mod'],row['N_mod']] for index,row in df.iterrows()]


#Replace 0 with 0.0001 for G-test calculations
control = [ list(map(lambda x: x if x != 0 else 0.0001, i)) for i in control ]
test = [list(map(lambda x: x if x != 0 else 0.0001, i)) for i in test]

#Perform G-test on each row and grab the gtest and pval
lst = [chi2_contingency(np.array([c,t]),lambda_ = 'log-likelihood')[0:2] for c,t in zip(control,test)]

#Save G and pval to variable and assign to column in dataframe
gtests = []
p_vals = []
for g,p in lst:
    gtests.append(round(g))
    p_vals.append(p)

df['G_test'] = gtests
df['p_val'] = p_vals

#Padj is obtained by multiplying pval by length of dataframe
df['padj'] = df['p_val'] * len(df)

output = create_output(output,input,'gTest')
df.to_csv(output,sep = '\t', index = False)