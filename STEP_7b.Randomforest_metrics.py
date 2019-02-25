# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 15:34:50 2018

@author: kncv078
"""

#import os
import time
import pandas as pd
from sklearn import metrics

start = time.time()

#input data dir
file_path = '//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/CrossValidation/crossValidation_RF_{}_analysis_sel8_Assays/'
sel8_IDs = '//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/splitter-sel8_assay_list.txt'

outdir = '//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/CrossValidation/'


#load list of assays
assay_list=[]
with open(sel8_IDs, 'r') as f:
    for line in f:
        assay_list.append(line.strip())
#print(assay_list)

fp_types = ['htsfp', 'ecfp', 'mesfp']
#colors_list = ['orange','royalblue','forestgreen']
m = 0    
Metric_dict = {}
for fp in fp_types:    
    print('Doing FP: {}'.format(fp))
    for assay in assay_list: 
        print('Assay: {}'.format(assay))
        for j in range(10):
        
            print('CV Fold: {}'.format(j))
            df = pd.read_csv('{}Assay_{}/{}-{}_predictions.txt'.format(file_path.format(fp),assay,assay,j), sep='\t').sort_values('predprob(A)', ascending=False).reset_index(drop=True)
            
            fpr1, tpr1, thresholds1 = metrics.roc_curve(df['label'], df['predprob(A)'], pos_label='A')
            roc_auc1 = metrics.auc(fpr1, tpr1)
            roc_auc1 = round(roc_auc1,3)
            prf1 = metrics.precision_recall_fscore_support(df['label'], df['pred'])
            CM1 = metrics.confusion_matrix(df['label'], df['pred'])
            kappa1 = metrics.cohen_kappa_score(df['label'], df['pred'])
            A_count1 = prf1[3][0]
            N_count1 = prf1[3][1]
            P1 = prf1[0][0]
            R1 = prf1[1][0]
            F11 = prf1[2][0]
            tp1 = CM1[0][0]
            fn1 = CM1[0][1]
            fp1 = CM1[1][0]
            tn1 = CM1[1][1]
            Metric_dict[m] = [fp, assay, j, roc_auc1, kappa1, P1, R1, F11, A_count1, N_count1, tp1, fn1, fp1, tn1]
            m += 1
    
print('saving dataframes containing all scores')
df = pd.DataFrame.from_dict(Metric_dict, orient='index')
df.columns = ['fp', 'assay', 'CV run','roc_auc', 'kappa', 'Precision', 'Recall', 'F1', 'A_count', 'N_count', 'tp', 'fn', 'fp', 'tn']
#dfav = df.mean()
df.to_csv(outdir+'All_metrics.csv', sep='\t')
#dfav.to_csv(outdir+'{}_{}_metrics_average.csv'.format(fp,assay,j), sep='\t')



