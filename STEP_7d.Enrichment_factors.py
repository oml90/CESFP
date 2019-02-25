# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 14:13:21 2018

@author: kncv078
"""

#import os
import time
import pandas as pd
#import matplotlib as mpl
#mpl.use('Agg') #
#import matplotlib.pyplot as plt
#from sklearn import metrics
#from matplotlib_venn import venn2, venn3
#import seaborn
#from rdkit import Chem
#from rdkit.Chem.Scaffolds import MurckoScaffold
#from scipy.sparse import coo_matrix

start = time.time()
FP_sel = str(0) # 0 or 75
EFx_list = [0.002,0.003,0.005,0.0075,0.01,0.02,0.03,0.05]# enrichment factor proportion

#input data dir
pred_files = '//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/CrossValidation/crossValidation_RF_{}_analysis_sel8_Assays/'
sel8Assays_IDs_file = '//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/splitter-sel8_assay_list.txt'
cid2smi_file = '//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/CID2smi_pubchem.txt'


#load list of assays
print('loading assay list..')
assay_list=[]
with open(sel8Assays_IDs_file, 'r') as f:
    for line in f:
        assay_list.append(line.strip())

# In[3] load all RF predictions for each of 9 new assays

labels1 = ['TP', 'TN', 'FP', 'FN', 'Recall', 'Precision', 'Accuracy', 'HR_all']
labels2 = ['EFx %','lenTopX', 'TP_2', 'FP_2', 'HR_topX', 'HR_all', 'EF']
ecfp_results = pd.DataFrame({'metric' : labels2})
htsfp_results  = pd.DataFrame({'metric' : labels2})
mesfp_results  = pd.DataFrame({'metric' : labels2})

print('beginning loop')
EF_dict = {}
for EFx in EFx_list: 
    print('---'*15)
    print('analysis of EF = {}%'.format(EFx*100))
    for assay in assay_list:
        print('Doing assay: {}'.format(assay))
        FP_list = ['ecfp', 'htsfp', 'mesfp']
        for fp_name in FP_list:
            for CV_run in range(10):
    
                FP_df = pd.read_csv('{}Assay_{}/{}-{}_predictions.txt'.format(pred_files.format(fp_name),assay,assay,CV_run), sep='\t').sort_values('predprob(A)', ascending=False).reset_index(drop=True)
#                df_htsfp = pd.read_csv('{}Assay_{}/{}-{}_predictions.txt'.format(pred_files.format('htsfp'),assay,assay,CV_run), sep='\t').sort_values('predprob(A)', ascending=False).reset_index(drop=True)
#                df_mesfp = pd.read_csv('{}Assay_{}/{}-{}_predictions.txt'.format(pred_files.format('mesfp'),assay,assay,CV_run), sep='\t').sort_values('predprob(A)', ascending=False).reset_index(drop=True)
            
            
                TP = len(FP_df[FP_df['label']=='A']) 
                count_all = len(FP_df)
                HR_all = TP/count_all
                        
                lenTopX = int(count_all*EFx) 
                Topx_df = FP_df.iloc[:lenTopX,:]
                TP_2 = len(Topx_df[Topx_df['label']=='A'])
                HR_topX = TP_2/lenTopX
                EF = HR_topX/HR_all
                
                EF_dict['{}_{}-{}_{}'.format(fp_name, assay, CV_run, (EFx*100))] = [EF, HR_topX, TP_2, lenTopX, HR_all, TP, count_all]
        
    
#print('{} \n {} \n {}'.format(ecfp_df, htsfp_df, mesfp_df))
#ecfp_df = ecfp_df.set_index('metric')
#htsfp_df = htsfp_df.set_index('metric')
#mesfp_df = mesfp_df.set_index('metric')


#ecfp_df.to_csv('/Random_Forest_hts0/Enrichment_ecfp.txt', sep='\t')
#htsfp_df.to_csv('/Random_Forest_hts0/Enrichment_htsfp.txt', sep='\t')
#mesfp_df.to_csv('/Random_Forest_hts0/Enrichment_mesfp.txt', sep='\t')

EF_df = pd.DataFrame.from_dict(EF_dict, orient='index')
EF_df.columns = ['EF', 'HR_topX', 'TP_topX', 'count_topX', 'HR_all', 'TP_all', 'count_all']
EF_averages = pd.DataFrame(columns=['EF', 'HR_topX', 'TP_topX', 'count_topX', 'HR_all', 'TP_all', 'count_all'])
for i in range(0,1920,10):
#    print(i)
    averages = EF_df.iloc[i:i+10,:].mean().values
    fp_type = FP_list[i%3]
    assay = assay_list[i//30%8]
    EFx = EFx_list[i//240]
    EF_averages.loc['{}_{}_{}_mean'.format(fp_type, assay, EFx*100)] = averages
    

EF_df.to_csv('//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/CrossValidation/Enrichment_factors.csv' ,sep='\t')
EF_averages.to_csv('//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/CrossValidation/Enrichment_factor_averages.csv' ,sep='\t')

