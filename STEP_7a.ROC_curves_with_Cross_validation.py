# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 13:49:21 2018

@author: kncv078
"""
print(__doc__)

import numpy as np
from scipy import interp
import matplotlib.pyplot as plt
#from itertools import cycle
#import os
import pandas as pd

#from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
#from sklearn.model_selection import StratifiedKFold
from sklearn import metrics

#CV_dir = '//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/CrossValidation/crossValidation_RF_{}_analysis_sel8_Assays/'
CV_file = '//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/CrossValidation/crossValidation_RF_{}_analysis_sel8_Assays/Assay_{}/{}-{}_predictions.txt'
sel8_assays_file = '//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/splitter-sel8_assay_list.txt'

# load assay list
sel8_assays = []
with open(sel8_assays_file, 'r')as f:
    for line in f:
        sel8_assays.append(line.strip())

fp_types = {'ECFP4':'ecfp','HTSFP':'htsfp','CESFP':'mesfp'}
my_colors = ['dodgerblue','darkorange','green']
# make plots

mean_fpr = np.linspace(0, 1, 100)
for Assay in sel8_assays:  #os.listdir(ECFP_cv_dir): 
    print(Assay)
    plt.figure(figsize=(14,8))
    j=0
    for fp_name, fp in fp_types.items():
        tprs = []
        aucs = []
        i=0
        Run_no = 0
        for CV_run in range(10):
            Run_no += 1
            df = pd.read_csv(CV_file.format(fp, Assay, Assay, CV_run), sep='\t')
            fpr, tpr, thresholds = metrics.roc_curve(df['label'], df['predprob(A)'], pos_label='A')
            tprs.append(interp(mean_fpr, fpr, tpr))
            tprs[-1][0] = 0.0
            roc_auc = auc(fpr, tpr)
            aucs.append(roc_auc)
            
#            plt.plot(fpr, tpr, lw=1, alpha=0.3, label='ROC fold %d (AUC = %0.2f)' % (i, roc_auc), color=my_colors[j])
            i += 1

        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = auc(mean_fpr, mean_tpr)
        std_auc = np.std(aucs)
        plt.plot(mean_fpr, mean_tpr, color=my_colors[j], label=r'{} mean (AUC: {} $\pm$ {})'.format(fp_name,'%.2f' % mean_auc, round(std_auc,2)), lw=3.0, alpha=.75)
        
        std_tpr = np.std(tprs, axis=0)
        tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
        tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
        plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color=my_colors[j], alpha=.2, label=r'{} $\pm$ 1 std. dev.'.format(fp_name))
        j+=1
    
    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='k', alpha=.8)
    
    plt.xlim([-0.01, 1.0])
    plt.ylim([0, 1.02])
    plt.xlabel('False Positive Rate', size = 22)
    plt.ylabel('True Positive Rate', size = 22)
    plt.title('AID:{}'.format(Assay), size=26)
    plt.legend(loc="lower right", fontsize = 17)
    plt.tick_params(labelsize = 18)
    plt.grid(axis='y')
    plt.savefig('CrossValidation_ROC_plot_s_{}'.format(Assay), bbox_inches='tight')
#    plt.show()
    


