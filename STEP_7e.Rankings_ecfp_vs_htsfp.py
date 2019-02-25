# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 11:02:22 2018

script to compare rankings of active compounds from the ECFP4 vs the HTSFP

@author: kncv078
"""

import pandas as pd
import matplotlib.pyplot as plt
import my_functions_local as mf

#method = 'P_score'
method = 'Ranking'
top_Number = 1000

#input data dir
#pred_files = '//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/MAT2A_test/RF_analysis_MAT2A_{}-fold_CV/MAT2A_{}_CV_predictions.txt'.format(CV_folds,{})
pred_files = '//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/CrossValidation/crossValidation_RF_{}_analysis_sel8_Assays/Assay_{}/{}_predictions_1.txt'
sel8_names_file = '//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/splitter-sel8_assay_list.txt'


sel8_names = mf.file_lines_to_list(sel8_names_file)
sel8_names=[int(x) for x in sel8_names]
sel8_names.sort()


for assay in sel8_names:
    print('doing AID:{}'.format(assay))
    if method == 'Ranking':
        mesfp = pd.read_csv(pred_files.format('mesfp',assay,assay), sep='\t').sort_values('predprob(A)', ascending=False).set_index('cmpd')
        ecfp = pd.read_csv(pred_files.format('ecfp',assay,assay), sep='\t').sort_values('predprob(A)', ascending=False).set_index('cmpd')
        htsfp = pd.read_csv(pred_files.format('htsfp',assay,assay), sep='\t').sort_values('predprob(A)', ascending=False).set_index('cmpd')
        
        mesfp['rank']=range(len(mesfp))
        ecfp['rank']=range(len(ecfp))
        htsfp['rank']=range(len(htsfp))
        
        sel_cmpds = mesfp[:top_Number]
        sel_A = sel_cmpds[sel_cmpds['label']=='A']
        sel_N = sel_cmpds[sel_cmpds['label']=='N']
         
        sC_ecfp_A = ecfp.loc[sel_A.index]
        sC_htsfp_A = htsfp.loc[sel_A.index]
        sC_ecfp_N = ecfp.loc[sel_N.index]
        sC_htsfp_N = htsfp.loc[sel_N.index]
        lowest_rank = max([max(sC_ecfp_A['rank']),max(sC_ecfp_N['rank']),max(sC_htsfp_A['rank']),max(sC_htsfp_N['rank'])])
        
        plt.figure(figsize=(12,8))
        ax = plt.gca()
        plt.title('AID:{}'.format(assay), size=26)
        plt.plot([0,lowest_rank],[0,lowest_rank], color='k', linewidth=2, alpha=0.8)
        plt.plot([1000,lowest_rank+100000],[1000,1000], color='k', linewidth=3, linestyle='--', alpha=0.4)
        plt.plot([1000,1000],[1000,lowest_rank+100000], color='k', linewidth=3, linestyle='--', alpha=0.4)
        plt.scatter(sC_ecfp_N['rank'],sC_htsfp_N['rank'], label='Inactives', color='orange', alpha=0.6)
        plt.scatter(sC_ecfp_A['rank'],sC_htsfp_A['rank'], label='Actives', color='green', alpha=0.85)
        plt.xlim([0.9,lowest_rank+100000])
        plt.ylim([0.9,lowest_rank+100000]) 
        plt.legend(loc="lower right", ncol=1, prop={'size': 22}, frameon=True, facecolor='lightgrey')  
        plt.xlabel('ECFP4 rank', size=24)
        plt.ylabel('HTSFP rank', size=24)
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.tick_params(axis='both', which='major', length=6, width=2.0)
        ax.tick_params(axis='both', which='minor', length=4, width=1.5)
        plt.rc('xtick', labelsize=20)
        plt.rc('ytick', labelsize=20)
        plt.grid(axis='both')
        plt.savefig('Rankings_plot_{}.png'.format(assay), bbox_inches='tight')
        plt.show()
#    break