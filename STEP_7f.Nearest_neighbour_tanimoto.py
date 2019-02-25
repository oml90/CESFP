# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 18:28:02 2018
DO pairwise tanimoto similarity for each compound
@author: kncv078
"""

import pandas as pd
import matplotlib.pyplot as plt
#import matplotlib as mpl
#mpl.use('agg') # this allows plt functions to work on scp even though no display possible
import my_functions_local as mf
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs as rdds
import seaborn as sns
import os
os.environ['QT_QPA_PLATFORM']='offscreen'
#plt.style.use('default')

# In[] variables

top_Number = 1000 # number of compounds to use in analysis

# In[] input data directories
#pred_files = '/projects/cc/kncv078/MESFP_projects/pubchem_dump_analysis/CrossValidation/crossValidation_RF_{}_analysis_sel8_Assays/Assay_{}/{}_predictions_1.txt'
#sel8_names_file = '/projects/cc/kncv078/MESFP_projects/pubchem_dump_analysis/splitter-sel8_assay_list.txt'
#cid2smi_file = '/projects/cc/kncv078/MESFP_projects/pubchem_dump_analysis/CID2smi_pubchem.txt'
pred_files = '//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/CrossValidation/crossValidation_RF_{}_analysis_sel8_Assays/Assay_{}/{}_predictions_1.txt'
sel8_names_file = '//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/splitter-sel8_assay_list.txt'
cid2smi_file = '//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/CID2smi_pubchem.txt'

# In[] load list of eight test assays 
sel8_names = mf.file_lines_to_list(sel8_names_file)
sel8_names=[int(x) for x in sel8_names]
sel8_names.sort()

FP_list = ['ecfp', 'htsfp', 'mesfp']
label_dict = {'htsfp':'HTSFP', 'ecfp':'ECFP4', 'mesfp':'CESFP'}
#alpha_dict = {'htsfp':0.9, 'ecfp':0.7, 'mesfp':0.6}
#color_list = {'htsfp':'darkorange', 'ecfp':'royalblue', 'mesfp':'springgreen'}

# In[]
for q, assay in enumerate(sel8_names):
    print('doing AID: {}'.format(assay))
    FP_tanimotos = {}
    plt.figure(q+1, figsize=(12,8))
    plt.title('AID: {}'.format(assay), size=24)

    plt.xlabel('Tanimoto Similarity ', size=22)
    plt.ylabel('Count', size=22)
    plt.grid(alpha=0.5)
    plt.xlim([0,1])
    plt.ylim([0,6])
    
    for fp in FP_list:
        print('doing FP type: {}'.format(fp))
        
        preds = pd.read_csv(pred_files.format(fp,assay,assay), sep='\t').sort_values('predprob(A)', ascending=False).set_index('cmpd')
        m1k = preds[:top_Number] # use all compounds (actives and inactives)
        
#        predsA = preds[preds['label'] == 'A']
#        m1k = predsA[:top_Number] # use only the true positive scaffolds
        
        m1k_cmpds = set(m1k.index)
        
        with open(cid2smi_file, 'r') as f:
            m1k_smiles = []
            i = 0
            for line in f:
                cid, smi = line.strip().split('\t')
                if int(cid) in m1k_cmpds: m1k_smiles.append(smi) ; i+=1
                if i == top_Number: break
        
        # make ECFP4 for all compounds in list    
        ecfp4s = []
        for cmpd1 in m1k_smiles:
            mol1 = Chem.MolFromSmiles(cmpd1)
            if not mol1: print('bad smiles. {}'.format(cmpd1))
            else: 
                ecfp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius=2, nBits=1024) 
                ecfp4s.append(ecfp1)
        
        # calculate Nearest Neighbour tanimoto similarity for the given compound set
        tanimotos = []
        NN = []
        for j, cmpd1 in enumerate(ecfp4s):
            cNN = []
            for cmpd2 in ecfp4s[j+1:]:
                similarity = rdds.FingerprintSimilarity(cmpd1,cmpd2)
                cNN.append(similarity)
            if len(cNN) < 1: break
            NN.append(max(cNN))
                
        FP_tanimotos[fp] = tanimotos
        sns.kdeplot(NN, label=label_dict[fp], bw=.015, linewidth=3)
#        plt.hist(NN, label=label_dict[fp])
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.legend(loc="upper right", ncol=1, prop={'size': 20}, frameon=True, facecolor='gainsboro')  
    plt.savefig('{}_Nearest_Neighbour_Tanimoto_Similarity.png'.format(assay), bbox_inches='tight')
    plt.show()
    
#for fp in FP_list:
#  sns.kdeplot(NN, label=label_dict[fp], bw=.015)
#plt.legend(loc="upper right", ncol=1, prop={'size': 16}, frameon=True, facecolor='gainsboro')  
#plt.figure(q+1, figsize=(12,8))
#plt.title('Nearest Neighbour Tanimoto Similarity {}'.format(assay), size=20)
#plt.xlabel('Tanimoto Similarity ', size=16)
#plt.ylabel('Count', size=16)
#plt.grid(alpha=0.5)
#plt.xlim([0,1])
#plt.ylim([0,6])