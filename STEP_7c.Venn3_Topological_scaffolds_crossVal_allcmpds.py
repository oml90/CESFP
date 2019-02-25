# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 13:59:38 2018

@author: kncv078
"""

import os
import time
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold as ms

start = time.time()

topXpercent = 0.2

# In[1] define input data dir
pred_files = '//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/CrossValidation/crossValidation_RF_{}_analysis_sel8_Assays/Assay_{}/{}_predictions_1.txt'

sel8Assays_IDs_file = '//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/splitter-sel8_assay_list.txt'
cid2smi_file = '//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/CID2smi_pubchem.txt'

# In[2] load files
#load list of assays
print('loading assay list..')
assay_list=[]
i = 0
with open(sel8Assays_IDs_file, 'r') as f:
    for line in f:
        i += 1
        assay_list.append(line.strip())

# load dict of Compound ID 2 smiles - required for RDkit to make scaffolds
print('loading smiles file...')
cid2smi = {}
i = 0
with open(cid2smi_file, 'r') as f:
    for line in f:
        i += 1
        cid, smi = line.strip().split('\t')
        cid2smi[cid] = smi

# In[3] load all RF predictions for each of 9 new assays
print('Top scoring {}% of predictions'.format(topXpercent))
i = 1
for assay in assay_list:
    outpath = '//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/CrossValidation/Venn_diagrams_Topological_scafs_top_{}%_allcmpds/'.format(topXpercent)
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    #load files and select top x percent of hits to DF
    print('Doing assay: {}'.format(assay))
    print('\t# g.scafs   /   # CIDs'.format(topXpercent))
    i += 1
    df_ecfp = pd.read_csv(pred_files.format('ecfp',assay,assay), sep='\t').sort_values('predprob(A)', ascending=False).reset_index(drop=True)
    numba = int(len(df_ecfp)*(topXpercent/100))
    TS_ec_df = df_ecfp.iloc[:numba,:].set_index('cmpd')
    del df_ecfp
    df_htsfp = pd.read_csv(pred_files.format('htsfp',assay,assay), sep='\t').sort_values('predprob(A)', ascending=False).reset_index(drop=True)
    TS_hts_df = df_htsfp.iloc[:numba,:].set_index('cmpd')
    del df_htsfp
    df_mesfp = pd.read_csv(pred_files.format('mesfp',assay,assay), sep='\t').sort_values('predprob(A)', ascending=False).reset_index(drop=True)
    TS_mes_df = df_mesfp.iloc[:numba,:].set_index('cmpd')
    del df_mesfp
    
    # remove the CIDs of the top x compounds...
    TPec = list(TS_ec_df.index)
    TPhts = list(TS_hts_df.index)
    TPmes = list(TS_mes_df.index)
    
    # caluculate smiles and Topological scaffold for each compound
    # analysing topological scaffolds...
    cmpd_lists = {'ecfp4':TPec, 'htsfp':TPhts, 'mesfp':TPmes}
    Generic_sets_dict = {}
    for FP_name, cmpds in cmpd_lists.items():
        gen_scaf_set = set()
        for cid in cmpds:
            if str(cid) in cid2smi:                      
                g_scaf = Chem.MolToSmiles(ms.MakeScaffoldGeneric(Chem.MolFromSmiles(cid2smi[str(cid)])))
                gen_scaf_set.add(g_scaf)
            else: 
                print('NA???, thats not meant to happen....')
                continue
        
        print('{}\t   {}     /    {}'.format(FP_name,len(gen_scaf_set),numba))
        Generic_sets_dict[FP_name] = gen_scaf_set
    
    ec_scafs = Generic_sets_dict['ecfp4']
    hts_scafs = Generic_sets_dict['htsfp']
    mes_scafs = Generic_sets_dict['mesfp']

    plt.figure(i, figsize=(6,5))
    plt.title('AID:'+assay, size=15)

    venndiag = venn3([ec_scafs, hts_scafs, mes_scafs], set_labels = ('ecfp4', 'htsfp', 'mesfp'))
    vennC = venn3_circles([ec_scafs, hts_scafs, mes_scafs], linewidth=8)
    vennC[0].set_color('orangered')
    vennC[1].set_color('limegreen')
    vennC[2].set_color('dodgerblue')
    vennC[0].set_alpha(0.3)
    vennC[1].set_alpha(0.3)
    vennC[2].set_alpha(0.3)
    for text in venndiag.set_labels:
        if text == None: continue
        text.set_fontsize(15)
    for text in venndiag.subset_labels:
        if text == None: continue
        text.set_fontsize(10)

    
    plt.savefig('{}{}-{}%_VennDiagram.png'.format(outpath,assay,topXpercent), bbox_inches='tight')
print('done')




