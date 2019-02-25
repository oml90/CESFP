#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 17:54:23 2018

@author: kncv078
"""
### NOTES: very inefficient! need more than 64 GB of RAM. otherwise: 'bus error'. possibly fix this.
# solution: use for x, y in zip(list1, list2) ?

import my_functions as mf
    
    
# 2 - define file paths
ecfp4_file   = '/home/kncv078/HTS_Project-Fingerprints/pubchem_dump_analysis/cid2ecfp4.txt'        
htsfp_fp_file  = '/home/kncv078/HTS_Project-Fingerprints/pubchem_dump_analysis/splitter-htsfp_other.txt'  

# In[1] - Load all required files.
print('\nLoading input file:')

# load CID 2 ecfp to Dictionary
print('-ECFP to dict')
ecfp, cmpd_list = mf.TwoColFile_to_dict(ecfp4_file,',')

# load htsfp old. this is used for training and test sets in RF
print('-Other htsfp to binary dict')
binaryhtsFP = mf.flagFP2Binary(htsfp_fp_file)

# In[4] - common cmpds. check to see that compound sets are the same (should be the same)
print('Identifying common cmpds between HTSFP and ECFP')
cmn_cmpds = list(set(binaryhtsFP.keys()) & set(ecfp.keys()))
print('htsFP_cmpds length:\t{}'.format(len(binaryhtsFP.keys()))) 
print('ecfp4_cmpds length:\t{}'.format(len(ecfp.keys())))
print('cmn_cmpds length:\t{}'.format(len(cmn_cmpds)))
if len(cmn_cmpds) != len(cmpd_list): # should not happen # should not happen
    print('WARNING!: compound set of ECFP and HTSFP is NOT the same')

# In[4] - Make new MESFP - combine htsfp with ecfp
print('generating mixed experimental and structural FP (MESFP) - concatinating FPs...\n')
mesFP = {}
for cmpd in cmpd_list:
    mesFP[cmpd] = binaryhtsFP[cmpd]+list(ecfp[cmpd])    

# In[5] - save MESFP to file and binary HTSFP to file
print('saving mesfp to file and saving binaryHTSFP file')
with open('MESFP.txt', 'w') as mesf:
    with open('binaryHTSFP.txt', 'w') as bhtsf:
        for cmpd in cmpd_list:
            mesf.write('{},{}\n'.format(cmpd, ''.join(mesFP[cmpd])))
            bhtsf.write('{},{}\n'.format(cmpd, ''.join(binaryhtsFP[cmpd])))

print('done')

