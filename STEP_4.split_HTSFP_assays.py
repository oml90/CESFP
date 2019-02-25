# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 13:31:19 2018

@author: kncv078
"""
#import csv
import time
import csv
start = time.time()

htsfp_file = '/projects/cc/kncv078/MESFP_projects/pubchem_dump_analysis/htsfp_gen_out/htsfp_t20000.txt'
htsfp_assays_file = '/projects/cc/kncv078/MESFP_projects/pubchem_dump_analysis/htsfp_gen_out/htsfp_t20000-Assay_list.txt'
ecfp_file = '/projects/cc/kncv078/MESFP_projects/pubchem_dump_analysis/cid2ecfp4.txt'
Assays8_file = '/projects/cc/kncv078/MESFP_projects/pubchem_dump_analysis/sel8_assays.txt'


# In[1] Load all required input files
print('\nLoading all input files:')
#load raw HTSFP into Dictionary
print('-htsfp to dict')
htsfpDict  = {}
k       = 0
with open(htsfp_file, 'r') as fp:    
    for line in fp:
        (key, val) = line.strip().split(',')
        htsfpDict[key] = list(val)
        
        k += 1
        if k % 1000 == 0:
            print('... {}%'.format(round(k/715000*100),0), end='\r')
            
# load ecfp4 cid list
print('\n-ECFP CIDs to list')
ecfp_CIDs = []
k       = 0
with open(ecfp_file, 'r') as f:    
    for line in f:
        Cid = line.strip().split(',')[0]
        ecfp_CIDs.append(Cid)
        
        k += 1
        if k % 1000 == 0:
            print('... {}%'.format(round(k/715000*100),0), end='\r')
    
#load list of Assays
print('\n-Assay list to list')
assays = []
with open(htsfp_assays_file, 'r') as f:
    for entry in f:
        assays.append(entry.strip())
   
#load list of 8 selected assays
print('-8 new assay list...')
Assays8=[]
with open(Assays8_file, 'r') as f:
    for line in f:
        Assays8.append(line.strip())

# In[2] make list of 8 assays index values according to full assay list
print('determining index values of 8 selected assays...')
sa8idx=[]
k=1
for i in Assays8:
    if i in assays:
        idx = assays.index(i)
        sa8idx.append(idx)
    else:
        print('assay:', i, 'not in fpDict')
        continue
sa8idx = sorted(sa8idx, reverse=True)
print('8 Assay index values:',sa8idx)

# In[3] MAIN PART
# split dictionary into new and old assays
print('beginning splitting of FP to old and 8new...')
htsfp_sa8       = {}
htsfp_other = {g: htsfpDict[g] for g in ecfp_CIDs}
del htsfpDict
assay_sa8_list  = []
assay_other_list  = assays.copy()
for i, a in enumerate(sa8idx):
    print('...{}/{}'.format(i+1,len(sa8idx)))
    assay_sa8_list.append(assay_other_list[a])
    del assay_other_list[a]
    for cmpd in ecfp_CIDs:
        if cmpd not in htsfp_sa8:
            htsfp_sa8[cmpd] = [htsfp_other[cmpd][a]]
        else:
            htsfp_sa8[cmpd].append(htsfp_other[cmpd][a])
        del htsfp_other[cmpd][a]
        

# In[4] Save output files 
print('saving output files...')
with open('splitter-other_assay_listQ.txt', 'w') as f:
    for item in assay_other_list:
        f.write(item+'\n')

with open('splitter-sel8_assay_listQ.txt', 'w') as f:
    for item in assay_sa8_list:
        f.write(item+'\n')

with open('splitter-htsfp_sel8Q.txt', 'w') as csv_file:
    writer = csv.writer(csv_file, lineterminator='\n')
    for key, value in htsfp_sa8.items():
       writer.writerow([key,''.join(value)])

with open('splitter-htsfp_otherQ.txt', 'w') as csv_file:
    writer = csv.writer(csv_file, lineterminator='\n')
    for key, value in htsfp_other.items():
       writer.writerow([key,''.join(value)])

print('done\ntime taken:', time.time()-start)
#print('newest 10:',newest10fp)
#print('d2:',old_assays_fp)






