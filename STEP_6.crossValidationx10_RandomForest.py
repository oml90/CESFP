# -*- coding: utf-8 -*-
"""
Created on Fri May 25 16:11:03 2018

@author: kncv078
"""
# NOTES:
# fingerprint files are very large, i should use sparse matrices


# 1 - import libraries

import os
import time
import math
import pandas as pd
import matplotlib as mpl
mpl.use('Agg') # this allows plt functions to work on scp even tough no display possible
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
import matplotlib.pyplot as plt
import my_functions as mf
#from sklearn.model_selection import train_test_split
#import seaborn
#from sklearn.svm import SVC
#from scipy.sparse import coo_matrix

def get_metrics(labels_all, pred_all, Pred_prob_A_all):
    (tp,fn),(fp,tn) = metrics.confusion_matrix(labels_all, pred_all)
    A_count = tp+fn
    N_count = tn+fp
    mcc = ((tp*tn)-(fp*fn))/math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    precision = tp/(tp+fp)
    recall = tp/(tp+fn)
    F1 = (2*tp)/(2*tp+fp+fn)
    kappa = metrics.cohen_kappa_score(labels_all,pred_all)
    print('\n\n#_cmpds:',len(compounds),'\ttotal A:',A_count,'\ttotal N:',N_count)
    print('True Positive:\t', tp)
    print('False Negative:\t', fn)
    print('False Positive:\t', fp)
    print('True Negative:\t', tn)
    print('MCC:\t\t{}'.format(mcc))
    print('Kappa:\t\t{}'.format(kappa))
    if (tp or fn or fp) > 0:
        print('Precision:\t', precision)
        print('Recall:\t\t', recall)
        print('F1 score:\t', F1)
    fpr, tpr, thresholds = metrics.roc_curve(labels_all, Pred_prob_A_all, pos_label='A')
    roc_auc = metrics.auc(fpr, tpr)
    print('ROC AUC:\t', roc_auc)
    return (tp, fn, tn, fp, mcc, precision, recall, F1, kappa, fpr, tpr, roc_auc)

def plot_roc(fpr, tpr, roc_auc, outpath, assay):
    plt.figure()
    plt.plot(fpr, tpr, lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1],[0, 1], lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic')
    plt.legend(loc="lower right")
    plt.savefig('{}{}_ROC_curve_3.png'.format(outpath,assay), bbox_inches='tight')
    plt.show()

start = time.time()
outdir = "CrossValidation/crossValidation_RF_{}_analysis_sel8_Assays/Assay_{}/"

# 2 - define paths
sa8_fp_file         = '/projects/cc/kncv078/MESFP_projects/pubchem_dump_analysis/splitter-htsfp_sel8.txt'
sa8_Assay_IDs_file  = '/projects/cc/kncv078/MESFP_projects/pubchem_dump_analysis/splitter-sel8_assay_list.txt'
ecfp4_file          = '/projects/cc/kncv078/MESFP_projects/pubchem_dump_analysis/cid2ecfp4.txt'        # size: 3.16 GB 
binaryhtsFP_file    = '/projects/cc/kncv078/MESFP_projects/pubchem_dump_analysis/binaryHTSFP.txt'      # size: 2.52 GB
mesfp_file          = '/projects/cc/kncv078/MESFP_projects/pubchem_dump_analysis/MESFP.txt'            # size: 5.65 GB

#sa8_fp_file         = '//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/splitter-htsfp_sel8.txt'
#sa8_Assay_IDs_file  = '//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/splitter-sel8_assay_list.txt'
#ecfp4_file          = '//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/cid2ecfp4.txt'        # size: 3.16 GB 
#binaryhtsFP_file    = '//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/binaryHTSFP.txt'      # size: 2.52 GB
#mesfp_file          = '//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/MESFP.txt'    

# In[1] - Load all required files.
print('\nLoading input file:')
            
# load FP of new 11 assays into dict. these are used as labels in RF
print('-New11 htsfp to df')
            
# load list of headers for newest assays
s8Assays_IDs = mf.file_lines_to_list(sa8_Assay_IDs_file)   

#s8Assays = mf.TwoColFile_to_dict_w_valList(sa8_fp_file,',')
s8Assays = pd.read_csv(sa8_fp_file, sep=',', header= None).set_index(0)
s8Assays = pd.DataFrame(s8Assays[1].apply(list).tolist()).set_index(s8Assays.index)
s8Assays.columns = s8Assays_IDs 

# In[] starting main part
Fingerprints = {"ecfp": ecfp4_file, "htsfp": binaryhtsFP_file, "mesfp": mesfp_file}
for FPname, FP_path in Fingerprints.items():
#    if FPname in ["ecfp","htsfp"]: continue     # TEMPORARY
    print('###'*12)
    print('Beginning analysis of: {}'.format(FPname))
    print('loading {} Fingerprint'.format(FPname))
#    FP = mf.TwoColFile_to_dict_w_valList(FP_path,',')
    FP = pd.read_csv(FP_path, sep=',', header=None).set_index(0)
    FP = pd.DataFrame(FP[1].apply(list).tolist()).set_index(FP.index)
    for assay in s8Assays_IDs:
        print('---'*12)
        print('Assay:',assay)
        
        outpath = outdir.format(FPname, assay)
        
        df_label = s8Assays.loc[s8Assays[assay] != 'x']
        cmpd_idx_list = list(df_label.index)
        # 10 - determine intersection of compound lists newAssay_cmpds and new11Assays
        print('number of compounds in FP:\t\t\t{}'.format(FP.shape[0]))
        print('number of compounds in assay - {}:\t\t{}'.format(assay, len(df_label)))
        
        ac_fp = FP.loc[cmpd_idx_list]
        
        # 14 - check to see if selected new assay has all flags
        x_count = df_label[df_label=='x'].count()
        A_count = df_label[df_label=='A'].count()
        length = df_label[df_label!='x'].count()
        XA2 = pd.concat([x_count, A_count, length],axis=1)
        XA2.columns = ['x_count', 'A_count', 'length']
        print(XA2)
        del x_count
        del A_count
        del XA2
        del length
    
        ## RANDOM FOREST ANALYSIS
        # 16 - define input matrix and labels for sklearn RF classifier
        #df_fp                            # training set        
        df_L = df_label.loc[:,assay]      # labels
        n_rows = len(ac_fp)
        idx_from    = 0
        idx_step    = int(n_rows/10)
        idx_to      = 0
        del df_label
        
        compounds = []
        labels_all = []
        Pred_prob_A_all = []
        pred_all = []
        
        for crossval in range(10):
            print('\nDoing crossvalidation step {} for assay: {} in FP: {}'.format(crossval, assay, FPname))
#            print('test set: rows {} - {}'.format)
            
            idx_to += idx_step
            if idx_to >= n_rows-50: idx_to = None
            fp_test = ac_fp.iloc[idx_from:idx_to]
            fp_train= ac_fp.drop(fp_test.index) 
            L_test  = df_L.iloc[idx_from:idx_to]
            L_train = df_L.drop(L_test.index)
            test_cpmds = list(df_L.iloc[idx_from:idx_to].index)
            print('rows {} - {}'.format(idx_from,idx_to))
            idx_from += idx_step        
            
            
            # 17 - run sklearn RF analysis
            print('Starting random forest')
            rf = RandomForestClassifier(n_estimators=100, class_weight='balanced', max_features='sqrt', min_samples_leaf=10, n_jobs=-1)
            print('fitting...', end='\t')
            rf.fit(fp_train,L_train)
            print('predicting...', end='\t')
            pred = rf.predict(fp_test)
            print('predicting with probability scores...')
            pred_prob = rf.predict_proba(fp_test)
            Feat_imp = rf.feature_importances_
            label = L_test.values
            
            with open('Feature_importance_CESFP_{}_{}.out'.format(assay,crossval), 'w') as f:
                for i in Feat_imp:
                    f.write('{}\n'.format(i))
        
            # 18 - analyse sklearn RF results
            labels_all += list(label)
            Pred_prob_A_all += (list(pred_prob[:,0]))
            pred_all += list(pred)
            compounds += list(test_cpmds)
            
        #save files to outpath directory
        if os.path.exists('{}{}_predictions_3.txt'.format(outpath,assay)) == True:
            print('file with name: {} already exists under this outpath, moving on to next assay'.format(assay))
            continue
        
        #create new dir if not already existing
        if not os.path.exists(outpath): 
            os.makedirs(outpath)
        print('fig saved')
        
        tp, fn, tn, fp, mcc, precision, recall, F1, kappa, fpr, tpr, roc_auc = get_metrics(labels_all, pred_all, Pred_prob_A_all)
#        plot_roc(fpr, tpr, roc_auc, outpath, assay)
        
        pp = pd.DataFrame({'cmpd': compounds,'label':labels_all,'pred':pred_all, 'PredProb(A)': Pred_prob_A_all})
        pp = pp.set_index('cmpd')
        pp = pp.sort_values('PredProb(A)', ascending=False)
#        pp.to_csv('{}{}_predictions_3.txt'.format(outpath,assay), sep='\t')
        print('fig and csv saved\n\n')
    
print('timetaken: {}'.format(time.time()-start))



