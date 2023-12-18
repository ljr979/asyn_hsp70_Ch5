
from enum import unique
import os, re
import threading
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
import math
results_output = f'results/4_end_middle/'
if not os.path.exists(results_output):
    os.makedirs(results_output)
paths = {
    #A8 no ATP
    ('Experiment 75-1', 'A8-noATP-1', 'old'): 'D:/Experiments/Experiment_75/python_results/1/',

    ('Experiment 75-2', 'A8-noATP-2', 'old'): 'D:/Experiments/Experiment_75/python_results/2/',


    #A8 + ATP
    ('Experiment 72-1', 'A8-ATP-1', 'old'): 'D:/Experiments/Experiment_72/python_results/1/',
    #('Experiment 72-2', 'A8-ATP-2'): 'D:/Experiments/Experiment_72/python_results/2/',
    ('Experiment 84-1', 'A8-ATP-3', 'old'): 'D:/Experiments/Experiment_84/python_results/FC1/',


    #JB1+A8+ATP

    #('Experiment 66-2','JB1-A8-1'):'D:',
    # ('Experiment 73-2','JB1-A8-2'):'D:/Experiments/Experiment_73/python_results/2/',
    # ('Experiment 77-2','JB1-A8-3'):'D:/Experiments/Experiment_77/python_results/2/',
    # ('Experiment 77-1','JB1-A8-4'):'D:/Experiments/Experiment_77/python_results/1/',
    # ('Experiment 76-1','JB1-A8-5'):'D:/Experiments/Experiment_76/python_results/1/',
    ('Experiment 84-2', 'JB1-A8-5', 'new'): 'D: / Experiments/Experiment_84/python_results/',

    ('Experiment 86-2', 'JB1-A8-1', 'new'): 'D:/Experiments/Experiment86-FC2-bleaches/python_results/',
    # ('Experiment 87-2','JB1-A8-7'):'D:/Experiments/Experiment87-FC2/python_results/',
    ('Experiment 90-1', 'JB1-A8-2', 'new'): 'D:/Experiments/Experiment_90-FC1/python_results/',
    ('Experiment 90-2', 'JB1-A8-9', 'new'): 'D:/Experiments/Experiment_90-FC2/python_results/',
    # ('Experiment 85-3','JB1-A8-0'):'D:/Experiments/Experiment85-FC3/python_results/FC3/coords_added/',
    # ('Experiment 92-1', 'JB1-A8-A'): 'D:/Experiments/Experiment_92-fibrils_SOD1_control/python_results/',
    # ('Experiment 96-1', 'JB1-A8-B'): 'D:/Experiments/Experiment96-SOD1-control/python_results/',
    #('Experiment 96-2', 'JB1-A8-3'): 'D:/Experiments/Experiment96-FC2-Hsp110-titration/python_results/FC2-1/subunits_per_foci/avg-t0-size-nobleach/',
    ('Experiment 98-1', 'JB1-A8-4', 'new'): 'D:/Experiments/Experiment98-SOD1/python_results/',


    #MIX DARK A8 2UM & LABELLED A8
    ('Experiment 86-1', 'darkA8-647A8-1', 'new'): 'D:/Experiments/Experiment86-FC1/python_results/',

    ('Experiment 87-1', 'darkA8-647A8-2', 'new'): 'D:/Experiments/Experiment87-FC1/python_results/',


    #JB1+A8+1nM hsp110-excessA8
    # ('Experiment 81-1','JB1-A8-110-1nM-1'):'D:/Experiments/Experiment_81/python_results/1/',
    ('Experiment 81-2', 'JB1-A8-110-1nM-1', 'old'): 'D:/Experiments/Experiment_81/python_results/2/',
    # ('Experiment 84-3','JB1-A8-110-1nM-3'):'D:/Experiments/Experiment_84/python_results/FC3/',
    ('Experiment 90-2', 'JB1-A8-110-1nM-2', 'new'): 'D:/Experiments/Experiment_90-FC2/python_results/',
    ('Experiment 100-2', 'JB1-A8-110-5nM-2', 'new'): 'D:/Experiments/Experiment-100_b/python_results/',
    ('Experiment 96-2', 'JB1-A8-110-5nM-3', 'new'): 'D:/Experiments/Experiment96-FC2-Hsp110-titration/python_results/',

    #JB1+A8+0.5uM hsp110- excess A8
    ('Experiment 81-1', 'JB1-A8-110-1nM-1', 'old'): 'D:/Experiments/Experiment_81/python_results/1/',

    ('Experiment 84-2', 'JB1-A8-110-0.5uM-1', 'new'): 'D:/Experiments/Experiment_84/python_results/',

    ('Experiment 86-2', 'JB1-A8-110-0.5uM-2', 'new'): 'D:/Experiments/Experiment86-FC2-bleaches/python_results/',
    # ('Experiment 90-2','JB1-A8-110-0.5uM-3'):'D:/Experiments/Experiment_90-FC2/python_results/',
    ('Experiment 100-1', 'JB1-A8-110-0.5uM-3', 'new'): 'D:/Experiments/Experiment100_hsp110-SOD-compare/python_results/',
    # ('Experiment 87-2', 'JB1-A8-110-0.5uM-4'): 'D:/Experiments/Experiment87-FC2/python_results/',
    ('Experiment 90-1', 'JB1-A8-110-0.5uM-5', 'new'): 'D:/Experiments/Experiment_90-FC1/python_results/',


    #JB1+A8+2uM hsp110- excess A8
    ('Experiment 90-2', 'JB1-A8-110-2uM-1', 'new'): 'D:/Experiments/Experiment_90-FC2/python_results/',
    ('Experiment 96-2', 'JB1-A8-110-2uM-2', 'new'): 'D:/Experiments/Experiment96-FC2-Hsp110-titration/python_results/',
    ('Experiment 90-1', 'JB1-A8-110-2uM-3', 'new'): 'D:/Experiments/Experiment_90-FC1/python_results/',

    #JB1+A8+0.5uM hsp110- NO excess A8
    # ('Experiment 87-2','JB1-110-0.5uM-1'):'D:/Experiments/Experiment87-FC2/python_results/',
    # ('Experiment 90-1','JB1-110-0.5uM-2'):'D:/Experiments/Experiment_90-FC1/python_results/',


    #JB1+A8+2uM hsp110- NO excess A8
    # ('Experiment 90-1','JB1-110-2uM-1'):'D:/Experiments/Experiment_90-FC1/python_results/',



    #jb1_a8+SOD1@ 0.5uM
    ('Experiment 92-1', 'JB1-A8-SOD-0.5uM-1', 'new'): 'D:/Experiments/Experiment_92-fibrils_SOD1_control/python_results/',
    # ('Experiment 96-1', 'JB1-A8-SOD-0.5uM-2'): 'D:/Experiments/Experiment96-SOD1-control/python_results/',
    ('Experiment 98-1', 'JB1-A8-SOD-0.5uM-3', 'new'): 'D:/Experiments/Experiment98-SOD1/python_results/',
    ('Experiment 100-1', 'JB1-A8-SOD-0.5uM-4', 'new'): 'D:/Experiments/Experiment100_hsp110-SOD-compare/python_results/',


}


files=[]
for key, value in paths.items():
    coloc_files =[[f'{root}/{filename}' for filename in files if 'end_vs_middle.csv' in filename] for root, dirs, files in os.walk(f'{value}')]
    
    coloc_files=[item for sublist in coloc_files for item in sublist if 'step3' not in item]
    
    df=pd.DataFrame(coloc_files, columns=['path'])
    df[['Experiment_number', 'Treatment-replicate', 'analysis-version']]=key
    files.append(df)

collated_collated_colocal=pd.concat(files)
old=collated_collated_colocal[collated_collated_colocal['analysis-version']=='old']
new=collated_collated_colocal[collated_collated_colocal['analysis-version']=='new']
#from collate colocalisation




old_collated=[]
paths_read=[]
for  version, df in collated_collated_colocal.groupby('analysis-version'):
     df
     version

     if version=='old':
        for (path,Treatment), row in df.groupby(['path', 'Treatment-replicate']):
            path
            row
            
            Treatment
            path=path.replace('\\','/')
            replicate=Treatment.split('-')[-1]

            treatment_from_file=path.split('/')[-1].split('_')[-4]
            treatment=Treatment.rstrip(Treatment[-1])
            Exp_num=list(row['Experiment_number'])[0]
            if path not in paths_read:
                print(Exp_num)
                coloc_old=pd.read_csv(path)
                paths_read.append(path)
                #fresh_df=pd.DataFrame(coloc['percent_colocalised'])
                coloc_old['Treatment_from_file']=treatment_from_file
                coloc_old['Experiment_number']=Exp_num
                coloc_old['path']=path
                coloc_old['replicate']=replicate
                coloc_old['Treatment']=treatment
                old_collated.append(coloc_old)
            #protein_coloc=list(coloc['proteins_colocalised'])[0]
old_collated=pd.concat(old_collated)
drops=['Unnamed: 0','Unnamed: 0.1','Class', 'Frame','EndX1, EndX2, EndY1, EndY2', 'Contour_ID']
old_collated.drop([col for col in old_collated.columns.tolist() if col in drops], axis=1, inplace=True)

new_collated=[]

for  version, df in collated_collated_colocal.groupby('analysis-version'):
     df
     version
     if version=='new':

        for (path,Treatment), row in df.groupby(['path', 'Treatment-replicate']):
            path
            


            path=path.replace('\\','/')
            
            replicate_t=Treatment
            replicate=replicate_t.split('-')[-1]

            treatment_from_file=path.split('/')[-1].split('_')[-4]
            treatment=Treatment.rstrip(replicate_t[-1])
            Exp_num=list(row['Experiment_number'])[0]

            row


            coloc_new=pd.read_csv(path)
            paths_read.append(path)
            fresh_df=coloc_new
            
            fresh_df['Treatment_from_file']=treatment_from_file
            fresh_df['Experiment_number']=Exp_num
            fresh_df['path']=path
            fresh_df['replicate']=replicate
            fresh_df['Treatment']=treatment
            
            new_collated.append(fresh_df)

new_collated=pd.concat(new_collated)




new_manual_filtered=[]
old_manual_filtered=[]


t = old_collated[old_collated['Experiment_number']
                      == 'Experiment 75-1']
t_filter = t[t['Treatment_from_file']!= 'collated']
old_manual_filtered.append(t_filter)

t = old_collated[old_collated['Experiment_number']
                      == 'Experiment 75-2']
old_manual_filtered.append(t)
#atp------------------------------------------------------
t = old_collated[old_collated['Experiment_number']
                      == 'Experiment 84-1']
t_filter=t[t['Treatment_from_file']=='84-FC1']
t_filter=t_filter[t_filter['concentration']=='flowout_A8-ATP']
old_manual_filtered.append(t_filter)


t = old_collated[old_collated['Experiment_number']
                      == 'Experiment 72-1']
t_filter=t[t['concentration']=='10nM']
old_manual_filtered.append(t_filter)

#----------------------------------------------------------
#JB1+A8
t = new_collated[new_collated['Experiment_number']
                      == 'Experiment 98-1']
#FILTER FOR JUST THE A8 AND JB1 HERE ANDDDDD colocalised v non-colocalised ()
t_filter = t[t['Treatment'] == 'JB1-A8-']
t_filter = t_filter[t_filter['Treatment_from_file'] == 'Experiment98-A8-B1-flowout']

new_manual_filtered.append(t_filter)

#--------
t=new_collated[new_collated['Experiment_number']== 'Experiment 90-1']
t_filter = t[t['Treatment'] == 'JB1-A8-']
t_filter = t_filter[t_filter['Treatment_from_file']== 'Experiment90-A8-B1-flowout']
new_manual_filtered.append(t_filter)



t = new_collated[new_collated['Experiment_number'] == 'Experiment 90-2']
t_filter = t[t['Treatment'] == 'JB1-A8-']
t_filter = t_filter[t_filter['Treatment_from_file'] == 'Experiment90-A8-B1']
t_filter['replicate'] = 4

new_manual_filtered.append(t_filter)


#now for exp 86 too lol send help
t = new_collated[new_collated['Experiment_number'] == 'Experiment 86-2']
t_filter = t[t['Treatment'] == 'JB1-A8-']
t_filter = t_filter[t_filter['Treatment_from_file'] == 'Experiment86-A8-B1-excess']
t_filter['replicate'] = 5

new_manual_filtered.append(t_filter)



#------------------------------------
#now onto the next treatment!!!!!!!!!!!!! JB1 + A8 + 110  @ low conc (sub A8)
t = new_collated[new_collated['Experiment_number']
                      == 'Experiment 96-2']
t_filter = t[t['Treatment'] == 'JB1-A8-110-5nM-']
t_filter = t_filter[t_filter['Treatment_from_file'] == '+HSP110-5nM']
new_manual_filtered.append(t_filter)


t = new_collated[new_collated['Experiment_number']
                      == 'Experiment 90-2']
t_filter = t[t['Treatment'] == 'JB1-A8-110-1nM-']
t_filter = t_filter[t_filter['Treatment_from_file'] == 'Experiment90-A8-B1-1nMhsp110']
t_filter['Treatment'] = 'JB1-A8-110-5nM-'

new_manual_filtered.append(t_filter)


t = old_collated[old_collated['Experiment_number'] == 'Experiment 81-1']
old_manual_filtered.append(t)


# t = new_collated[new_collated['Experiment_number']
#                       == 'Experiment 100-2']

# t_filter = t[t['Treatment'] == 'JB1-A8-110-5nM-']
# manual_filtered.append(t_filter)
#------------------------------------
#now onto the next treatment!!!!!!!!!!!!! JB1 + A8 + 110  @ middle conc (above A8)

t = new_collated[new_collated['Experiment_number']
                      == 'Experiment 100-1']
t_filter = t[t['Treatment_from_file'] == 'Experiment100-FC1-1_HSPA8_A8-B1-110']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-110-0.5uM-']
new_manual_filtered.append(t_filter)


t = new_collated[new_collated['Experiment_number']
                      == 'Experiment 90-1']
t_filter = t[t['Treatment_from_file'] == 'A8-B1-110']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-110-0.5uM-']
new_manual_filtered.append(t_filter)


t = new_collated[new_collated['Experiment_number']
                      == 'Experiment 86-2']
t_filter = t[t['Treatment_from_file'] == 'A8-B1-0.5HSP110-t1h']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-110-0.5uM-']
new_manual_filtered.append(t_filter)

# t = new_collated[new_collated['Experiment_number']
#                       == 'Experiment 84-2']
# t_filter = t[t['Treatment_from_file'] ==
#              'Experiment84_A8-JB1-110-0.5uM_flowout']
# t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-110-0.5uM-']
# new_manual_filtered.append(t_filter)

#-------------------------------------------
#now onto the next treatment!!!!!!!!!!!!! JB1 + A8 + 110  @ high conc (well above A8)

t = new_collated[new_collated['Experiment_number']
                      == 'Experiment 96-2']
t_filter = t[t['Treatment_from_file'] == '+HSP110-2uM']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-110-2uM-']
new_manual_filtered.append(t_filter)


t = new_collated[new_collated['Experiment_number']
                      == 'Experiment 90-1']
t_filter = t[t['Treatment_from_file'] == 'Experiment90-B1-2uMhsp110']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-110-2uM-']
new_manual_filtered.append(t_filter)


t = new_collated[new_collated['Experiment_number']
                      == 'Experiment 90-2']
t_filter = t[t['Treatment_from_file'] == 'Experiment90-A8-B1-2uMhsp110']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-110-2uM-']
new_manual_filtered.append(t_filter)


#__________________________________________________________
#finally, onto SOD1
t = new_collated[new_collated['Experiment_number']
                      == 'Experiment 100-1']
t_filter = t[t['Treatment_from_file'] == 'Experiment100-FC2-1_HSPA8_A8-B1-SOD1']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-SOD-0.5uM-']

new_manual_filtered.append(t_filter)


t = new_collated[new_collated['Experiment_number']
                      == 'Experiment 92-1']
t_filter = t[t['Treatment_from_file'] == 'A8-B1-SOD1']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-SOD-0.5uM-']

new_manual_filtered.append(t_filter)


t = new_collated[new_collated['Experiment_number']
                      == 'Experiment 98-1']
t_filter = t[t['Treatment_from_file'] == 'Experiment98-A8-B1-SOD1']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-SOD-0.5uM-']


new_manual_filtered.append(t_filter)



#now we want to concat everything in these lists
old_manual_filtered=pd.concat(old_manual_filtered)
new_manual_filtered=pd.concat(new_manual_filtered)

drops=['Unnamed: 0',
 'Unnamed: 0.1',
 'Frame',
 'Contour_ID',
 'EndX1, EndX2, EndY1, EndY2', 
 'Class']
new_manual_filtered.drop([col for col in new_manual_filtered.columns.tolist() if col in drops], axis=1, inplace=True)


collect=old_manual_filtered.append(new_manual_filtered)
collect.to_csv(f'{results_output}filtered_end_vs_middle.csv')