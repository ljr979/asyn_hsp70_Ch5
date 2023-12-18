
from enum import unique
import os
import re
import threading
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
import math


output_folder='data/DNAJB1/2_gather_filter_endmid/'
#defining a dictionary which tells us the experiment number in an easy word as the key, and the PATH to that repo as the value matching the key to call on later

results_output = f'results/DNAJB1/2_gather_filter_endmid/'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
if not os.path.exists(results_output):
    os.makedirs(results_output)

paths = {
    #JB1 ALONE
    ('Experiment 74-2', 'JB1-alone-1', 'one_colour'): 'D:/Experiments/Experiment_64_20220809_DNAJB1titration/python_results/20221221_1/',

    ('Experiment 64-2', 'JB1-alone-2', 'one_colour'): 'D:/Experiments/Experiment64_2/python_results/',
    #for exp 74 it is a subfolder of exp 64 which is weird. have this added on '20221221' or 20221221_1, to access in imagej results folder or python results folder I Think



    #JB1+A8 3 COLOUR
    ('Experiment 84-2', 'JB1-A8-1', 'three_colour'): 'D:/Experiments/Experiment_84/python_results/',

    ('Experiment 85-3', 'JB1-A8-2', 'three_colour'): 'D:/Experiments/Experiment85-FC3/python_results/',
    # ('Experiment 87-2','JB1-A8-7'):'D:/Experiments/Experiment87-FC2/python_results/',
    #('Experiment 66-1', 'JB1-A8-3', 'three_colour'): 'D:/Experiments/Experiment_66/python_results/',


    #jb1_a8+hsp110 @ 0.5uM
    ('Experiment 84-2', 'JB1-A8-110-0.5uM-1', 'three_colour'): 'D:/Experiments/Experiment_84/python_results/',
    ('Experiment 85-3', 'JB1-A8-110-0.5uM-2', 'three_colour'): 'D:/Experiments/Experiment_85-FC3/python_results/',



}
#collate all the % coloc data
#here we collate all the files and make a dataframe with the filepath and the key info we defined above!
#could better automate this by adding a line to sort by old and new BEFORE reading int the data, as this is where it is different. it needs to , for the new ones, loop deeper into the python_results folder to FIND the num_foci and then the treatment folder before looking for percentage colocalisation.for now, hardcoded the folder and can later adjust

#first, gather colocalisation files

files=[]
for key, value in paths.items():
    coloc_files =[[f'{root}/{filename}' for filename in files if 'end_vs_middle.csv' in filename] for root, dirs, files in os.walk(f'{value}')]
    
    coloc_files=[item for sublist in coloc_files for item in sublist if 'step3' not in item]
    coloc_files=[item for item in coloc_files if 'collated' not in item]
    coloc_files=[item for item in coloc_files if 'fibrils-647' not in item]
    df=pd.DataFrame(coloc_files, columns=['path'])
    df[['Experiment_number', 'Treatment-replicate', 'analysis-version']]=key
    files.append(df)

collated_collated_colocal=pd.concat(files)

#from collate colocalisation




alls=[]
paths_read=[]

for (path,Treatment), row in collated_collated_colocal.groupby(['path', 'Treatment-replicate']):
    path
    row
    
    Treatment
    path=path.replace('\\','/')
    replicate=Treatment.split('-')[-1]

    treatment_from_file=path.split('/')[-2]
    treatment=Treatment.rstrip(Treatment[-1])
    Exp_num=list(row['Experiment_number'])[0]

    if path not in paths_read:
        

        coloc_old=pd.read_csv(path)
        paths_read.append(path)
        print(Exp_num, coloc_old.columns.tolist())
        #fresh_df=pd.DataFrame(coloc['percent_colocalised'])
        coloc_old['Treatment_from_file']=treatment_from_file
        coloc_old['Experiment_number']=Exp_num
        coloc_old['path']=path
        coloc_old['replicate']=replicate
        coloc_old['Treatment']=treatment
        alls.append(coloc_old)
    #protein_coloc=list(coloc['proteins_colocalised'])[0]
alls=pd.concat(alls)
drops=['Unnamed: 0','Unnamed: 0.1','Class', 'Frame','EndX1, EndX2, EndY1, EndY2', 'Contour_ID']
alls.drop([col for col in alls.columns.tolist() if col in drops], axis=1, inplace=True)

filtered=[
]
c=alls[alls['Experiment_number']=='Experiment 64-2']
c=c[c['concentration']=='1_JB1-alone-5nMDNAJB1']
c['protein']='JB1'
filtered.append(c)

c=alls[alls['Experiment_number']=='Experiment 74-2']
c=c[c['concentration']=='5nM']
c=c[c['Treatment_from_file']=='DNAJB1']
filtered.append(c)

c=alls[(alls['Experiment_number']=='Experiment 84-2') & (alls['Treatment']=='JB1-A8-')]
c=c[c['concentration']=='A8JB1']
c=c[c['protein']=='JB1']
filtered.append(c)

#----------start here tomorrow
c=alls[(alls['Experiment_number']=='Experiment 85-3') & (alls['Treatment']=='JB1-A8-')]
c=c[c['concentration']=='A8-JB1-excess']
c=c[c['Treatment_from_file']=='dnajb1']
filtered.append(c)



c=alls[(alls['Experiment_number']=='Experiment 84-2') & (alls['Treatment']=='JB1-A8-')]
c=c[c['concentration']=='A8-JB1-110-0.5uM']
c=c[c['protein']=='JB1']
c=c[c['path'].str.contains('flowout')]
c['Treatment']='JB1-A8-110-0.5uM-'
filtered.append(c)

c=alls[(alls['Experiment_number']=='Experiment 85-3') & (alls['Treatment']=='JB1-A8-')]
c=c[c['concentration']=='A8-B1-0.5HSP110-t1h']
c['Treatment']= 'JB1-A8-110-0.5uM-'
c['Treatment_from_file']= 'dnajb1'
filtered.append(c)



filtered=pd.concat(filtered)

filtered.to_csv(f'{output_folder}filtered_end_vs_middle.csv')

