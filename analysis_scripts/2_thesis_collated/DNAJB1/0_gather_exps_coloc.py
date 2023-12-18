
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


output_folder='data/DNAJB1/0_gather_filter/'
#defining a dictionary which tells us the experiment number in an easy word as the key, and the PATH to that repo as the value matching the key to call on later

results_output = f'results/DNAJB1/0_gather_filter/'
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
    ('Experiment 66-1', 'JB1-A8-3', 'three_colour'): 'D:/Experiments/Experiment_66/python_results/',


    #jb1_a8+hsp110 @ 0.5uM
    ('Experiment 84-2', 'JB1-A8-110-0.5uM-1', 'three_colour'): 'D:/Experiments/Experiment_84/python_results/',
    ('Experiment 85-3', 'JB1-A8-110-0.5uM-2', 'three_colour'): 'D:/Experiments/Experiment_85-FC3/python_results/',



}
#collate all the % coloc data
#here we collate all the files and make a dataframe with the filepath and the key info we defined above!
#could better automate this by adding a line to sort by old and new BEFORE reading int the data, as this is where it is different. it needs to , for the new ones, loop deeper into the python_results folder to FIND the colocalisation and then the treatment folder before looking for percentage colocalisation.for now, hardcoded the folder and can later adjust

#first, gather colocalisation files
files = []
for key, value in paths.items():
    coloc_files = [[f'{root}/{filename}' for filename in files if 'percentage_colocalisation' in filename]
                   for root, dirs, files in os.walk(f'{value}')]

    coloc_files = [item for sublist in coloc_files for item in sublist]

    df = pd.DataFrame(coloc_files, columns=['path'])
    df[['Experiment_number', 'Treatment-replicate', 'analysis-version']] = key
    files.append(df)

collated_info = pd.concat(files)

colocalisation=[]

for (path, Treatment), row in collated_info.groupby(['path', 'Treatment-replicate']):
    path
    row
    Treatment
    path = path.replace('\\', '/')
    replicate = Treatment.split('-')[-1]

    treatment_from_file = path.split('/')[-3]
    treatment = Treatment.rstrip(Treatment[-1])
    Exp_num = list(row['Experiment_number'])[0]
    coloc = pd.read_csv(path)

    if 'percent_fibs_colocalised' in coloc.columns.tolist():
        coloc.rename(
            columns={'percent_fibs_colocalised': 'percent_colocalised'}, inplace=True)
    if 'percent_hspa8_colocalised' in coloc.columns.tolist():
        coloc.rename(
            columns={'percent_hspa8_colocalised': 'percent_colocalised'}, inplace=True)
    if 'concentration' in coloc.columns.tolist():
        coloc=coloc[coloc['concentration']=='5nM']
    
    coloc=coloc.drop([col for col in coloc.columns.tolist() if 'concentration' in col], axis=1)
    #fresh_df=pd.DataFrame(coloc['percent_colocalised'])
    coloc['Treatment_from_file'] = treatment_from_file
    coloc['Experiment_number'] = Exp_num
    coloc['path'] = path
    coloc['replicate'] = replicate
    coloc['Treatment'] = treatment
    colocalisation.append(coloc)
colocalisation = pd.concat(colocalisation)

filtered=[
]
c=colocalisation[colocalisation['Experiment_number']=='Experiment 64-2']
c=c[c['Treatment_from_file']=='python_results']
c['proteins_colocalised']='JB1_fibril'
filtered.append(c)

c=colocalisation[colocalisation['Experiment_number']=='Experiment 74-2']
c.drop([col for col in c.columns.tolist() if 'treatment' in col], axis=1, inplace=True)
c['proteins_colocalised']='JB1_fibril'
filtered.append(c)

c=colocalisation[(colocalisation['Experiment_number']=='Experiment 84-2') & (colocalisation['Treatment']=='JB1-A8-')]
c=c[c['Treatment_from_file']=='Colocalisation']
c=c[c['path'].str.contains('excess')]
filtered.append(c)


c=colocalisation[(colocalisation['Experiment_number']=='Experiment 85-3') & (colocalisation['Treatment']=='JB1-A8-')]
c=c[c['Treatment_from_file']=='FC3']
c=c[c['path'].str.contains('excess')]
filtered.append(c)

c=colocalisation[(colocalisation['Experiment_number']=='Experiment 84-2')& (colocalisation['Treatment']=='JB1-A8-110-0.5uM-')]

c=c[c['Treatment_from_file']=='Colocalisation']
c=c[c['path'].str.contains('0.5uM')]

filtered.append(c)

c=colocalisation[(colocalisation['Experiment_number']=='Experiment 85-3') & (colocalisation['Treatment']=='JB1-A8-')]
c=c[c['Treatment_from_file']=='FC3']
c=c[c['path'].str.contains('0.5')]
c['Treatment']= 'JB1-A8-110-0.5uM-'
filtered.append(c)



filtered=pd.concat(filtered)

filtered.to_csv(f'{output_folder}filtered_coloc_data.csv')


melted=pd.melt(filtered, id_vars=['Treatment_from_file', 'Experiment_number', 'path', 'Treatment', 'proteins_colocalised'], value_vars=['proteins_colocalised','percent_colocalised'], var_name=['stuff'])


Fig, ax = plt.subplots()
ax = sns.barplot(x='Treatment', y='percent_colocalised',
                 data=filtered, palette='BuPu', hue='proteins_colocalised', alpha=0.45, edgecolor='black')
sns.stripplot(x='Treatment', y='percent_colocalised', data=filtered, palette='BuPu', hue='proteins_colocalised',dodge=True, alpha=0.7, size=8, edgecolor='Black', linewidth=1, ax=ax)

#little loop to abbreviate replicates with their experiment number
# Exp_abbrev = colocalisation['Experiment_number'].str.split('t').str[-1].tolist()
# for i, txt in enumerate(Exp_abbrev):
#     ax.annotate(txt, (melted.Treatment[i], melted.value[i]))
plt.legend(bbox_to_anchor=(1.1, 0), loc='upper left', ncol=1)
ax.set_xlabel('HSPA8 treatment')
ax.set_ylabel(f'Percentage of fibrils bound')
ax.set_title('Percentage of fibrils bound by HSPA8')
plt.xticks(rotation=90)
ax.set_ylim(0,100)
plt.tight_layout()
plt.savefig(f'{results_output}/percentage_fibrils_colocalised.png')
plt.savefig(f'{results_output}/percentage_fibrils_colocalised.svg')
plt.show()