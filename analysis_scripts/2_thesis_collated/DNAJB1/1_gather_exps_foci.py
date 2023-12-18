
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


output_folder='data/DNAJB1/1_gather_filter_foci/'
#defining a dictionary which tells us the experiment number in an easy word as the key, and the PATH to that repo as the value matching the key to call on later

results_output = f'results/DNAJB1/1_gather_filter_foci/'
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
#could better automate this by adding a line to sort by old and new BEFORE reading int the data, as this is where it is different. it needs to , for the new ones, loop deeper into the python_results folder to FIND the num_foci and then the treatment folder before looking for percentage colocalisation.for now, hardcoded the folder and can later adjust

#first, gather colocalisation files
files = []
for key, value in paths.items():
    coloc_files = [[f'{root}/{filename}' for filename in files if 'foci_per_length_unit.csv' in filename]
                   for root, dirs, files in os.walk(f'{value}')]

    coloc_files = [item for sublist in coloc_files for item in sublist]

    df = pd.DataFrame(coloc_files, columns=['path'])
    df[['Experiment_number', 'Treatment-replicate', 'analysis-version']] = key
    files.append(df)

collated_info = pd.concat(files)


#collate all the % coloc data
num_foci=[]

for key, value in paths.items():
    Experiment=key[0]
    treatment=key[1]
    top_path=f'{value}'
    replicate=treatment.split('-')[-1]
    top_path=f'{value}'

    pix_files =[[f'{root}/{filename}' for filename in files if 'foci_per_length_unit.csv' in filename] for root, dirs, files in os.walk(f'{top_path}')]
    pix_files=[item for sublist in pix_files for item in sublist ]


    for fileo in pix_files: 
        pix= pd.read_csv(f'{fileo}')
        pix['replicate']=replicate
        pix['Experiment']=Experiment
        pix['treatment']=treatment.rstrip(treatment[-1])
        pix['path'] = fileo
        num_foci.append(pix)

num_foci=pd.concat(num_foci)
num_foci.drop([col for col in num_foci.columns.tolist() if 'log_foci' in col], axis=1, inplace=True)


filtered=[
]
c=num_foci[num_foci['Experiment']=='Experiment 64-2']
c=c[c['concentration']=='5nM']
c['protein']='JB1'
filtered.append(c)

c=num_foci[num_foci['Experiment']=='Experiment 74-2']
c=c[c['concentration']=='5nM']

filtered.append(c)

c=num_foci[(num_foci['Experiment']=='Experiment 84-2') & (num_foci['treatment']=='JB1-A8-')]
c=c[c['concentration']=='A8JB1']
c=c[c['protein']=='JB1']
filtered.append(c)

#----------start here tomorrow
c=num_foci[(num_foci['Experiment']=='Experiment 85-3') & (num_foci['treatment']=='JB1-A8-')]
c=c[c['concentration']=='A8-JB1-excess']
c=c[c['protein']=='JB1']
filtered.append(c)



c=num_foci[(num_foci['Experiment']=='Experiment 84-2')& (num_foci['treatment']=='JB1-A8-110-0.5uM-')]
c=c[c['concentration']=='A8-JB1-110-0.5uM']
c=c[c['protein']=='JB1']
filtered.append(c)

c=num_foci[(num_foci['Experiment']=='Experiment 85-3') & (num_foci['treatment']=='JB1-A8-')]
c=c[c['concentration']=='A8-B1-0.5HSP110-t1h']
c['treatment']= 'JB1-A8-110-0.5uM-'
c=c[c['protein']=='JB1']
filtered.append(c)



filtered=pd.concat(filtered)

filtered.to_csv(f'{output_folder}filtered_foci_data.csv')

#plots unnormalised, I think this is more accurate actually.
fig, ax = plt.subplots(figsize=(5,5))
ax=sns.boxplot(x='treatment', y='foci_per_pixel', data=filtered, palette='BuPu')
plt.title(f'foci per length unit fibril (pixel)')
plt.ylabel('Treatment')
plt.ylabel('normalised (to max) number of foci')
plt.ylim(0,0.4)
plt.xlabel(f'Treatment')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f'{results_output}non-norm_JB1_foci.svg')
plt.savefig(f'{results_output}non-norm_JB1_foci.png')

#not normalised between experiments, but log10 them because that's what I did for HSpa8?
filtered['log10_foci_per_pix*100']=np.log10(filtered['foci_per_pixel']*100)
filtered.to_csv(f'{results_output}add_log10_100_foci_perPix.csv')
fig, ax = plt.subplots(figsize=(5,5))
ax=sns.boxplot(x='treatment', y='log10_foci_per_pix*100', data=filtered, palette='BuPu')
plt.title(f'foci per length unit fibril (pixel)')
plt.ylabel('Treatment')
plt.ylabel('log10(Foci per pixel * 100 )')
#plt.ylim(0,0.4)
plt.xlabel(f'Treatment')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f'{results_output}log10_JB1_foci.svg')
plt.savefig(f'{results_output}log10_JB1_foci.png')

#maybe I should be normalising this normalised foci per pixel by the EXPERIMENT i.e. to account for differences in the number of fibrils? 
# or divided by the total number of fibrils available?
norm_per_exp=[]
for experiment, df in filtered.groupby('Experiment'):
       maxim=max(df['foci_per_pixel'])
       df['normed']=df['foci_per_pixel']/maxim
       num_fibs=len(df)
       med_foci_100=np.median(df['foci_per_pixel'])*100
       print(f'{experiment, num_fibs, med_foci_100}')
       norm_per_exp.append(df)

new_normed=pd.concat(norm_per_exp)

for experiment, df in filtered.groupby('Experiment'):
       for treatment, df1 in df.groupby('treatment'):
        num_fibs=len(df1)
        num_available_pixels=sum(df1['lengths (pixel)'])
        med_foci_100=np.median(df1['foci_per_pixel'])*100

        print(f'{experiment, treatment, num_fibs, num_available_pixels, med_foci_100}')


new_normed.to_csv(f'{results_output}filtered_normalised_by_exp.csv')

filtered['exp_foci']=np.exp(filtered['foci_per_pixel'])
filtered['maxi_exp']=max(filtered['exp_foci'])
filtered['norm_exp_foci']=filtered['exp_foci']/filtered['maxi_exp']



fig, ax = plt.subplots(figsize=(10,10))
ax=sns.boxplot(x='treatment', y='normed', data=new_normed, palette='BuPu')
plt.title(f'foci per length unit fibril (pixel)')
plt.ylabel('Treatment')
plt.ylabel('normalised (to max) number of foci')
plt.ylim(0,1)
plt.xlabel(f'Treatment')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f'{results_output}/norm_by_foci_per_length_unit_colocalised_boxplot.png')

plt.savefig(f'{results_output}/norm_by_foci_per_length_unit_colocalised_boxplot.svg')



#trying to normalise to the NUMBER OF AVAILABLE pixels to bind
test=[]
for experiment, df in filtered.groupby('Experiment'):
       for treatment, df1 in df.groupby('treatment'):
        num_fibs=len(df1)
        num_available_pixels=sum(df1['lengths (pixel)'])
        med_foci_100=np.median(df1['foci_per_pixel'])*100
        df1['norm_crazy']= df1['foci_per_pixel']/num_available_pixels
        df1['cbrt_norm_crazy']=np.cbrt(df1['norm_crazy'])
        print(f'{experiment, treatment, num_fibs, num_available_pixels, med_foci_100}')
        test.append(df1)


test=pd.concat(test)


#trying to normalise to the NUMBER OF AVAILABLE pixels to bind
test=[]
for experiment, df in filtered.groupby('Experiment'):
       for treatment, df1 in df.groupby('treatment'):
        num_fibs=len(df1)
        num_available_pixels=sum(df1['lengths (pixel)'])
        med_foci_100=np.median(df1['foci_per_pixel'])*100
        df1['norm_crazy']= df1['foci_per_pixel']*num_available_pixels
        #df1['ccrazymulitple_norm_crazy']=df1['norm_crazy']*100000
        df1['log_norm_crazy']=np.log10(df1['norm_crazy'])
        print(f'{experiment, treatment, num_fibs, num_available_pixels, med_foci_100}')
        test.append(df1)


test=pd.concat(test)