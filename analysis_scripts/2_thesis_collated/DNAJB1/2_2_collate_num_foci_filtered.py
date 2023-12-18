
from enum import unique
import os, re
import threading
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
import math

#output_folder='data/0_colocalisation/'
#defining a dictionary which tells us the experiment number in an easy word as the key, and the PATH to that repo as the value matching the key to call on later

results_output = f'results/2_num_foci/'
if not os.path.exists(results_output):
    os.makedirs(results_output)
paths = {
    #A8 no ATP
    ('Experiment 75-1', 'A8-noATP-1'): 'D:/Experiments/Experiment_75/python_results/1/',

    ('Experiment 75-2', 'A8-noATP-2'): 'D:/Experiments/Experiment_75/python_results/2/',


    #A8 + ATP
    ('Experiment 72-1', 'A8-ATP-1'): 'D:/Experiments/Experiment_72/python_results/1/',
    #('Experiment 72-2', 'A8-ATP-2'): 'D:/Experiments/Experiment_72/python_results/2/',
    ('Experiment 84-1', 'A8-ATP-3'): 'D:/Experiments/Experiment_84/python_results/',


    #JB1+A8+ATP

    #('Experiment 66-2','JB1-A8-1'):'D:',
    # ('Experiment 73-2','JB1-A8-2'):'D:/Experiments/Experiment_73/python_results/2/',
    # ('Experiment 77-2','JB1-A8-3'):'D:/Experiments/Experiment_77/python_results/2/',
    # ('Experiment 77-1','JB1-A8-4'):'D:/Experiments/Experiment_77/python_results/1/',
    # ('Experiment 76-1','JB1-A8-5'):'D:/Experiments/Experiment_76/python_results/1/',
    ('Experiment 84-2', 'JB1-A8-5'): 'D: / Experiments/Experiment_84/python_results/',

    ('Experiment 86-2', 'JB1-A8-1'): 'D:/Experiments/Experiment86-FC2-bleaches/python_results/',
    # ('Experiment 87-2','JB1-A8-7'):'D:/Experiments/Experiment87-FC2/python_results/',
    ('Experiment 90-1', 'JB1-A8-2'): 'D:/Experiments/Experiment_90-FC1/python_results/',
    ('Experiment 90-2','JB1-A8-9'):'D:/Experiments/Experiment_90-FC2/python_results/',
    # ('Experiment 85-3','JB1-A8-0'):'D:/Experiments/Experiment85-FC3/python_results/FC3/coords_added/',
    # ('Experiment 92-1', 'JB1-A8-A'): 'D:/Experiments/Experiment_92-fibrils_SOD1_control/python_results/',
    # ('Experiment 96-1', 'JB1-A8-B'): 'D:/Experiments/Experiment96-SOD1-control/python_results/',
    #('Experiment 96-2', 'JB1-A8-3'): 'D:/Experiments/Experiment96-FC2-Hsp110-titration/python_results/FC2-1/subunits_per_foci/avg-t0-size-nobleach/',
    ('Experiment 98-1', 'JB1-A8-4'): 'D:/Experiments/Experiment98-SOD1/python_results/',


    #MIX DARK A8 2UM & LABELLED A8
    ('Experiment 86-1', 'darkA8-647A8-1'): 'D:/Experiments/Experiment86-FC1/python_results/',

    ('Experiment 87-1', 'darkA8-647A8-2'): 'D:/Experiments/Experiment87-FC1/python_results/',


    #JB1+A8+1nM hsp110-excessA8
    # ('Experiment 81-1','JB1-A8-110-1nM-1'):'D:/Experiments/Experiment_81/python_results/1/',
    ('Experiment 81-2', 'JB1-A8-110-1nM-1'): 'D:/Experiments/Experiment_81/python_results/2/',
    # ('Experiment 84-3','JB1-A8-110-1nM-3'):'D:/Experiments/Experiment_84/python_results/FC3/',
    ('Experiment 90-2', 'JB1-A8-110-1nM-2'): 'D:/Experiments/Experiment_90-FC2/python_results/',
    ('Experiment 100-2','JB1-A8-110-5nM-2'):'D:/Experiments/Experiment-100_b/python_results/',
    ('Experiment 96-2', 'JB1-A8-110-5nM-3'): 'D:/Experiments/Experiment96-FC2-Hsp110-titration/python_results/',

    #JB1+A8+0.5uM hsp110- excess A8
    ('Experiment 81-1', 'JB1-A8-110-1nM-1'): 'D:/Experiments/Experiment_81/python_results/1/',

    ('Experiment 84-2','JB1-A8-110-0.5uM-1'):'D:/Experiments/Experiment_84/python_results/',

    ('Experiment 86-2', 'JB1-A8-110-0.5uM-2'): 'D:/Experiments/Experiment86-FC2-bleaches/python_results/',
    # ('Experiment 90-2','JB1-A8-110-0.5uM-3'):'D:/Experiments/Experiment_90-FC2/python_results/',
    ('Experiment 100-1', 'JB1-A8-110-0.5uM-3'): 'D:/Experiments/Experiment100_hsp110-SOD-compare/python_results/',
    # ('Experiment 87-2', 'JB1-A8-110-0.5uM-4'): 'D:/Experiments/Experiment87-FC2/python_results/',
    ('Experiment 90-1', 'JB1-A8-110-0.5uM-5'): 'D:/Experiments/Experiment_90-FC1/python_results/',


    #JB1+A8+2uM hsp110- excess A8
    ('Experiment 90-2', 'JB1-A8-110-2uM-1'): 'D:/Experiments/Experiment_90-FC2/python_results/',
    ('Experiment 96-2', 'JB1-A8-110-2uM-2'): 'D:/Experiments/Experiment96-FC2-Hsp110-titration/python_results/',
    ('Experiment 90-1', 'JB1-A8-110-2uM-3'): 'D:/Experiments/Experiment_90-FC1/python_results/',

    #JB1+A8+0.5uM hsp110- NO excess A8
    # ('Experiment 87-2','JB1-110-0.5uM-1'):'D:/Experiments/Experiment87-FC2/python_results/',
    # ('Experiment 90-1','JB1-110-0.5uM-2'):'D:/Experiments/Experiment_90-FC1/python_results/',


    #JB1+A8+2uM hsp110- NO excess A8
    # ('Experiment 90-1','JB1-110-2uM-1'):'D:/Experiments/Experiment_90-FC1/python_results/',



    #jb1_a8+SOD1@ 0.5uM
    ('Experiment 92-1', 'JB1-A8-SOD-0.5uM-1'): 'D:/Experiments/Experiment_92-fibrils_SOD1_control/python_results/',
    # ('Experiment 96-1', 'JB1-A8-SOD-0.5uM-2'): 'D:/Experiments/Experiment96-SOD1-control/python_results/',
    ('Experiment 98-1', 'JB1-A8-SOD-0.5uM-3'): 'D:/Experiments/Experiment98-SOD1/python_results/',
    ('Experiment 100-1', 'JB1-A8-SOD-0.5uM-4'): 'D:/Experiments/Experiment100_hsp110-SOD-compare/python_results/',


}

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

protein_filt=['HSPA8']
pix_filt=num_foci[num_foci['protein'].isin(protein_filt)]
no_include = ['30nM', '3nM', '1nM', '0.3nM',
              'A8-dark1uM-A8-647-10nM', 'mixed-100min']

filter_treatments=['noATP-10nM', '10nM', 'darkJB1','mixed',
        'Experiment87-A8-B1-excess', '+hsp110-2h',
       'Experiment90-A8-B1-excess', 'Experiment90-A8-B1-flowout',
       'Experiment90-A8-B1-0.5hsp110', 'Experiment90-A8-B1-2uMhsp110',
       'Experiment90-A8-B1', 'Experiment90-A8-extra',
       'Experiment90-A8-B1-1nMhsp110', 'A8-dark1uM-A8-647-10nM',
       'A8-647-flowin-excess', 'HSPA8-JB1-110', 'A8-JB1-110-1nM-excess',
       'A8-JB1-110-1nM-flowout', 'A8JB1', 'A8-JB1-110-0.5uM',
       'Experiment86-A8-B1-excess', 'A8-B1-0.5HSP110-t0',
       'A8-B1-0.5HSP110-t1h']


pix_filt=pix_filt[~pix_filt['concentration'].isin(no_include)]

data_output = 'data/2_num_foci/'

if not os.path.exists(data_output):
         os.makedirs(data_output)

pix_filt.to_csv(f'{data_output}all_foci_per_length.csv')


#now do the manual filter bullshit
manual_filtered=[]
#noatp------------------------------------------------------
t=pix_filt[pix_filt['Experiment']=='Experiment 75-1']
manual_filtered.append(t)
t=pix_filt[pix_filt['Experiment']=='Experiment 75-2']
manual_filtered.append(t)
#atp------------------------------------------------------
t=pix_filt[pix_filt['Experiment']=='Experiment 84-1']
t_filter=t[t['concentration']=='excess_A8']
manual_filtered.append(t_filter)

t=pix_filt[pix_filt['Experiment']=='Experiment 72-1']
manual_filtered.append(t)

#----------------------------------------------------------
#JB1+A8
t=pix_filt[pix_filt['Experiment']=='Experiment 98-1']
#FILTER FOR JUST THE A8 AND JB1 HERE ANDDDDD colocalised v non-colocalised ()
t_filter = t[t['treatment'] == 'JB1-A8-']
t_filter = t_filter[t_filter['concentration'] == 'Experiment98-A8-B1-flowout']
manual_filtered.append(t_filter)

#this mol counts is only the + hsp110 ones! need to read in the 'average initial intensity file and concat with the others in subunits foci all then filter again. this portion of code is literally just turning the dataframe into a matching one to those that were output from py4bleaching so i can put them all together later.
t=pix_filt[pix_filt['Experiment']=='Experiment 90-1']
#FILTER FOR JUST THE A8 AND JB1 HERE ANDDDDD colocalised v non-colocalised ()
t_filter = t[t['treatment'] == 'JB1-A8-']
t_filter = t_filter[t_filter['concentration'] == 'Experiment90-A8-B1-excess']
manual_filtered.append(t_filter)




t=pix_filt[pix_filt['Experiment']=='Experiment 90-2']
t_filter = t[t['treatment'] == 'JB1-A8-']
t_filter = t_filter[t_filter['concentration'] == 'Experiment90-A8-B1']

manual_filtered.append(t_filter)


#now for exp 86 too lol send help


t=pix_filt[pix_filt['Experiment']=='Experiment 86-2']
t_filter = t[t['treatment'] == 'JB1-A8-']
t_filter = t_filter[t_filter['concentration'] == 'Experiment86-A8-B1-excess']

manual_filtered.append(t_filter)

#------------------------------------
#now onto the next treatment!!!!!!!!!!!!! JB1 + A8 + 110  @ low conc (sub A8)
t = pix_filt[pix_filt['Experiment'] == 'Experiment 96-2']
t_filter = t[t['treatment'] == 'JB1-A8-110-5nM-']
t_filter = t_filter[t_filter['concentration'] == '+HSP110-5nM']
manual_filtered.append(t_filter)


t = pix_filt[pix_filt['Experiment'] == 'Experiment 90-2']
t_filter = t[t['concentration'] == 'Experiment90-A8-B1-1nMhsp110']
t_filter = t_filter[t_filter['treatment'] == 'JB1-A8-110-1nM-']
manual_filtered.append(t_filter)


t = pix_filt[pix_filt['Experiment'] == 'Experiment 81-1']

t_filter = t[t['treatment'] == 'JB1-A8-110-1nM-']
manual_filtered.append(t_filter)


t = pix_filt[pix_filt['Experiment']== 'Experiment 100-2']

t_filter = t[t['treatment'] == 'JB1-A8-110-5nM-']
manual_filtered.append(t_filter)
#------------------------------------
#now onto the next treatment!!!!!!!!!!!!! JB1 + A8 + 110  @ middle conc (above A8)

t = pix_filt[pix_filt['Experiment']
                      == 'Experiment 100-1']
t_filter = t[t['concentration'] == 'A8-B1-110']
t_filter = t_filter[t_filter['treatment'] == 'JB1-A8-110-0.5uM-']

manual_filtered.append(t_filter)


t = pix_filt[pix_filt['Experiment']
                      == 'Experiment 90-1']
t_filter = t[t['concentration'] == 'Experiment90-A8-B1-0.5hsp110']
t_filter = t_filter[t_filter['treatment'] == 'JB1-A8-110-0.5uM-']
manual_filtered.append(t_filter)


# t = pix_filt[pix_filt['Experiment']== 'Experiment 87-2']
# t_filter = t[t['variable1_variable2'] == 'HSPA8-a8-b1+hsp110-2h']
# t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-110-0.5uM-']
# manual_filtered.append(t_filter)


t = pix_filt[pix_filt['Experiment']
                      == 'Experiment 86-2']
t_filter = t[t['concentration'] == 'A8-B1-0.5HSP110-t1h']
t_filter = t_filter[t_filter['treatment'] == 'JB1-A8-110-0.5uM-']
manual_filtered.append(t_filter)

t = pix_filt[pix_filt['Experiment'] == 'Experiment 84-2']

t_filter = t[t['concentration'] == 'A8-JB1-110-0.5uM']
t_filter = t_filter[t_filter['treatment'] == 'JB1-A8-110-0.5uM-']
manual_filtered.append(t_filter)

#-------------------------------------------
#now onto the next treatment!!!!!!!!!!!!! JB1 + A8 + 110  @ high conc (well above A8)

t = pix_filt[pix_filt['Experiment']
                      == 'Experiment 96-2']
t_filter = t[t['concentration'] == '+HSP110-2uM']
t_filter = t_filter[t_filter['treatment'] == 'JB1-A8-110-2uM-']
manual_filtered.append(t_filter)


t = pix_filt[pix_filt['Experiment']
                      == 'Experiment 90-1']
t_filter = t[t['concentration'] == 'Experiment90-A8-B1-2uMhsp110']
t_filter = t_filter[t_filter['treatment'] == 'JB1-A8-110-2uM-']
manual_filtered.append(t_filter)


t = pix_filt[pix_filt['Experiment']
                      == 'Experiment 90-2']
t_filter = t[t['concentration'] == 'Experiment90-A8-B1-2uMhsp110']
t_filter = t_filter[t_filter['treatment'] == 'JB1-A8-110-2uM-']
manual_filtered.append(t_filter)


#__________________________________________________________
#finally, onto SOD1
t = pix_filt[pix_filt['Experiment']
                      == 'Experiment 100-1']
t_filter = t[t['concentration'] == 'A8-B1-SOD1']
t_filter = t_filter[t_filter['treatment'] == 'JB1-A8-SOD-0.5uM-']

manual_filtered.append(t_filter)


t = pix_filt[pix_filt['Experiment']
                      == 'Experiment 92-1']
t_filter = t[t['concentration'] == '+SOD1-2h']
t_filter = t_filter[t_filter['treatment'] == 'JB1-A8-SOD-0.5uM-']

manual_filtered.append(t_filter)


t = pix_filt[pix_filt['Experiment']
                      == 'Experiment 98-1']
t_filter = t[t['treatment'] == 'JB1-A8-SOD-0.5uM-']
t_filter = t_filter[t_filter['concentration'] == 'Experiment98-A8-B1-SOD1']

manual_filtered.append(t_filter)


updated_filtered_all=pd.concat(manual_filtered)

updated_filtered_all.to_csv(f'{data_output}filtered_foci_pixel.csv')

#updated_filtered_all['norm_max_exp']=updated_filtered_all['exp_foci']/updated_filtered_all['maximum']

output='results/2_num_foci/'
if not os.path.exists(output):
         os.makedirs(output)


data_output='data/DNAJB1/1_gather_filter_foci/'
updated_filtered_all=pd.read_csv(f'{data_output}filtered_foci_pixel.csv')


#maybe I should be normalising this normalised foci per pixel by the EXPERIMENT i.e. to account for differences in the number of fibrils? 
# or divided by the total number of fibrils available?
norm_per_exp=[]
for experiment, df in updated_filtered_all.groupby('Experiment'):
       maxim=max(df['foci_per_pixel'])
       df['normed']=df['foci_per_pixel']/maxim
       num_fibs=len(df)
       med_foci_100=np.median(df['foci_per_pixel'])*100
       print(f'{experiment, num_fibs, med_foci_100}')
       norm_per_exp.append(df)

new_normed=pd.concat(norm_per_exp)

for experiment, df in new_normed.groupby('Experiment'):
       for treatment, df1 in df.groupby('treatment'):
        num_fibs=len(df1)
        med_foci_100=np.median(df1['foci_per_pixel'])*100
        print(f'{experiment, treatment, num_fibs, med_foci_100}')


new_normed.to_csv(f'{output}filtered_normalised_by_exp.csv')

updated_filtered_all['exp_foci']=np.exp(updated_filtered_all['foci_per_pixel'])
updated_filtered_all['maxi_exp']=max(updated_filtered_all['exp_foci'])
updated_filtered_all['norm_exp_foci']=updated_filtered_all['exp_foci']/updated_filtered_all['maxi_exp']





fig, ax = plt.subplots(figsize=(10,10))
ax=sns.boxplot(x='treatment', y='normed', data=new_normed, palette='PuOr')
plt.title(f'foci per length unit fibril (pixel)')
plt.ylabel('Treatment')
plt.ylabel('normalised (to max) number of foci')
plt.ylim(0,0.8)
plt.xlabel(f'Treatment')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f'{output}/norm_by_exp_foci_per_length_unit_colocalised_boxplot.png')

plt.savefig(f'{output}/norm_by_exp_foci_per_length_unit_colocalised_boxplot.svg')



fig, ax = plt.subplots(figsize=(10,10))
ax=sns.boxplot(x='treatment', y='exp_foci', data=updated_filtered_all, palette='PuOr')
plt.title(f'foci per length unit fibril (pixel)(exponential of foci)')
plt.ylabel('Treatment')
plt.ylabel('(exp) number of foci')
#plt.ylim(0.45,0.7)
plt.xlabel(f'Treatment')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(f'{output}/exp_foci_per_length_unit_colocalised_boxplot.png')

plt.savefig(f'{output}/exp_foci_per_length_unit_colocalised_boxplot.svg')




ax=sns.violinplot(x='treatment', y='norm_foci', data=updated_filtered_all, scale='width', palette='PuOr')
plt.title(f'foci per length unit fibril (pixel)')
plt.ylabel('Treatment')
plt.ylabel('normalised_#foci')
#plt.ylim(1,1.4)
plt.xlabel(f'treatment')
plt.xticks(rotation=45)
plt.tight_layout()

#for thesis :)

conditions_list=['A8-noATP-','A8-ATP-','JB1-A8-','JB1-A8-110-5nM-', 'JB1-A8-110-0.5uM-','JB1-A8-110-2uM-','JB1-A8-SOD-0.5uM-']

select=new_normed[new_normed['treatment'].isin(conditions_list)]
select.to_csv(f'{output}norm_per_experiment_foci_per_pixel.csv')

fig, ax = plt.subplots(figsize=(10,10))
plt.rcParams.update({'font.size': 24})
ax = sns.boxplot(x='treatment', y='normed', data=select,
                 order=conditions_list, palette='RdYlGn')
plt.title(f'Treatment')
plt.ylabel('# foci/pixel fibril length')
plt.ylim(0,1)
plt.xlabel(f'# of foci/pixel')
plt.xticks(rotation=90)

plt.tight_layout()
plt.savefig(f'{output}/norm_exp_filtered_foci_per_length_unit.png')
plt.savefig(f'{output}/norm_exp_filtered_foci_per_length_unit.svg')
plt.show()





#stats





















