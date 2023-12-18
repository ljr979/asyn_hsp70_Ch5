
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

results_output = f'results/3_mols_per_foci/'
if not os.path.exists(results_output):
    os.makedirs(results_output)
paths = {
    #A8 no ATP
    ('Experiment 75-1', 'A8-noATP-1'): 'D:/Experiments/Experiment_75/python_results/1/',

    ('Experiment 75-2', 'A8-noATP-2'): 'D:/Experiments/Experiment_75/python_results/2/',


    #A8 + ATP
    ('Experiment 72-1', 'A8-ATP-1'): 'D:/Experiments/Experiment_72/python_results/1/',
    #('Experiment 72-2', 'A8-ATP-2'): 'D:/Experiments/Experiment_72/python_results/2/',
    ('Experiment 84-1', 'A8-ATP-3'): 'D:/Experiments/Experiment_84/python_results/py4bleaching/FC1/A8_excess/',


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

files = []
for key, value in paths.items():
    coloc_files = [[f'{root}/{filename}' for filename in files if 'molecule_counts.csv' in filename] for root, dirs, files in os.walk(f'{value}')]


    coloc_files = [
        item for sublist in coloc_files for item in sublist if 'all_combined' not in item]

    df = pd.DataFrame(coloc_files, columns=['path'])
    df[['Experiment_number', 'Treatment-replicate']] = key
    files.append(df)

collated_counts = pd.concat(files)

subunits_foci_all = []
for (path, Treatment), row in collated_counts.groupby(['path', 'Treatment-replicate']):
    path
    path = path.replace('\\', '/')
    row
    replicate_t = Treatment
    replicate = replicate_t.split('-')[-1]

    treatment_from_file = path.split('/')[-3]
    treatment = Treatment.rstrip(replicate_t[-1])
    Exp_num = list(row['Experiment_number'])[0]

    coloc = pd.read_csv(path)
    #coloc['variable1_variable2']=coloc['variable1']+'-'+coloc['variable2']
    fresh_df = pd.DataFrame(coloc['all_small_mol_count'])
    fresh_df['variable1_variable2'] = coloc['variable1']+'-'+coloc['variable2']
    fresh_df['Treatment_from_file'] = treatment_from_file
    fresh_df['Experiment_number'] = Exp_num
    fresh_df['path'] = path
    fresh_df['replicate'] = replicate
    fresh_df['Treatment'] = treatment
    fresh_df['Intensity'] = coloc['max_fluorescence']
    subunits_foci_all.append(fresh_df)

subunits_foci_all = pd.concat(subunits_foci_all)


subunits_foci_all['Experiment_number'].unique().tolist()
manual_filtered=[]
#noatp------------------------------------------------------
t=subunits_foci_all[subunits_foci_all['Experiment_number']=='Experiment 75-1']
manual_filtered.append(t)
t=subunits_foci_all[subunits_foci_all['Experiment_number']=='Experiment 75-2']
manual_filtered.append(t)
#atp------------------------------------------------------
t=subunits_foci_all[subunits_foci_all['Experiment_number']=='Experiment 84-1']
manual_filtered.append(t)
t=subunits_foci_all[subunits_foci_all['Experiment_number']=='Experiment 72-1']
manual_filtered.append(t)

#----------------------------------------------------------
#JB1+A8
t=subunits_foci_all[subunits_foci_all['Experiment_number']=='Experiment 98-1']
#FILTER FOR JUST THE A8 AND JB1 HERE ANDDDDD colocalised v non-colocalised ()
t_filter = t[t['Treatment'] == 'JB1-A8-']
t_filter = t_filter[t_filter['Treatment_from_file'] == 'coloc']
t_filter = t_filter[t_filter['variable1_variable2'] == 'a8-b1-flowout-HSPA8']
manual_filtered.append(t_filter)

#this mol counts is only the + hsp110 ones! need to read in the 'average initial intensity file and concat with the others in subunits foci all then filter again. this portion of code is literally just turning the dataframe into a matching one to those that were output from py4bleaching so i can put them all together later.
test = pd.read_csv(
    'D:/Experiments/Experiment96-FC2-Hsp110-titration/python_results/FC2-1/subunits_per_foci/avg-t0-size-nobleach/avg-intensity-nobleach.csv')
new = pd.DataFrame(test['colocalisation'])
new['all_small_mol_count'] = test['approx_mol_size']
new['variable1_variable2'] = 'a8-b1-flowout-HSPA8'
new['Experiment_number'] = 'Experiment 96-2'
new['path']='D:/Experiments/Experiment96-FC2-Hsp110-titration/python_results/FC2-1/subunits_per_foci/avg-t0-size-nobleach/avg-intensity-nobleach.csv'
new['Treatment_from_file']='avg-initial'
new['Treatment']='JB1-A8-'
new['Intensity']=test['Avg-intensity']
new['replicate']=2
new.to_csv(f'data/3_mols_per_foci/avg_intensity_molsize_E96_2.csv')
#now filter for coloc only, and put into our new list to become a dataframe later.
new.rename(columns={'last_step_mol_size': 'all_small_mol_count'}, inplace=True)
new.drop([col for col in new.columns.tolist() if 'Treatment_from_file' in col], axis=1, inplace=True)
new.rename(columns={'colocalisation': 'Treatment_from_file'}, inplace=True)
t_filter=new[new['Treatment_from_file']=='coloc']
manual_filtered.append(t_filter)



#same situation with experiment 90 FC1 molecules. read in average intensity for a8+b1 situation
test = pd.read_csv(
    'D:/Experiments/Experiment_90-FC1/python_results/subunits_per_foci/avg-nobleach/avg-intensity-nobleach.csv')
new = pd.DataFrame(test['colocalisation'])
new['all_small_mol_count'] = test['approx_mol_size']
new['Intensity'] = test['Avg-intensity']
new['variable1_variable2'] = 'a8-b1-flowout-HSPA8'
new['Experiment_number'] = 'Experiment 90-1'
new['path'] = 'D:/Experiments/Experiment_90-FC1/python_results/subunits_per_foci/avg-nobleach/avg-intensity-nobleach.csv'

#new['Treatment_from_file'] = 'avg-initial'
new['Treatment'] = 'JB1-A8-'
new['replicate'] = 3
new.rename(columns={'colocalisation': 'Treatment_from_file'}, inplace=True)
t_filter = new[new['Treatment_from_file'] == 'Experiment90-FC1-1_HSPA8_A8-B1']
t_filter.to_csv(f'data/3_mols_per_foci/avg_intensity_molsize_E90_1.csv')

manual_filtered.append(t_filter)

t=subunits_foci_all[subunits_foci_all['Experiment_number']=='Experiment 90-2']
t_filter = t[t['Treatment_from_file'] == 'FC2-1']
t_filter['replicate']=4
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-']
manual_filtered.append(t_filter)


#now for exp 86 too lol send help
test=pd.read_csv('D:/Experiments/Experiment86-FC2-bleaches/python_results/FC2/subunits_per_foci/avg-t0-size-nobleach/avg-intensity-nobleach.csv')
new = pd.DataFrame(test['colocalisation'])
new['Intensity'] = test['Avg-intensity']
new['all_small_mol_count'] = test['approx_mol_size']
new['variable1_variable2'] = 'a8-b1-flowout-HSPA8'
new['Experiment_number'] = 'Experiment 86-1'
new['path'] = 'D:/Experiments/Experiment86-FC2-bleaches/python_results/FC2/subunits_per_foci/avg-t0-size-nobleach/avg-intensity-nobleach.csv'

#new['Treatment_from_file'] = 'avg-initial'
new['Treatment'] = 'JB1-A8-'
new['replicate'] = 5
new.rename(columns={'colocalisation': 'Treatment_from_file'}, inplace=True)
new.to_csv(f'data/3_mols_per_foci/avg_intensity_molsize_E86_1.csv')

manual_filtered.append(new)


t = subunits_foci_all[subunits_foci_all['Experiment_number'] == 'Experiment 84-2']

t['FC'] = t.path.str[57:60]
t=t[t['FC']=='FC2']
t=t.drop([col for col in t.columns.tolist() if 'FC' in col], axis=1)
t_filter = t[t['variable1_variable2'] =='hspa8-a8-jb1-excess']
t_filter['replicate'] = 6
t_filter['Treatment'] = 'JB1-A8-'
manual_filtered.append(t_filter)

#------------------------------------
#now onto the next treatment!!!!!!!!!!!!! JB1 + A8 + 110  @ low conc (sub A8)
t = subunits_foci_all[subunits_foci_all['Experiment_number'] == 'Experiment 96-2']
t_filter = t[t['Treatment'] == 'JB1-A8-110-5nM-']
t_filter = t_filter[t_filter['variable1_variable2'] == '+HSP110-5nM-Coloc']
manual_filtered.append(t_filter)


t = subunits_foci_all[subunits_foci_all['Experiment_number'] == 'Experiment 90-2']
t_filter = t[t['Treatment_from_file'] == 'FC2-3']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-110-1nM-']
manual_filtered.append(t_filter)


t = subunits_foci_all[subunits_foci_all['Experiment_number']
                      == 'Experiment 81-1']

t_filter = t[t['Treatment'] == 'JB1-A8-110-1nM-']
manual_filtered.append(t_filter)


t = subunits_foci_all[subunits_foci_all['Experiment_number']
                      == 'Experiment 100-2']

t_filter = t[t['Treatment'] == 'JB1-A8-110-5nM-']
manual_filtered.append(t_filter)
#------------------------------------
#now onto the next treatment!!!!!!!!!!!!! JB1 + A8 + 110  @ middle conc (above A8)

t = subunits_foci_all[subunits_foci_all['Experiment_number']
                      == 'Experiment 100-1']
t_filter = t[t['Treatment_from_file'] == 'coloc']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-110-0.5uM-']
t_filter = t_filter[t_filter['variable1_variable2'] == 'a8-b1-110-HSPA8']
manual_filtered.append(t_filter)


t = subunits_foci_all[subunits_foci_all['Experiment_number']
                      == 'Experiment 90-1']
t_filter = t[t['Treatment_from_file'] == 'FC1-3']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-110-0.5uM-']
manual_filtered.append(t_filter)


t = subunits_foci_all[subunits_foci_all['Experiment_number']
                      == 'Experiment 87-2']
t_filter = t[t['variable1_variable2'] == 'HSPA8-a8-b1+hsp110-2h']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-110-0.5uM-']
manual_filtered.append(t_filter)


t = subunits_foci_all[subunits_foci_all['Experiment_number']
                      == 'Experiment 86-2']
t_filter = t[t['variable1_variable2'] == 't1h-t1h-a8-jb1-110-0.5um']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-110-0.5uM-']
manual_filtered.append(t_filter)

t = subunits_foci_all[subunits_foci_all['Experiment_number'] == 'Experiment 84-2']

t['FC'] = t.path.str[57:60]
t = t[t['FC'] == 'FC2']
t = t.drop([col for col in t.columns.tolist() if 'FC' in col], axis=1)
t_filter = t[t['variable1_variable2'] == 'a8-jb1-110-0.5um-flowout']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-110-0.5uM-']
manual_filtered.append(t_filter)

#-------------------------------------------
#now onto the next treatment!!!!!!!!!!!!! JB1 + A8 + 110  @ high conc (well above A8)

t = subunits_foci_all[subunits_foci_all['Experiment_number']
                      == 'Experiment 96-2']
t_filter = t[t['variable1_variable2'] == '+HSP110-2uM-Coloc']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-110-2uM-']
manual_filtered.append(t_filter)


t = subunits_foci_all[subunits_foci_all['Experiment_number']
                      == 'Experiment 90-1']
t_filter = t[t['variable1_variable2'] == 'a8-b1-2umhsp110-1h-HSPA8']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-110-2uM-']
manual_filtered.append(t_filter)


t = subunits_foci_all[subunits_foci_all['Experiment_number']
                      == 'Experiment 90-2']
t_filter = t[t['variable1_variable2'] == 'a8-10nm-atp-2umhsp110-HSPA8']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-110-2uM-']
manual_filtered.append(t_filter)


#__________________________________________________________
#finally, onto SOD1
t = subunits_foci_all[subunits_foci_all['Experiment_number']
                      == 'Experiment 100-1']
t_filter = t[t['variable1_variable2'] == 'a8-b1-sod1-HSPA8']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-SOD-0.5uM-']
t_filter = t_filter[t_filter['Treatment_from_file'] == 'coloc']
manual_filtered.append(t_filter)


t = subunits_foci_all[subunits_foci_all['Experiment_number']
                      == 'Experiment 92-1']
t_filter = t[t['variable1_variable2'] == 'a8-b1-+sod1-2h-HSPA8']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-SOD-0.5uM-']

manual_filtered.append(t_filter)


t = subunits_foci_all[subunits_foci_all['Experiment_number']
                      == 'Experiment 98-1']
t_filter = t[t['Treatment_from_file'] == 'coloc']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-SOD-0.5uM-']
t_filter = t_filter[t_filter['variable1_variable2'] == 'a8-b1-sod1-HSPA8']

manual_filtered.append(t_filter)


updated_filtered_all=pd.concat(manual_filtered)

data_output = 'data/3_mols_per_foci/'

if not os.path.exists(data_output):
    os.makedirs(data_output)

updated_filtered_all.to_csv(f'{data_output}filter_mols_per_foci.csv')

#little plot how we used to


ax = sns.boxplot(
    x="Treatment",
    y="all_small_mol_count",
    data=updated_filtered_all,
    #order=plot_order,
    palette='Oranges')


# sns.stripplot(x="Treatment",
#               y="all_small_mol_count",
#               data=updated_filtered_all,
#               #order=plot_order,
#               palette='Oranges',
#               jitter=True,
#               alpha=0.1, ax=ax)

ax.set_ylabel('# of subunits')
ax.set_xlabel('Molecule count')
plt.ylim(-2,20)
plt.xticks(rotation=90)
plt.title(f'# of HSPA8 molecules per foci')
plt.tight_layout()
plt.savefig(f'{results_output}boxplot_number_subunits_per_foci.png')
plt.savefig(f'{results_output}boxplot_number_subunits_per_foci.svg')
plt.show()


#and if we log the data?
updated_filtered_all['log_count']=np.log(updated_filtered_all['all_small_mol_count'])
ax = sns.boxplot(
    x="Treatment",
    y="log_count",
    data=updated_filtered_all,
    #order=plot_order,
    palette='Oranges')


# sns.stripplot(x="Treatment",
#               y="all_small_mol_count",
#               data=updated_filtered_all,
#               #order=plot_order,
#               palette='Oranges',
#               jitter=True,
#               alpha=0.1, ax=ax)

ax.set_ylabel('# of subunits')
ax.set_xlabel('Molecule count')
#plt.ylim(-2, 20)
plt.xticks(rotation=90)
plt.title(f'log10 # of HSPA8 molecules per foci')
plt.tight_layout()
plt.savefig(f'{results_output}boxplot_log10_number_subunits_per_foci.png')
plt.savefig(f'{results_output}boxplot_log10_number_subunits_per_foci.svg')
plt.show()

#the above df, is now the data i will work off from here on out. 
#---------------------------------------------------------------------------------
#now i just want to try and do a little thing where I can see what proportion of molecules are 'big' in each treatment
bigs=[]
Experiments=[]
subunits_counts_proportion=[]
threshold=10
for treatment, df in updated_filtered_all.groupby('Treatment'):
    treatment
    df
    for exp, df1 in df.groupby('Experiment_number'):
        exp
        df1
        num_molecules1=len(df1)
        print(exp)
        print(num_molecules1)
        bigs_df1=df1[df1['all_small_mol_count']>threshold]
        num_bigs1=len(bigs_df1)
        proportion_bigs1=num_bigs1/num_molecules1
        print (proportion_bigs1)
        info=[treatment, exp, proportion_bigs1]
        Experiments.append(info)

        num_subunits=sum(df1['all_small_mol_count'])
        num_subunits_bigsdf1=sum(bigs_df1['all_small_mol_count'])
        proportion_subunits_big=num_subunits_bigsdf1/num_subunits
        subunits=[treatment, exp, num_subunits, proportion_subunits_big]
        subunits_counts_proportion.append(subunits)



experiment_split=pd.DataFrame(Experiments, columns=['Treatment', 'Experiment', 'Proportion'])
experiment_split['percentage']=experiment_split['Proportion']*100
experiment_split.to_csv(f'{results_output}threshold{threshold}_proportion-big_per_exp.csv')


subunits_above_thresh = pd.DataFrame(subunits_counts_proportion, columns=['Treatment', 'Experiment', 'total_number_subunits', 'proportion of total subunits in big mols'])
subunits_above_thresh['percentage']=subunits_above_thresh['proportion of total subunits in big mols']*100
subunits_above_thresh.to_csv(f'{results_output}subunits_above_threshold_{threshold}.csv')


melted=pd.melt(experiment_split, id_vars=['Treatment', 'Proportion', 'Experiment'], value_vars=['percentage'], var_name=['stuff'])

melted

Fig,ax=plt.subplots()
ax=sns.barplot(data=melted, x='Treatment', y='value', palette='RdYlGn', saturation=1, capsize=0.3, errcolor='0.4',linewidth=2,edgecolor='0.5',alpha=0.5)
sns.scatterplot(data=melted, x='Treatment', y='value', legend=False, color='k', ax=ax)
Exp_abbrev = subunits_above_thresh['Experiment'].str.split('t').str[-1].tolist()
for i, txt in enumerate(Exp_abbrev):
    ax.annotate(txt, (melted.Treatment[i], melted.value[i]))
ax.set_ylabel('HSPA8 above threshold (%)')
ax.set_xlabel('Treatment')
plt.ylim(0,40)
plt.xticks(rotation=90)
plt.title(f'Proportion of HSPA8 molecules per foci > {threshold}')
plt.tight_layout()
plt.savefig(f'{results_output}percentage_{threshold}-2.svg')
plt.savefig(f'{results_output}percentage_{threshold}-2.pdf')
plt.savefig(f'{results_output}percentage_{threshold}-2.png')
#-------------------------------------------------------------------



#now, I also want to count up the number of TOTAL A8 molecules in each treatment that are colocalised with fibrils, and what proportion of these are sequestered into these large molecules. I've done this in 'subunits above thresh' df. so just want to plot it here. 

Fig, ax = plt.subplots()
ax = sns.barplot(data=subunits_above_thresh, x='Treatment', y='percentage', palette='PuRd', saturation=1, capsize=0.3, errcolor='0.4', linewidth=2, edgecolor='0.5', alpha=0.5)
sns.scatterplot(data=subunits_above_thresh, x='Treatment', y='percentage',
                legend=False, color='k', ax=ax)
Exp_abbrev = subunits_above_thresh['Experiment'].str.split('t').str[-1].tolist()
for i, txt in enumerate(Exp_abbrev):
    ax.annotate(txt, (subunits_above_thresh.Treatment[i], subunits_above_thresh.percentage[i]))
ax.set_ylabel('HSPA8 subunits within molecules above threshold (%)')
ax.set_xlabel('Treatment')
plt.ylim(0,100)
plt.xticks(rotation=90)
plt.title(f'Proportion of HSPA8 molecules within molecules > {threshold}')
plt.tight_layout()
plt.savefig(f'{results_output}percentage_subunits_big_{threshold}.svg')
plt.savefig(f'{results_output}percentage_subunits_big_{threshold}.pdf')
plt.savefig(f'{results_output}percentage_subunits_big_{threshold}.png')


#and now the total # subunits colocalised at each treatment
Fig, ax = plt.subplots()
#test=subunits_above_thresh[subunits_above_thresh['total_number_subunits']<= 10000]
ax = sns.barplot(data=subunits_above_thresh, x='Treatment', y='total_number_subunits', palette='mako',
                 saturation=1, capsize=0.3, errcolor='0.4', linewidth=2, edgecolor='0.5', alpha=0.5)

sns.scatterplot(data=subunits_above_thresh, x='Treatment', y='total_number_subunits',
                legend=False, color='k', ax=ax)
Exp_abbrev = subunits_above_thresh['Experiment'].str.split(
    't').str[-1].tolist()
for i, txt in enumerate(Exp_abbrev):
    ax.annotate(txt, (subunits_above_thresh.Treatment[i], subunits_above_thresh.total_number_subunits[i]))
ax.set_ylabel('colocalised HSPA8 subunits total')
ax.set_xlabel('Treatment')
#plt.ylim(0,45)
plt.xticks(rotation=90)
plt.title(f'Total # HSPA8 molecules colocalised with fibrils')
plt.savefig(f'{results_output}number_total_subunits_coloc.svg')
plt.savefig(f'{results_output}number_total_subunits_coloc.pdf')
plt.savefig(f'{results_output}number_total_subunits_coloc.png')



#-----------------------------------
#now I want to find the fibril density in each experiment, for each treatment. this is going to be compared to the number of SUBUNITS (the raw data un-normalised) on each of these days
foci_fibril_data = pd.read_csv(f'data/2_num_foci/all_foci_per_length.csv')
foci_fibril_data = foci_fibril_data.drop([col for col in foci_fibril_data.columns.tolist() if 'Unnamed: 0' in col], axis=1)

experiment_split=pd.read_csv(f'results/3_mols_per_foci/subunits_above_threshold_8.csv')
experiment_split = experiment_split.drop([col for col in experiment_split.columns.tolist() if 'Unnamed: 0' in col], axis=1)

fibrils = []
for treatment, data in foci_fibril_data.groupby('treatment'):
    treatment
    data
    for e, t in data.groupby('Experiment'):
        e
        t
        num_fibs = len(t)
        info = [treatment, e, num_fibs]
        fibrils.append(info)


fibrils_together = pd.DataFrame(
    fibrils, columns=['Treatment', 'Experiment', 'num_fibs'])
fibrils_together.loc[fibrils_together["Treatment"] ==
                     "JB1-110-0.5uM-", "Treatment"] = 'JB1-A8-110-0.5uM-'
fibrils_together.loc[fibrils_together["Treatment"] ==
                     "JB1-110-2uM-", "Treatment"] = 'JB1-A8-110-2uM-'
#now I am making a temporary column which is the two columns smooshed together so I can map the number of fibrils to these things separately
fibrils_together['temp'] = fibrils_together['Treatment'] + \
    fibrils_together['Experiment']
experiment_split['temp'] = experiment_split['Treatment'] + \
    experiment_split['Experiment']


d = dict(zip(fibrils_together['temp'], fibrils_together['num_fibs']))
experiment_split['num_fibs'] = experiment_split['temp'].map(d)

#now drop the temporary column
fibrils_together = fibrils_together.drop(
    [col for col in fibrils_together.columns.tolist() if 'temp' in col], axis=1)
experiment_split = experiment_split.drop(
    [col for col in experiment_split.columns.tolist() if 'temp' in col], axis=1)


#plot this to see any correlation?
#plotted separately 
for treatment, df in experiment_split.groupby('Treatment'):
    Fig, ax = plt.subplots()
    ax = sns.scatterplot(x='num_fibs', y='percentage', data=df, hue='Treatment', legend='brief', palette='magma')
    plt.ylim(0,100)

Fig, ax= plt.subplots()
ax = sns.scatterplot(x='num_fibs', y='percentage', data=experiment_split, hue='Experiment', legend='brief', palette='viridis')
plt.ylim(20,80)
plt.legend(bbox_to_anchor=(1.25, 1), borderaxespad=0)




#plotting on the same graph
Fig, (ax1, ax2)= plt.subplots(nrows=1, ncols=2, figsize=(16,6))
ax1.set_title ('by_experiment')
sns.scatterplot(x='num_fibs', y='percentage', data=experiment_split, hue='Experiment', legend='brief', palette='Blues', ax=ax1)
ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.0),
          ncol=4, fancybox=True, shadow=True)
ax1.set_ylim(0,100)

ax2.set_title('by_treatment')
sns.scatterplot(x='num_fibs', y='percentage', data=experiment_split, hue='Treatment', legend='brief', palette='Greens', ax=ax2)
ax2.legend(loc='upper center', bbox_to_anchor=(0.5, 1.0),
           ncol=2, fancybox=True, shadow=True)
ax2.set_ylim(0, 100)
plt.savefig('results/3_mols_per_foci/proportion_A8_big8_vs_fibs.png')



