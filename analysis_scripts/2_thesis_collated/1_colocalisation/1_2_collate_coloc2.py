
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

#output_folder='data/0_colocalisation/'
#defining a dictionary which tells us the experiment number in an easy word as the key, and the PATH to that repo as the value matching the key to call on later, we also put a word in the tuple to define the 'version' of analysis these are coming from as it will dfeine how we wrangle the data
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


#output_folder='data/0_colocalisation/'
#defining a dictionary which tells us the experiment number in an easy word as the key, and the PATH to that repo as the value matching the key to call on later

results_output = f'results/0_colocalisation/'
if not os.path.exists(results_output):
    os.makedirs(results_output)
paths = {
    #A8 no ATP
    ('Experiment 75-1', 'A8-noATP-1', 'old'): 'D:/Experiments/Experiment_75/python_results/1/',

    ('Experiment 75-2', 'A8-noATP-2', 'old'): 'D:/Experiments/Experiment_75/python_results/2/',


    #A8 + ATP
    ('Experiment 72-1', 'A8-ATP-1', 'old'): 'D:/Experiments/Experiment_72/python_results/1/',
    #('Experiment 72-2', 'A8-ATP-2'): 'D:/Experiments/Experiment_72/python_results/2/',
    ('Experiment 84-1', 'A8-ATP-3', 'old'): 'D:/Experiments/Experiment_84/python_results/',


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
#collate all the % coloc data
#here we collate all the files and make a dataframe with the filepath and the key info we defined above!
#could better automate this by adding a line to sort by old and new BEFORE reading int the data, as this is where it is different. it needs to , for the new ones, loop deeper into the python_results folder to FIND the colocalisation and then the treatment folder before looking for percentage colocalisation.for now, hardcoded the folder and can later adjust
files = []
for key, value in paths.items():
    coloc_files = [[f'{root}/{filename}' for filename in files if 'percentage_colocalisation' in filename]
                   for root, dirs, files in os.walk(f'{value}')]

    coloc_files = [item for sublist in coloc_files for item in sublist]

    df = pd.DataFrame(coloc_files, columns=['path'])
    df[['Experiment_number', 'Treatment-replicate', 'analysis-version']] = key
    files.append(df)

collated_info = pd.concat(files)

colocalisation = []


new_collated = []
old_collated = []
for version, df in collated_info.groupby('analysis-version'):
    df
    version
    if version == 'old':
        for (path, Treatment), row in df.groupby(['path', 'Treatment-replicate']):
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

            if 'proteins_colocalised' not in coloc.columns.tolist():
                coloc['proteins_colocalised'] = 'HSPA8_fibril'

            if 'concentration' not in coloc.columns.tolist():
                coloc['concentration'] = '10nM'

            #fresh_df=pd.DataFrame(coloc['percent_colocalised'])
            coloc['Treatment_from_file'] = treatment_from_file
            coloc['Experiment_number'] = Exp_num
            coloc['path'] = path
            coloc['replicate'] = replicate
            coloc['Treatment'] = treatment
            old_collated.append(coloc)
old_collated = pd.concat(old_collated)
            #protein_coloc=list(coloc['proteins_colocalised'])[0]
for version, df in collated_info.groupby('analysis-version'):
    df
    version
    if version == 'new':
        for (path, Treatment), row in df.groupby(['path', 'Treatment-replicate']):
            path
            path = path.replace('\\', '/')

            replicate_t = Treatment
            replicate = replicate_t.split('-')[-1]

            treatment_from_file = path.split('/')[-3]
            treatment = Treatment.rstrip(replicate_t[-1])
            Exp_num = list(row['Experiment_number'])[0]

            row
            if 'fibrils-647' in path:
                coloc = pd.read_csv(path)

                fresh_df = pd.DataFrame(coloc['percent_colocalised'])
                fresh_df['Treatment_from_file'] = treatment_from_file
                fresh_df['Experiment_number'] = Exp_num
                fresh_df['path'] = path
                fresh_df['replicate'] = replicate
                fresh_df['Treatment'] = treatment
                new_collated.append(fresh_df)

new_collated = pd.concat(new_collated)


#now in old data just need to drop the JB1 data

output = 'data/0_colocalisation/'
if not os.path.exists(output):
      os.makedirs(output)


new_manual_filtered=[]
old_manual_filtered=[]
t = old_collated[old_collated['Experiment_number']
                      == 'Experiment 75-1']
old_manual_filtered.append(t)

t = old_collated[old_collated['Experiment_number']
                      == 'Experiment 75-2']
old_manual_filtered.append(t)
#atp------------------------------------------------------
t = old_collated[old_collated['Experiment_number']
                      == 'Experiment 84-1']
t_filter=t[t['Treatment_from_file']=='FC1']
t_filter=t[t['treatment']=='flowout-A8-ATP']
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
t_filter = t[t['Treatment_from_file'] == 'Experiment90-B1-0.5hsp110']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-110-0.5uM-']
new_manual_filtered.append(t_filter)


t = new_collated[new_collated['Experiment_number']
                      == 'Experiment 86-2']
t_filter = t[t['Treatment_from_file'] == 'Experiment86_A8-B1-0.5HSP110-t1h']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-110-0.5uM-']
new_manual_filtered.append(t_filter)

t = new_collated[new_collated['Experiment_number']
                      == 'Experiment 84-2']
t_filter = t[t['Treatment_from_file'] ==
             'Experiment84_A8-JB1-110-0.5uM_flowout']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-110-0.5uM-']
new_manual_filtered.append(t_filter)

#-------------------------------------------
#now onto the next treatment!!!!!!!!!!!!! JB1 + A8 + 110  @ high conc (well above A8)

t = new_collated[new_collated['Experiment_number']
                      == 'Experiment 96-2']
t_filter = t[t['Treatment_from_file'] == 'Experiment96_+HSP110-2uM']
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
t_filter = t[t['Treatment_from_file'] == 'Experiment92_+SOD1-2h']
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

#------------------
#filter then drop concentration column

old_manual_filtered.drop([col for col in old_manual_filtered.columns.tolist()
                  if 'concentration' in col], axis=1, inplace=True)


#drop 'treatment'
old_manual_filtered.drop([col for col in old_manual_filtered.columns.tolist()
                  if 'treatment' in col], axis=1, inplace=True)
#'drop proteins_colocalised'

old_manual_filtered.drop([col for col in old_manual_filtered.columns.tolist()
                  if 'proteins_colocalised' in col], axis=1, inplace=True)


colocalisation = pd.concat([old_manual_filtered, new_manual_filtered])
colocalisation.drop([col for col in colocalisation.columns.tolist()
                     if 'Unnamed: 0' in col], axis=1, inplace=True)



data_output = 'data/0_colocalisation/'

if not os.path.exists(data_output):
    os.makedirs(data_output)

colocalisation.to_csv(f'{data_output}filter_colocalisation.csv')


#___________________________________________________________________


#now to plot the filtered data

results_output = 'results/0_colocalisation/'
if not os.path.exists(results_output):
    os.makedirs(results_output)



melted=pd.melt(colocalisation, id_vars=['Treatment_from_file', 'Experiment_number', 'path', 'Treatment'], value_vars=['percent_colocalised'], var_name=['stuff'])


Fig, ax = plt.subplots()
ax = sns.barplot(x='Treatment', y='value',
                 data=melted, palette='Purples', alpha=0.45, edgecolor='black')
sns.stripplot(x='Treatment', y='value', data=melted, palette='Purples',
              dodge=True, alpha=0.7, size=8, edgecolor='Purple', linewidth=1, ax=ax)

#little loop to abbreviate replicates with their experiment number
# Exp_abbrev = colocalisation['Experiment_number'].str.split('t').str[-1].tolist()
# for i, txt in enumerate(Exp_abbrev):
#     ax.annotate(txt, (melted.Treatment[i], melted.value[i]))

ax.set_xlabel('HSPA8 treatment')
ax.set_ylabel(f'Percentage of fibrils bound')
ax.set_title('Percentage of fibrils bound by HSPA8')
plt.xticks(rotation=90)
ax.set_ylim(0,100)
plt.tight_layout()
plt.savefig(f'{results_output}/percentage_fibrils_colocalised.png')
plt.savefig(f'{results_output}/percentage_fibrils_colocalised.svg')
plt.show()


#PLOTTING THEM ALL SEPARATELY
conditions_list = ['A8-noATP-']
treat = 'A8_noATP'
condition = colocalisation[colocalisation['Treatment'].isin(conditions_list)]

ax = sns.barplot(x='Treatment', y='percent_colocalised',
                 data=condition, palette='Purples', alpha=0.45, edgecolor='black')
sns.stripplot(x='Treatment', y='percent_colocalised', data=condition, palette='Purples',
              dodge=True, alpha=0.7, size=8, edgecolor='Purple', linewidth=1, ax=ax)
ax.set_xlabel('HSPA8 treatment')
ax.set_ylabel(f'Percentage of fibrils bound')
ax.set_title('Percentage of fibrils bound by HSPA8')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(f'{results_output}/{treat}_percentage_fibrils_colocalised.png')
plt.show()

#+atp
conditions_list = ['A8-noATP-', 'A8-ATP-']
treat = '+A8_ATP'
condition = colocalisation[colocalisation['Treatment'].isin(conditions_list)]
order = []
ax = sns.barplot(x='Treatment', y='percent_colocalised', data=condition,
                 order=conditions_list, palette='Purples', alpha=0.45, edgecolor='black')
sns.stripplot(x='Treatment', y='percent_colocalised', data=condition, order=conditions_list,
              palette='Purples', dodge=True, alpha=0.7, size=8, edgecolor='Purple', linewidth=1, ax=ax)
ax.set_xlabel('HSPA8 treatment')
ax.set_ylabel(f'Percentage of fibrils bound')
ax.set_title('Percentage of fibrils bound by HSPA8')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(f'{results_output}/{treat}_percentage_fibrils_colocalised.png')
plt.show()

#+jb1
conditions_list = ['A8-noATP-', 'A8-ATP-', 'JB1-A8-']
treat = 'b1_a8_atp'
condition = colocalisation[colocalisation['Treatment'].isin(conditions_list)]

ax = sns.barplot(x='Treatment', y='percent_colocalised', data=condition,
                 order=conditions_list, palette='Purples', alpha=0.45, edgecolor='black')
sns.stripplot(x='Treatment', y='percent_colocalised', data=condition, order=conditions_list,
              palette='Purples', dodge=True, alpha=0.7, size=8, edgecolor='Purple', linewidth=1, ax=ax)
ax.set_xlabel('HSPA8 treatment')
ax.set_ylabel(f'Percentage of fibrils bound')
ax.set_title('Percentage of fibrils bound by HSPA8')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(f'{results_output}/{treat}_percentage_fibrils_colocalised.png')
plt.show()

#+110 1nM
conditions_list = ['A8-noATP-', 'A8-ATP-', 'JB1-A8-','JB1-A8-110-1nM-']
treat = 'b1_a8__110_1nm_atp'
condition = colocalisation[colocalisation['Treatment'].isin(conditions_list)]

ax = sns.barplot(x='Treatment', y='percent_colocalised', data=condition,
                 order=conditions_list, palette='Purples', alpha=0.45, edgecolor='black')
sns.stripplot(x='Treatment', y='percent_colocalised', data=condition, order=conditions_list,
              palette='Purples', dodge=True, alpha=0.7, size=8, edgecolor='Purple', linewidth=1, ax=ax)
ax.set_xlabel('HSPA8 treatment')
ax.set_ylabel(f'Percentage of fibrils bound')
ax.set_title('Percentage of fibrils bound by HSPA8')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(f'{results_output}/{treat}_percentage_fibrils_colocalised.png')
plt.show()

#+110 500nM
conditions_list = ['A8-noATP-', 'A8-ATP-', 'JB1-A8-','JB1-A8-110-1nM-', 'JB1-A8-110-0.5uM-']
treat = 'b1_a8__110_0.5nm_atp'
condition = colocalisation[colocalisation['Treatment'].isin(conditions_list)]

ax = sns.barplot(x='Treatment', y='percent_colocalised', data=condition,
                 order=conditions_list, palette='Purples', alpha=0.45, edgecolor='black')
sns.stripplot(x='Treatment', y='percent_colocalised', data=condition, order=conditions_list,
              palette='Purples', dodge=True, alpha=0.7, size=8, edgecolor='Purple', linewidth=1, ax=ax)
ax.set_xlabel('HSPA8 treatment')
ax.set_ylabel(f'Percentage of fibrils bound')
ax.set_title('Percentage of fibrils bound by HSPA8')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(f'{results_output}/{treat}_percentage_fibrils_colocalised.png')
plt.savefig(f'{results_output}/percentage_fibrils_colocalised.svg')
plt.show()

#--------------------------

#for thesis

colocalisation=pd.read_csv(f'data/0_colocalisation/filter_colocalisation.csv')
colocalisation.loc[(colocalisation["Experiment_number"] == 'Experiment 90-2') & (colocalisation["Treatment"] == 'JB1-A8-110-1nM-'), 'Treatment'] = 'JB1-A8-110-5nM-'
conditions_list = ['A8-ATP-', 'JB1-A8-',
                   'JB1-A8-110-5nM-', 'JB1-A8-110-0.5uM-', 'JB1-A8-110-2uM-','JB1-A8-SOD-0.5uM-']

condition = colocalisation[colocalisation['Treatment'].isin(conditions_list)]
fig, ax = plt.subplots(figsize=(10,10))
ax = sns.barplot(x='Treatment', y='percent_colocalised', data=condition,
                 order=conditions_list, palette='RdYlGn', alpha=0.7, capsize=0.3,edgecolor='black')
sns.stripplot(x='Treatment', y='percent_colocalised', data=condition, order=conditions_list,
              palette='RdYlGn', dodge=True, size=8, edgecolor='grey', linewidth=1, ax=ax)
ax.set_xlabel('HSPA8 treatment')
ax.set_ylabel(f'Fibrils bound by HSPA8 (%)')
ax.set_title('Percentage of fibrils bound by HSPA8')
plt.xticks(rotation=90)
plt.ylim(0, 100)
#plt.tight_layout()
plt.savefig(f'{results_output}/percentage_fibrils_colocalised.png')
plt.savefig(f'{results_output}/percentage_fibrils_colocalised.svg')
plt.show()
