
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
#defining a dictionary which tells us the experiment number in an easy word as the key, and the PATH to that repo as the value matching the key to call on later, we also put a word in the tuple to define the 'version' of analysis these are coming from as it will dfeine how we wrangle the data

paths = {
#A8 no ATP
('Experiment 75-1','A8-noATP-1', 'old'):'D:/Experiments/Experiment_75/python_results/1/',

('Experiment 75-2','A8-noATP-2', 'old'):'D:/Experiments/Experiment_75/python_results/2/',

#A8 + ATP
('Experiment 72-1','A8-ATP-1', 'old'):'D:/Experiments/Experiment_72/python_results/1/',
('Experiment 72-2','A8-ATP-2', 'old'):'D:/Experiments/Experiment_72/python_results/2/',
('Experiment 84-1','A8-ATP-3','old'):'D:/Experiments/Experiment_84/python_results/FC1/',


#JB1+A8+ATP
#('Experiment 66-2','JB1-A8-1'):'D:',
('Experiment 73-2','JB1-A8-2', 'old'):'D:/Experiments/Experiment_73/python_results/2/',
('Experiment 77-2','JB1-A8-3', 'old'):'D:/Experiments/Experiment_77/python_results/2/',
('Experiment 77-1','JB1-A8-4', 'old'):'D:/Experiments/Experiment_77/python_results/1/',
('Experiment 76-1','JB1-A8-5', 'old'):'D:/Experiments/Experiment_76/python_results/1/',
('Experiment 76-2','JB1-A8-6', 'old'):'D:/Experiments/Experiment_76/python_results/2/',
('Experiment 87-2','JB1-A8-7','new'):'D:/Experiments/Experiment87-FC2/python_results/Colocalisation/FC2/Experiment87-A8-B1-excess/',
('Experiment 90-1','JB1-A8-8','new'):'D:/Experiments/Experiment_90-FC1/python_results/Colocalisation/Experiment90-A8-B1-flowout/',
('Experiment 90-2','JB1-A8-9','new'):'D:/Experiments/Experiment_90-FC2/python_results/FC2/Colocalisation/Experiment90-A8-B1/',
('Experiment 85-3','JB1-A8-0','new'):'D:/Experiments/Experiment85-FC3/python_results/Colocalisation/FC3/Experiment85_A8-JB1-excess/',
('Experiment 92-1','JB1-A8-A', 'new'):'D:/Experiments/Experiment_92-fibrils_SOD1_control/python_results/Colocalisation/FC1/Experiment92-A8-B1-flowout/',
('Experiment 96-1', 'JB1-A8-B', 'new'): 'D:/Experiments/Experiment96-SOD1-control/python_results/Colocalisation/Experiment96_A8-B1/',
('Experiment 96-2', 'JB1-A8-C', 'new'): 'D:/Experiments/Experiment96-FC2-Hsp110-titration/python_results/Colocalisation/Experiment96-A8-B1-flowout/',
('Experiment 98-1', 'JB1-A8-D', 'new'): 'D:/Experiments/Experiment98-SOD1/python_results/Colocalisation/Experiment98-A8-B1-flowout/',



#MIX DARK A8 2UM & LABELLED A8
('Experiment 86-1','darkA8-647A8-1','new'):'D:/Experiments/Experiment86-FC1/python_results/Colocalisation/FC1/Experiment86_A8-dark1uM-A8-647-10nM/',

('Experiment 87-1','darkA8-647A8-2','new'):'D:/Experiments/Experiment87-FC1/python_results/Colocalisation/FC1/Experiment87_A8-dark1uM-A8-647-10nM/',


#JB1+A8+1nM hsp110-excessA8
('Experiment 81-1','JB1-A8-110-1nM-1', 'old'):'D:/Experiments/Experiment_81/python_results/1/',
('Experiment 81-2','JB1-A8-110-1nM-2', 'old'):'D:/Experiments/Experiment_81/python_results/2/',
('Experiment 84-3','JB1-A8-110-1nM-3','new'):'D:/Experiments/Experiment_84/python_results/FC3/Colocalisation/Experiment84_A8-JB1-110-1nM-flowout/',
('Experiment 90-2','JB1-A8-110-1nM-4','new'):'D:/Experiments/Experiment_90-FC2/python_results/FC2/Colocalisation/Experiment90-A8-B1-1nMhsp110/',
('Experiment 96-2', 'JB1-A8-110-5nM-5', 'new'): 'D:/Experiments/Experiment96-FC2-Hsp110-titration/python_results/Colocalisation/Experiment96_+HSP110-5nM/',
#('Experiment 66-3','JB1-A8-110-1nM-5'):'D:',


#JB1+A8+0.5uM hsp110- excess A8
('Experiment 84-2','JB1-A8-110-0.5uM-1','new'):'D:/Experiments/Experiment_84/python_results/FC2/Colocalisation/Experiment84_A8-JB1-110-0.5uM_flowout/',
('Experiment 86-2','JB1-A8-110-0.5uM-2','new'):'D:/Experiments/Experiment86-FC2-bleaches/python_results/Colocalisation/FC2/Experiment86_A8-B1-0.5HSP110-t1h/',
('Experiment 85-3','JB1-A8-110-0.5uM-3','new'):'D:/Experiments/Experiment85-FC3/python_results/Colocalisation/FC3/Experiment85_A8-B1-0.5HSP110-t1h/',
('Experiment 100-1', 'JB1-A8-110-0.5uM-4', 'new'): 'D:/Experiments/Experiment100_hsp110-SOD-compare/python_results/Colocalisation/Experiment100-FC1-1_HSPA8_A8-B1-110/',


#JB1+A8+2uM hsp110- excess A8
('Experiment 90-2','JB1-A8-110-2uM-1','new'):'D:/Experiments/Experiment_90-FC2/python_results/FC2/Colocalisation/Experiment90-A8-B1-2uMhsp110/',
('Experiment 96-2', 'JB1-A8-110-2uM-2', 'new'): 'D:/Experiments/Experiment96-FC2-Hsp110-titration/python_results/Colocalisation/Experiment96_+HSP110-2uM/',

#JB1+A8+0.5uM hsp110- NO excess A8
('Experiment 87-2','JB1-110-0.5uM-1','new'):'D:/Experiments/Experiment87-FC2/python_results/Colocalisation/FC2/Experiment87_+hsp110-2h/',

('Experiment 90-1','JB1-110-0.5uM-2','new'):'D:/Experiments/Experiment_90-FC1/python_results/Colocalisation/Experiment90-B1-0.5hsp110/',


#JB1+A8+2uM hsp110- NO excess A8
('Experiment 90-1','JB1-110-2uM-1','new'):'D:/Experiments/Experiment_90-FC1/python_results/Colocalisation/Experiment90-B1-2uMhsp110/',


#jb1_a8+SOD1@ 0.5uM 
('Experiment 92-1', 'JB1-A8-SOD-0.5uM-1', 'new'): 'D:/Experiments/Experiment_92-fibrils_SOD1_control/python_results/Colocalisation/FC1/Experiment92_+SOD1-2h/',
('Experiment 96-1', 'JB1-A8-SOD-0.5uM-2', 'new'): 'D:/Experiments/Experiment96-SOD1-control/python_results/Colocalisation/Experiment96_+SOD1-2h/',
('Experiment 98-1', 'JB1-A8-SOD-0.5uM-3', 'new'): 'D:/Experiments/Experiment98-SOD1/python_results/Colocalisation/Experiment98-A8-B1-SOD1/',
('Experiment 100-1', 'JB1-A8-SOD-0.5uM-4', 'new'): 'D:/Experiments/Experiment100_hsp110-SOD-compare/python_results/Colocalisation/Experiment100-FC2-1_HSPA8_A8-B1-SOD1/',

}
#collate all the % coloc data
#here we collate all the files and make a dataframe with the filepath and the key info we defined above!
#could better automate this by adding a line to sort by old and new BEFORE reading int the data, as this is where it is different. it needs to , for the new ones, loop deeper into the python_results folder to FIND the colocalisation and then the treatment folder before looking for percentage colocalisation.for now, hardcoded the folder and can later adjust
files=[]
for key, value in paths.items():
    coloc_files =[[f'{root}/{filename}' for filename in files if 'lengths_data.csv' in filename] for root, dirs, files in os.walk(f'{value}')]
    
    coloc_files=[item for sublist in coloc_files for item in sublist ]
 
    df=pd.DataFrame(coloc_files, columns=['path'])
    df[['Experiment_number', 'Treatment-replicate', 'analysis-version']]=key
    files.append(df)

collated_info=pd.concat(files)

colocalisation=[]


new_collated=[]
old_collated=[]
for  version, df in collated_info.groupby('analysis-version'):
     df
     version
     if version=='old':
        for (path,Treatment), row in df.groupby(['path', 'Treatment-replicate']):
            path
            row
            Treatment
            path=path.replace('\\','/')
            replicate=Treatment.split('-')[-1]

            treatment_from_file=path.split('/')[-3]
            treatment=Treatment.rstrip(Treatment[-1])
            Exp_num=list(row['Experiment_number'])[0]
            coloc=pd.read_csv(path)
            if 'percent_fibs_colocalised' in coloc.columns.tolist():
                coloc.rename(columns={'percent_fibs_colocalised':'percent_colocalised'}, inplace=True)

            if 'proteins_colocalised' not in coloc.columns.tolist():
                coloc['proteins_colocalised']='HSPA8_fibril'

            if 'concentration' not in coloc.columns.tolist():
                coloc['concentration']='10nM'

            #fresh_df=pd.DataFrame(coloc['percent_colocalised'])
            coloc['Treatment_from_file']=treatment_from_file
            coloc['Experiment_number']=Exp_num
            coloc['path']=path
            coloc['replicate']=replicate
            coloc['Treatment']=treatment
            old_collated.append(coloc)
            #protein_coloc=list(coloc['proteins_colocalised'])[0]
            

     if version=='new':
        for (path,Treatment), row in df.groupby(['path', 'Treatment-replicate']):
            path
            path=path.replace('\\','/')

            replicate_t=Treatment
            replicate=replicate_t.split('-')[-1]

            treatment_from_file=path.split('/')[-3]
            treatment=Treatment.rstrip(replicate_t[-1])
            Exp_num=list(row['Experiment_number'])[0]

            row
            if 'fibrils-647' in path:
                coloc=pd.read_csv(path)

                fresh_df=pd.DataFrame(coloc['percent_colocalised'])
                fresh_df['Treatment_from_file']=treatment_from_file
                fresh_df['Experiment_number']=Exp_num
                fresh_df['path']=path
                fresh_df['replicate']=replicate
                fresh_df['Treatment']=treatment
                new_collated.append(fresh_df)

new_collated=pd.concat(new_collated)
old_collated=pd.concat(old_collated)

#now in old data just need to drop the JB1 data
protein_filt=['JB1_fibril']
old_collated=old_collated[~old_collated['proteins_colocalised'].isin(protein_filt)]

#filter then drop concentration column
conc_filt=['10nM','darkJB1', 'noATP-10nM','mixed', 'mixed-100min']

old_collated=old_collated[old_collated['concentration'].isin(conc_filt)]

old_collated.drop([col for col in old_collated.columns.tolist() if 'concentration' in col], axis=1, inplace=True)

#filter old treatment column
conc_filt=['flowout-A8-noATP','flowout-A8-ATP']

old_collated=old_collated[~old_collated['treatment'].isin(conc_filt)]
#drop 'treatment'
old_collated.drop([col for col in old_collated.columns.tolist() if 'treatment' in col], axis=1, inplace=True)
#'drop proteins_colocalised'

old_collated.drop([col for col in old_collated.columns.tolist() if 'proteins_colocalised' in col], axis=1, inplace=True)

          
colocalisation=pd.concat([old_collated, new_collated])



output='data/0_colocalisation/'
if not os.path.exists(output):
         os.makedirs(output)



colocalisation.to_csv(f'{output}/all_colocalisation.csv')


results_output = 'results/0_colocalisation/'
if not os.path.exists(results_output):
    os.makedirs(results_output)
ax=sns.barplot(x='Treatment', y='percent_colocalised', data=colocalisation, palette='Purples', alpha=0.45, edgecolor='black')
sns.stripplot(x='Treatment', y='percent_colocalised', data=colocalisation, palette='Purples', dodge=True, alpha=0.7, size=8, edgecolor='Purple', linewidth=1, ax=ax)
ax.set_xlabel('HSPA8 treatment')
ax.set_ylabel(f'Percentage of fibrils bound')
ax.set_title('Percentage of fibrils bound by HSPA8')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(f'{results_output}/percentage_fibrils_colocalised.png')
plt.savefig(f'{results_output}/percentage_fibrils_colocalised.svg')
plt.show()


#PLOTTING THEM ALL SEPARATELY
conditions_list=['A8-noATP-']
treat='A8_noATP'
condition=colocalisation[colocalisation['Treatment'].isin(conditions_list)]

ax=sns.barplot(x='Treatment', y='percent_colocalised', data=condition, palette='Purples', alpha=0.45, edgecolor='black')
sns.stripplot(x='Treatment', y='percent_colocalised', data=condition, palette='Purples', dodge=True, alpha=0.7, size=8, edgecolor='Purple', linewidth=1, ax=ax)
ax.set_xlabel('HSPA8 treatment')
ax.set_ylabel(f'Percentage of fibrils bound')
ax.set_title('Percentage of fibrils bound by HSPA8')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(f'{results_output}/{treat}_percentage_fibrils_colocalised.png')
plt.show()

#+atp
conditions_list=['A8-noATP-','A8-ATP-']
treat='+A8_ATP'
condition=colocalisation[colocalisation['Treatment'].isin(conditions_list)]
order=[]
ax=sns.barplot(x='Treatment', y='percent_colocalised', data=condition, order=conditions_list, palette='Purples', alpha=0.45, edgecolor='black')
sns.stripplot(x='Treatment', y='percent_colocalised', data=condition, order=conditions_list, palette='Purples', dodge=True, alpha=0.7, size=8, edgecolor='Purple', linewidth=1, ax=ax)
ax.set_xlabel('HSPA8 treatment')
ax.set_ylabel(f'Percentage of fibrils bound')
ax.set_title('Percentage of fibrils bound by HSPA8')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(f'{results_output}/{treat}_percentage_fibrils_colocalised.png')
plt.show()

#+jb1
conditions_list=['A8-noATP-','A8-ATP-','JB1-A8-']
treat='b1_a8_atp'
condition=colocalisation[colocalisation['Treatment'].isin(conditions_list)]

ax=sns.barplot(x='Treatment', y='percent_colocalised', data=condition, order=conditions_list, palette='Purples', alpha=0.45, edgecolor='black')
sns.stripplot(x='Treatment', y='percent_colocalised', data=condition,order=conditions_list,  palette='Purples', dodge=True, alpha=0.7, size=8, edgecolor='Purple', linewidth=1, ax=ax)
ax.set_xlabel('HSPA8 treatment')
ax.set_ylabel(f'Percentage of fibrils bound')
ax.set_title('Percentage of fibrils bound by HSPA8')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(f'{results_output}/{treat}_percentage_fibrils_colocalised.png')
plt.show()

#+110 1nM
conditions_list=['A8-noATP-','A8-ATP-','JB1-A8-','JB1-A8-110-1nM-']
treat='b1_a8__110_1nm_atp'
condition=colocalisation[colocalisation['Treatment'].isin(conditions_list)]

ax=sns.barplot(x='Treatment', y='percent_colocalised', data=condition, order=conditions_list, palette='Purples', alpha=0.45, edgecolor='black')
sns.stripplot(x='Treatment', y='percent_colocalised', data=condition, order=conditions_list, palette='Purples', dodge=True, alpha=0.7, size=8, edgecolor='Purple', linewidth=1, ax=ax)
ax.set_xlabel('HSPA8 treatment')
ax.set_ylabel(f'Percentage of fibrils bound')
ax.set_title('Percentage of fibrils bound by HSPA8')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(f'{results_output}/{treat}_percentage_fibrils_colocalised.png')
plt.show()

#+110 500nM
conditions_list=['A8-noATP-','A8-ATP-','JB1-A8-','JB1-A8-110-1nM-', 'JB1-A8-110-0.5uM-']
treat='b1_a8__110_0.5nm_atp'
condition=colocalisation[colocalisation['Treatment'].isin(conditions_list)]

ax=sns.barplot(x='Treatment', y='percent_colocalised', data=condition, order=conditions_list, palette='Purples', alpha=0.45, edgecolor='black')
sns.stripplot(x='Treatment', y='percent_colocalised', data=condition,order=conditions_list,  palette='Purples', dodge=True, alpha=0.7, size=8, edgecolor='Purple', linewidth=1, ax=ax)
ax.set_xlabel('HSPA8 treatment')
ax.set_ylabel(f'Percentage of fibrils bound')
ax.set_title('Percentage of fibrils bound by HSPA8')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(f'{results_output}/{treat}_percentage_fibrils_colocalised.png')
plt.savefig(f'{results_output}/percentage_fibrils_colocalised.svg')
plt.show()


#FOR EMBO POSTER
results_output = 'D:/presentations_conferences/2023/EMBO/figs_from_fibrils_summary_repo/'

colocalisation.to_csv(f'{results_output}/all_colocalisation.csv')

conditions_list = ['A8-ATP-', 'JB1-A8-',
                   'JB1-A8-110-1nM-', 'JB1-A8-110-0.5uM-']
treat = 'b1_a8__110_0.5nm_atp_for_poster'
condition = colocalisation[colocalisation['Treatment'].isin(conditions_list)]

ax = sns.barplot(x='Treatment', y='percent_colocalised', data=condition,
                 order=conditions_list, palette='RdYlGn', alpha=0.7, edgecolor='black')
sns.stripplot(x='Treatment', y='percent_colocalised', data=condition, order=conditions_list,
              palette='RdYlGn', dodge=True, size=8, edgecolor='grey', linewidth=1, ax=ax)
ax.set_xlabel('HSPA8 treatment')
ax.set_ylabel(f'Percentage of fibrils bound')
ax.set_title('Percentage of fibrils bound by HSPA8')
plt.xticks(rotation=90)
plt.ylim(0,100)
#plt.tight_layout()
plt.savefig(f'{results_output}/{treat}_percentage_fibrils_colocalised.png')
plt.savefig(f'{results_output}/percentage_fibrils_colocalised.svg')
plt.show()
