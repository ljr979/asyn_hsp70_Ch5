
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
paths = {
#A8 no ATP
('Experiment 75-1','A8-noATP-1'):'D:/Experiments/Experiment_75/python_results/1/',

('Experiment 75-2','A8-noATP-2'):'D:/Experiments/Experiment_75/python_results/2/',


#A8 + ATP
('Experiment 72-1','A8-ATP-1'):'D:/Experiments/Experiment_72/python_results/1/',
('Experiment 72-2','A8-ATP-2'):'D:/Experiments/Experiment_72/python_results/2/',
('Experiment 84-1','A8-ATP-3'):'D:/Experiments/Experiment_84/python_results/py4bleaching/FC1/A8_excess/',


#JB1+A8+ATP

#('Experiment 66-2','JB1-A8-1'):'D:',
# ('Experiment 73-2','JB1-A8-2'):'D:/Experiments/Experiment_73/python_results/2/',
# ('Experiment 77-2','JB1-A8-3'):'D:/Experiments/Experiment_77/python_results/2/',
# ('Experiment 77-1','JB1-A8-4'):'D:/Experiments/Experiment_77/python_results/1/',
# ('Experiment 76-1','JB1-A8-5'):'D:/Experiments/Experiment_76/python_results/1/',
('Experiment 86-2','JB1-A8-1'):'D:/Experiments/Experiment86-FC2-bleaches/python_results/',
# ('Experiment 87-2','JB1-A8-7'):'D:/Experiments/Experiment87-FC2/python_results/',
('Experiment 90-1','JB1-A8-2'):'D:/Experiments/Experiment_90-FC1/python_results/',
# ('Experiment 90-2','JB1-A8-9'):'D:/Experiments/Experiment_90-FC2/python_results/',
# ('Experiment 85-3','JB1-A8-0'):'D:/Experiments/Experiment85-FC3/python_results/FC3/coords_added/',
# ('Experiment 92-1', 'JB1-A8-A'): 'D:/Experiments/Experiment_92-fibrils_SOD1_control/python_results/',
# ('Experiment 96-1', 'JB1-A8-B'): 'D:/Experiments/Experiment96-SOD1-control/python_results/',
('Experiment 96-2', 'JB1-A8-3'): 'D:/Experiments/Experiment96-FC2-Hsp110-titration/python_results/',
('Experiment 98-1', 'JB1-A8-4'): 'D:/Experiments/Experiment98-SOD1/python_results/',


#MIX DARK A8 2UM & LABELLED A8
('Experiment 86-1','darkA8-647A8-1'):'D:/Experiments/Experiment86-FC1/python_results/',

('Experiment 87-1','darkA8-647A8-2'):'D:/Experiments/Experiment87-FC1/python_results/',


#JB1+A8+1nM hsp110-excessA8
# ('Experiment 81-1','JB1-A8-110-1nM-1'):'D:/Experiments/Experiment_81/python_results/1/',
 ('Experiment 81-2','JB1-A8-110-1nM-1'):'D:/Experiments/Experiment_81/python_results/2/',
# ('Experiment 84-3','JB1-A8-110-1nM-3'):'D:/Experiments/Experiment_84/python_results/FC3/',
('Experiment 90-2','JB1-A8-110-1nM-2'):'D:/Experiments/Experiment_90-FC2/python_results/',
#('Experiment 66-3','JB1-A8-110-1nM-5'):'D:',
('Experiment 96-2', 'JB1-A8-110-5nM-3'): 'D:/Experiments/Experiment96-FC2-Hsp110-titration/python_results/',

#JB1+A8+0.5uM hsp110- excess A8
    ('Experiment 81-1', 'JB1-A8-110-1nM-1'): 'D:/Experiments/Experiment_81/python_results/1/',
# ('Experiment 84-2','JB1-A8-110-0.5uM-1'):'D:/Experiments/Experiment_84/python_results/py4bleaching/FC2/A8-JB1-110-0.5uM_flowout/',
('Experiment 86-2','JB1-A8-110-0.5uM-2'):'D:/Experiments/Experiment86-FC2-bleaches/python_results/',
# ('Experiment 90-2','JB1-A8-110-0.5uM-3'):'D:/Experiments/Experiment_90-FC2/python_results/',
('Experiment 100-1', 'JB1-A8-110-0.5uM-3'): 'D:/Experiments/Experiment100_hsp110-SOD-compare/python_results/',
('Experiment 87-2','JB1-A8-110-0.5uM-4'):'D:/Experiments/Experiment87-FC2/python_results/',
('Experiment 90-1','JB1-A8-110-0.5uM-5'):'D:/Experiments/Experiment_90-FC1/python_results/',


#JB1+A8+2uM hsp110- excess A8
('Experiment 90-2','JB1-A8-110-2uM-1'):'D:/Experiments/Experiment_90-FC2/python_results/',
('Experiment 96-2', 'JB1-A8-110-2uM-2'): 'D:/Experiments/Experiment96-FC2-Hsp110-titration/python_results/',
('Experiment 90-1','JB1-A8-110-2uM-3'):'D:/Experiments/Experiment_90-FC1/python_results/',

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

files=[]
for key, value in paths.items():
    coloc_files =[[f'{root}/{filename}' for filename in files if 'molecule_counts.csv' in filename] for root, dirs, files in os.walk(f'{value}')]
    
    coloc_files=[item for sublist in coloc_files for item in sublist if 'all_combined' not in item]
 
    df=pd.DataFrame(coloc_files, columns=['path'])
    df[['Experiment_number', 'Treatment-replicate']]=key
    files.append(df)

collated_counts=pd.concat(files)

subunits_foci_all=[]
for (path,Treatment), row in collated_counts.groupby(['path', 'Treatment-replicate']):
    path
    path=path.replace('\\','/')
    row
    replicate_t=Treatment
    replicate=replicate_t.split('-')[-1]

    treatment_from_file=path.split('/')[-3]
    treatment=Treatment.rstrip(replicate_t[-1])
    Exp_num=list(row['Experiment_number'])[0]

    coloc=pd.read_csv(path)
    #coloc['variable1_variable2']=coloc['variable1']+'-'+coloc['variable2']
    fresh_df=pd.DataFrame(coloc['all_small_mol_count'])
    fresh_df['variable1_variable2']=coloc['variable1']+'-'+coloc['variable2']
    fresh_df['Treatment_from_file']=treatment_from_file
    fresh_df['Experiment_number']=Exp_num
    fresh_df['path']=path
    fresh_df['replicate']=replicate
    fresh_df['Treatment']=treatment
    subunits_foci_all.append(fresh_df)

subunits_foci_all=pd.concat(subunits_foci_all)



treatments_to_filter=list(subunits_foci_all['variable1_variable2'].unique())

filter_list= [
'a8-b1-110-HSPA8',
 'a8-b1-sod1-HSPA8',
 'a8-jb1-+hsp110-1h',
't1h-t1h-a8-jb1-110-0.5um',
'HSPA8-a8-647-flowin-excess',
'HSPA8-a8-b1',
'HSPA8-a8-b1+hsp110-2h',
'+HSP110-2uM-Coloc',
'+HSP110-5nM-Coloc',
'a8-b1-+sod1-2h-HSPA8',
'a8-b1-+hsp110-0.5um-HSPA8',
'a8-b1-flowout-HSPA8',
'hspa8-a8-b1-flowout',
'hspa8-a8-b1-sod1',
'hspa8-10nm',
'hspa8-darkjb1',
'hspa8-10nm-noatp',
'hspa8-pre-mixed-darkjb1',
'hspa8-647-110-100min',
'hspa8-647-110-mixed',
'hspa8-a8-atp-excess',
'a8-jb1-110-0.5um-flowout',
'a8-b1-0.5umhsp110-1h-HSPA8',
'a8-b1-2umhsp110-1h-HSPA8',
'a8-b1-HSPA8',
'a8-10nm-atp-HSPA8',
'a8-10nm-atp-1nmhsp110-HSPA8',
'a8-10nm-atp-2umhsp110-HSPA8',
 ]

subunits_foci_filter=subunits_foci_all[subunits_foci_all['variable1_variable2'].isin(filter_list)]

data_output='data/3_mols_per_foci/'

if not os.path.exists(data_output):
         os.makedirs(data_output)
subunits_foci_filter.to_csv(f'{data_output}filter_mols_per_foci.csv')

results_output=f'results/3_mols_per_foci/'
if not os.path.exists(results_output):
         os.makedirs(results_output)


# ax = sns.violinplot(x="treatment",
#  y="all_small_mol_count", 
#  data=mols_per_foci,
#   scale='width',
#    palette='Oranges')
plot_order=[
'darkA8-647A8-',
 'A8-noATP-',
 'A8-ATP-',
 'JB1-A8-',
 'JB1-A8-110-1nM-',
 'JB1-A8-110-0.5uM-',
 'JB1-A8-SOD-0.5uM-',
 'JB1-A8-110-2uM-',
 'JB1-110-0.5uM-',
 'JB1-110-2uM-',
 ]



ax=sns.violinplot(
    x="Treatment",
       y="all_small_mol_count",
        data=subunits_foci_filter,
        order=plot_order,
        palette='Oranges',
        scale='width'
)



ax = sns.boxplot(
     x="Treatment",
    y="all_small_mol_count", 
    data=subunits_foci_filter,
    order=plot_order,
    palette='Oranges')



sns.stripplot(x="Treatment",
 y="all_small_mol_count", 
  data=subunits_foci_filter,
  order=plot_order,
    palette='Oranges', 
    jitter=True,
     alpha=0.1, ax=ax)

ax.set_ylabel('# of subunits')
ax.set_xlabel('Molecule count')
#plt.ylim(0,10)
plt.xticks(rotation=90)
plt.title(f'# of HSPA8 molecules per foci')
plt.tight_layout()
plt.savefig(f'{results_output}boxplot_number_subunits_per_foci.png')
plt.show()





EXCESS_A8_OUTPUT='results/3_mols_per_foci/A8-B1-110-all/'
if not os.path.exists(EXCESS_A8_OUTPUT):
         os.makedirs(EXCESS_A8_OUTPUT)
#filter for condition A8+ATP, B1+A8 & +110
new_filter=['A8-ATP-','JB1-A8-','JB1-A8-110-1nM-','JB1-A8-110-0.5uM-','JB1-A8-110-2uM-']
plot_order=['A8-ATP-','JB1-A8-','JB1-A8-110-1nM-','JB1-A8-110-0.5uM-','JB1-A8-110-2uM-']
select=subunits_foci_filter[subunits_foci_filter['Treatment'].isin(new_filter)]
ax=sns.boxplot(x='all_small_mol_count', y='Treatment', data=select, order= plot_order, palette='BuGn_r')
plt.title(f'# HSPA8 subunits / foci')
plt.ylabel('Treatment')
plt.xlim(0,20)
plt.xlabel(f'# subunits')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f'{EXCESS_A8_OUTPUT}/foci_per_length_unit_colocalised_boxplot.png')




no_EXCESS_A8_OUTPUT='results/3_mols_per_foci/B1-110-all/'
if not os.path.exists(no_EXCESS_A8_OUTPUT):
         os.makedirs(no_EXCESS_A8_OUTPUT)
#filter for condition A8+ATP, B1+A8 & +110
new_filter=['A8-ATP-','JB1-A8-','JB1-110-0.5uM-','JB1-110-2uM-']
plot_order=['A8-ATP-','JB1-A8-','JB1-110-0.5uM-','JB1-110-2uM-']
select=subunits_foci_filter[subunits_foci_filter['Treatment'].isin(new_filter)]
ax=sns.boxplot(x='all_small_mol_count', y='Treatment', data=select, order= plot_order, palette='BuGn_r')
plt.title(f'# HSPA8 subunits / foci')
plt.ylabel('Treatment')
plt.xlim(0,20)
plt.xlabel(f'# subunits')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f'{no_EXCESS_A8_OUTPUT}/subunits_per_foci_colocalised_boxplot.png')



output='results/3_mols_per_foci/'

#plotting separately
conditions_list=['A8-noATP-']
treat='A8_noATP'
condition=subunits_foci_filter[subunits_foci_filter['Treatment'].isin(conditions_list)]

ax=sns.boxplot(x='all_small_mol_count', y='Treatment', data=subunits_foci_filter, order= conditions_list, palette='BuGn_r')
plt.title(f'# HSPA8 subunits / foci')
plt.ylabel('Treatment')
plt.xlim(0,20)
plt.xlabel(f'# subunits')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f'{output}/{treat}_foci_per_length_unit.png')
plt.show()





conditions_list=['A8-noATP-','A8-ATP-']
treat='+A8_ATP'
condition=subunits_foci_filter[subunits_foci_filter['Treatment'].isin(conditions_list)]
ax=sns.boxplot(x='all_small_mol_count', y='Treatment', data=subunits_foci_filter, order= conditions_list, palette='BuGn_r')
plt.title(f'# HSPA8 subunits / foci')
plt.ylabel('Treatment')
plt.xlim(0,20)
plt.xlabel(f'# subunits')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f'{output}/{treat}_foci_per_length_unit.png')
plt.show()



#+jb1
conditions_list=['A8-noATP-','A8-ATP-','JB1-A8-']
treat='b1_a8_atp'
condition=subunits_foci_filter[subunits_foci_filter['Treatment'].isin(conditions_list)]
ax=sns.boxplot(x='all_small_mol_count', y='Treatment', data=subunits_foci_filter, order= conditions_list, palette='BuGn_r')
plt.title(f'# HSPA8 subunits / foci')
plt.ylabel('Treatment')
plt.xlim(0,20)
plt.xlabel(f'# subunits')
plt.xticks(rotation=45)
plt.tight_layout()


plt.savefig(f'{output}/{treat}_foci_per_length_unit.png')
plt.show()


#+110 1nM
conditions_list=['A8-noATP-','A8-ATP-','JB1-A8-','JB1-A8-110-1nM-']
treat='b1_a8__110_1nm_atp'
condition=subunits_foci_filter[subunits_foci_filter['Treatment'].isin(conditions_list)]
ax=sns.boxplot(x='all_small_mol_count', y='Treatment', data=subunits_foci_filter, order= conditions_list, palette='BuGn_r')
plt.title(f'# HSPA8 subunits / foci')
plt.ylabel('Treatment')
plt.xlim(0,20)
plt.xlabel(f'# subunits')
plt.xticks(rotation=45)
plt.tight_layout()


plt.savefig(f'{output}/{treat}_foci_per_length_unit.png')
plt.show()


#+110 500nM
conditions_list=['A8-noATP-','A8-ATP-','JB1-A8-','JB1-A8-110-1nM-', 'JB1-A8-110-0.5uM-']
treat='b1_a8__110_0.5nm_atp'
condition=subunits_foci_filter[subunits_foci_filter['Treatment'].isin(conditions_list)]
ax=sns.boxplot(x='all_small_mol_count', y='Treatment', data=subunits_foci_filter, order= conditions_list, palette='BuGn_r')
plt.title(f'# HSPA8 subunits / foci')
plt.ylabel('Treatment')
plt.xlim(0,20)
plt.xlabel(f'# subunits')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f'{output}/{treat}_foci_per_length_unit.png')
plt.show()


#for EMBO conference

results_output = 'D:/presentations_conferences/2023/EMBO/figs_from_fibrils_summary_repo/'
subunits_foci_filter.to_csv(f'{results_output}filter_mols_per_foci.csv')
conditions_list=['A8-ATP-','JB1-A8-','JB1-A8-110-1nM-', 'JB1-A8-110-0.5uM-']
treat='b1_a8__110_0.5nm_atp_for_poster'
condition=subunits_foci_filter[subunits_foci_filter['Treatment'].isin(conditions_list)]
ax = sns.violinplot(x='all_small_mol_count', y="Treatment", 
                    data=condition, palette='RdYlGn', order=conditions_list, gridsize=1400, inner='quartile', linewidth=1, saturation=0.8, orient='h')
#sns.stripplot(x="Treatment", y='all_small_mol_count',
            #    data=condition, palette='RdYlGn', order=conditions_list, alpha=0.3, edgecolor='grey',ax=ax)

plt.xlim(0, 20)

plt.ylabel('HSPA8 subunits per foci')

plt.xlabel(f'Treatment')
plt.xticks(rotation=45)
#plt.tight_layout()
plt.savefig(f'{results_output}/{treat}_subunits_per_foci_poster.png')
plt.savefig(f'{results_output}/{treat}_subunits_per_foci_poster.svg')
plt.show()

#A BUNCH OF STUFF TO HELP IN TROUBLESHOOTING AND FIGURING OUT THE DATA I HAVE IE HOW MANY OF EACH TREATMENT ARE IN EACH FILE ETC.

#get medians instead, as there is less data to fiure it out from but it's calculated from thesame dataset as molecule counts
files=[]
for key, value in paths.items():
    coloc_files =[[f'{root}/{filename}' for filename in files if 'molecule_counts.csv' in filename] for root, dirs, files in os.walk(f'{value}')]
    
    coloc_files=[item for sublist in coloc_files for item in sublist ]
 
    df=pd.DataFrame(coloc_files, columns=['path'])
    df[['Experiment_number', 'Treatment-replicate']]=key
    files.append(df)

collated_info=pd.concat(files)

test=[]
for (path,Treatment, Experiment_number), df in collated_info.groupby(['path','Treatment-replicate', 'Experiment_number']):
    path
    Treatment
    df
    meds=pd.read_csv(path)
    one=len(meds['variable2'].unique())

    two=len(meds['variable2'])

    meds['treatmentnum1']=one
    meds['treatmentnum2']=two
    meds['Treatment']=Treatment
    meds['Experiment_number']=Experiment_number

    test.append(meds)
test=pd.concat(test)


multiple_treatments=[]
for (length, Experiment, Treatment), df in test.groupby(['treatmentnum1', 'Experiment_number', 'Treatment']):
     length
     Experiment
     df

     if length>1:
          l=df['variable2'].unique().tolist()
          m=df['variable1'].unique().tolist()
          t=pd.DataFrame(l, columns=['variable2'])
          t['Experiment_number']=Experiment
          t['Treatment']=Treatment
          t['variable1']=m
          t['Number_treatments']=length

          multiple_treatments.append(t)
multiple_treatments=pd.concat(multiple_treatments)



