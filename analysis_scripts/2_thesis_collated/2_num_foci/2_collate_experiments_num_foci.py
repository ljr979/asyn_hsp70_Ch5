
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
('Experiment 92-1', 'JB1-A8-A', 'new'): 'D:/Experiments/Experiment_92-fibrils_SOD1_control/python_results/Colocalisation/FC1/Experiment92-A8-B1-flowout/',
('Experiment 96-1', 'JB1-A8-B', 'new'): 'D:/Experiments/Experiment96-SOD1-control/python_results/Colocalisation/Experiment96_A8-B1/',
('Experiment 96-2', 'JB1-A8-C', 'new'): 'D:/Experiments/Experiment96-FC2-Hsp110-titration/python_results/Colocalisation/Experiment96-A8-B1-flowout/',
('Experiment 98-1', 'JB1-A8-D', 'new'): 'D:/Experiments/Experiment98-SOD1/python_results/Colocalisation/Experiment98-A8-B1-flowout/',


#MIX DARK A8 2UM & LABELLED A8
('Experiment 86-1','darkA8-647A8-1','new'):'D:/Experiments/Experiment86-FC1/python_results/Colocalisation/FC1/Experiment86_A8-dark1uM-A8-647-10nM/',

('Experiment 87-1','darkA8-647A8-2','new'):'D:/Experiments/Experiment87-FC1/python_results/Colocalisation/FC1/Experiment87_A8-dark1uM-A8-647-10nM/',


#JB1+A8+1nM hsp110-excessA8
('Experiment 81-1','JB1-A8-110-1nM-1', 'old'):'D:/Experiments/Experiment_81/python_results/1/',
('Experiment 81-2','JB1-A8-110-1nM-2', 'new'):'D:/Experiments/Experiment_81/python_results/2/Colocalisation/Experiment81-2_HSPA8-JB1-110_mixed/',
('Experiment 84-3','JB1-A8-110-1nM-3','new'):'D:/Experiments/Experiment_84/python_results/FC3/Colocalisation/Experiment84_A8-JB1-110-1nM-flowout/',
('Experiment 90-2','JB1-A8-110-1nM-4','new'):'D:/Experiments/Experiment_90-FC2/python_results/FC2/Colocalisation/Experiment90-A8-B1-1nMhsp110/',
#('Experiment 66-3','JB1-A8-110-1nM-5'):'D:',
('Experiment 96-2', 'JB1-A8-110-5nM-5', 'new'): 'D:/Experiments/Experiment96-FC2-Hsp110-titration/python_results/Colocalisation/Experiment96_+HSP110-5nM/',

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

output='results/2_num_foci/'
if not os.path.exists(output):
         os.makedirs(output)

#ax = sns.violinplot(x='treatment', y='foci_per_pixel', data=pix_filt, scale='width', palette='coolwarm_r')
#sns.stripplot(x='treatment', y='foci_per_pixel', data=pix_filt, palette='coolwarm_r', jitter=True,alpha=0.1, ax=ax)



ax=sns.boxplot(x='foci_per_pixel', y='treatment', data=pix_filt, palette='PuOr')
plt.title(f'foci per length unit fibril (pixel)')
plt.ylabel('Treatment')
#plt.ylim(0,0.6)
plt.xlabel(f'Foci per pixel')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f'{output}/foci_per_length_unit_colocalised_boxplot.png')




EXCESS_A8_OUTPUT='results/2_num_foci/A8-B1-110-all/'
if not os.path.exists(EXCESS_A8_OUTPUT):
         os.makedirs(EXCESS_A8_OUTPUT)
#filter for condition A8+ATP, B1+A8 & +110
new_filter=['A8-ATP-','JB1-A8-','JB1-A8-110-1nM-','JB1-A8-110-0.5uM-','JB1-A8-110-2uM-']
plot_order=['A8-ATP-','JB1-A8-','JB1-A8-110-1nM-','JB1-A8-110-0.5uM-','JB1-A8-110-2uM-']
select=pix_filt[pix_filt['treatment'].isin(new_filter)]
ax=sns.boxplot(x='foci_per_pixel', y='treatment', data=select, order= plot_order, palette='OrRd')
plt.title(f'# foci per unit fibril length (pixels)')
plt.ylabel('Treatment')
plt.xlim(0,0.6)
plt.xlabel(f'# of foci/pixel')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f'{EXCESS_A8_OUTPUT}/foci_per_length_unit_colocalised_boxplot.png')



no_EXCESS_A8_OUTPUT='results/2_num_foci/B1-110-all/'
if not os.path.exists(no_EXCESS_A8_OUTPUT):
         os.makedirs(no_EXCESS_A8_OUTPUT)
#filter for condition A8+ATP, B1+A8 & +110
new_filter=['A8-ATP-','JB1-A8-','JB1-110-0.5uM-','JB1-110-2uM-']
plot_order=['A8-ATP-','JB1-A8-','JB1-110-0.5uM-','JB1-110-2uM-']
select=pix_filt[pix_filt['treatment'].isin(new_filter)]
ax=sns.boxplot(x='foci_per_pixel', y='treatment', data=select, order= plot_order, palette='OrRd')
plt.title(f'# foci per unit fibril length (pixels)')
plt.ylabel('Treatment')
plt.xlim(0,0.6)
plt.xlabel(f'# of foci/pixel')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f'{no_EXCESS_A8_OUTPUT}/foci_per_length_unit_colocalised_boxplot.png')














#plotting separately
select=pix_filt
select=select.rename(columns={'treatment':'Treatment'})
conditions_list=['A8-noATP-']
treat='A8_noATP'
condition=select[select['Treatment'].isin(conditions_list)]

ax=sns.boxplot(x='foci_per_pixel', y='Treatment', data=select, order= conditions_list, palette='OrRd')
plt.title(f'# foci per unit fibril length (pixels)')
plt.ylabel('Treatment')
plt.xlim(0,0.6)
plt.xlabel(f'# of foci/pixel')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f'{output}/{treat}_foci_per_length_unit.png')
plt.show()

#+atp
conditions_list=['A8-noATP-','A8-ATP-']
treat='+A8_ATP'
condition=select[select['Treatment'].isin(conditions_list)]
order=[]

ax=sns.boxplot(x='foci_per_pixel', y='Treatment', data=select, order= conditions_list, palette='OrRd')
plt.title(f'# foci per unit fibril length (pixels)')
plt.ylabel('Treatment')
plt.xlim(0,0.6)
plt.xlabel(f'# of foci/pixel')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f'{output}/{treat}_foci_per_length_unit.png')
plt.show()

#+jb1
conditions_list=['A8-noATP-','A8-ATP-','JB1-A8-']
treat='b1_a8_atp'
condition=select[select['Treatment'].isin(conditions_list)]

ax=sns.boxplot(x='foci_per_pixel', y='Treatment', data=select, order= conditions_list, palette='OrRd')
plt.title(f'# foci per unit fibril length (pixels)')
plt.ylabel('Treatment')
plt.xlim(0,0.6)
plt.xlabel(f'# of foci/pixel')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f'{output}/{treat}_foci_per_length_unit.png')
plt.show()

#+110 1nM
conditions_list=['A8-noATP-','A8-ATP-','JB1-A8-','JB1-A8-110-1nM-']
treat='b1_a8__110_1nm_atp'
condition=select[select['Treatment'].isin(conditions_list)]

ax=sns.boxplot(x='foci_per_pixel', y='Treatment', data=select, order= conditions_list, palette='OrRd')
plt.title(f'# foci per unit fibril length (pixels)')
plt.ylabel('Treatment')
plt.xlim(0,0.6)
plt.xlabel(f'# of foci/pixel')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f'{output}/{treat}_foci_per_length_unit.png')
plt.show()

#+110 500nM
conditions_list=['A8-noATP-','A8-ATP-','JB1-A8-','JB1-A8-110-1nM-', 'JB1-A8-110-0.5uM-']
treat='b1_a8__110_0.5nm_atp'
condition=select[select['Treatment'].isin(conditions_list)]

ax=sns.boxplot(x='foci_per_pixel', y='Treatment', data=select, order= conditions_list, palette='OrRd')
plt.title(f'# foci per unit fibril length (pixels)')
plt.ylabel('Treatment')
plt.xlim(0,0.6)
plt.xlabel(f'# of foci/pixel')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f'{output}/{treat}_foci_per_length_unit.png')
plt.show()


#for EMBO conference

results_output = 'D:/presentations_conferences/2023/EMBO/figs_from_fibrils_summary_repo/'
pix_filt.to_csv(f'{results_output}all_foci_per_length.csv')
conditions_list=['A8-ATP-','JB1-A8-','JB1-A8-110-1nM-', 'JB1-A8-110-0.5uM-']
treat='b1_a8__110_0.5nm_atp_for_poster'
select=pix_filt[pix_filt['treatment'].isin(conditions_list)]

ax = sns.boxplot(x='foci_per_pixel', y='treatment', data=select,
                 order=conditions_list, palette='RdYlGn')
plt.title(f'# foci per unit fibril length (pixels)')
plt.ylabel('Treatment')
plt.xlim(0,0.6)
plt.xlabel(f'# of foci/pixel')
plt.xticks(rotation=45)
#plt.tight_layout()
plt.savefig(f'{results_output}/{treat}_foci_per_length_unit.png')
plt.savefig(f'{results_output}/{treat}_foci_per_length_unit.svg')
plt.show()













