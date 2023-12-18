
from enum import unique
import os, re
import threading
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
import math
input_folder='results/3_mols_per_foci/normalising/'
#now i just want to try and do a little thing where I can see what proportion of molecules are 'big' in each treatment
updated_filtered_all=pd.read_csv(f'{input_folder}normalised_intensities.csv')

smalls = []
Experiments = []

threshold = 0.1
for treatment, dif in updated_filtered_all.groupby('Treatment'):
    treatment
    dif
    for exp, df1 in dif.groupby('Experiment_number'):
        exp
        df1
        num_molecules1 = len(df1)
        print(exp)
        print(num_molecules1)
        smalls_df1 = df1[df1['rel_intens'] < threshold]
        num_smalls1 = len(smalls_df1)
        proportion_smalls1 = num_smalls1/num_molecules1
        print(proportion_smalls1)
        info = [treatment, exp, proportion_smalls1]
        Experiments.append(info)


experiment_split = pd.DataFrame(
    Experiments, columns=['Treatment', 'Experiment', 'Proportion'])
experiment_split['percentage'] = experiment_split['Proportion']*100

melted=pd.melt(experiment_split, id_vars=['Treatment', 'Proportion', 'Experiment'], value_vars=['percentage'], var_name=['stuff'])

melted
melted = melted[melted['Treatment'] != 'JB1-A8-110-1nM-']
melted = melted[melted['Experiment'] != 'Experiment 90-2']
melted=melted[melted['Experiment']!='Experiment 84-2']
#dropping replicates I don't want
melted = melted.drop(melted[(
    melted['Treatment'] == 'JB1-A8-SOD-0.5uM-') & (melted['Experiment'] == 'Experiment 98-1')].index)
Fig,ax=plt.subplots()
ax=sns.barplot(data=melted, x='Treatment', y='value', palette='RdYlGn', saturation=1, capsize=0.3, errcolor='0.4',linewidth=2,edgecolor='0.5',alpha=0.5)
sns.scatterplot(data=melted, x='Treatment', y='value', legend=False, color='k', ax=ax)
# Exp_abbrev = subunits_above_thresh['Experiment'].str.split('t').str[-1].tolist()
# for i, txt in enumerate(Exp_abbrev):
#     ax.annotate(txt, (melted.Treatment[i], melted.value[i]))
ax.set_ylabel('HSPA8 above threshold (%)')
ax.set_xlabel('Treatment')
#plt.ylim(0,40)
plt.xticks(rotation=90)
plt.title(f'Proportion of HSPA8 molecules per foci < {threshold}')
plt.tight_layout()
plt.savefig(f'{input_folder}percentage_below_{threshold}.svg')
plt.savefig(f'{input_folder}percentage_below_{threshold}.pdf')
plt.savefig(f'{input_folder}percentage_below_{threshold}.png')


melted.to_csv(f'{input_folder}norm_relative_intensity_above_{threshold}_small.csv')