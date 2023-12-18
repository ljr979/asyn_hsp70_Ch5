
from enum import unique
import os, re
import threading
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
import math
from scipy.stats import ttest_ind

input_folder='results/'

results_output='results/3_mols_per_foci/normalising/num_subs_not_rel_intens/'
threshold=0.3
df=pd.read_csv(f'{input_folder}3_mols_per_foci/normalising/num_subs_not_rel_intens/norm_relative_intensity_above_{threshold}.csv')
#can also do this for NOT the proportions
df=df.drop([col for col in df.columns.tolist() if 'Unnamed: 0' in col], axis=1)
df=df[df['Experiment']!='Experiment 84-2']

#df = df[df['Experiment'] != 'Experiment 90-2']
df.loc[df["Treatment"] == "JB1-110-0.5uM-", "Treatment"] = 'JB1-A8-110-0.5uM-'
#df.loc[df["Experiment"] == "Experiment 100-2", "Treatment"] = 'JB1-A8-110-0.5uM-'
#need to replace the b1 & A8 from exp 98 to say exp 100, as these were essentiall the same experiment and were always meant to be used for comparison
df.loc[(df["Experiment"] == 'Experiment 98-1') & (df["Treatment"] == 'JB1-A8-'), 'Experiment'] = 'Experiment 100-1'
df.loc[(df["Experiment"] == 'Experiment 90-2') & (df["Treatment"] == 'JB1-A8-110-1nM-'), 'Treatment'] = 'JB1-A8-110-5nM-'
#dropping sod from exp 98-1 because heath told me to
#result = df.drop(df[(df['Treatment'] == 'JB1-A8-SOD-0.5uM-') & (df['Experiment'] == 'Experiment 98-1')].index)
result = df.drop(df[(df['Treatment'] == 'JB1-A8-110-1nM-') & (df['Experiment'] == 'Experiment 81-1')].index)





x=result.sort_values('Experiment')

comparison_name='JB1-A8-'
ratios_after_A8_B1=[]
for experiment, d in x.groupby('Experiment'):
    experiment, d 
    z=d.Treatment.tolist()
    if comparison_name in z:
        initial=d['percentage'][d['Treatment']==comparison_name].tolist()[0]
        othersdf=d[d['Treatment']!= comparison_name]
        othersdf['ratio']=othersdf['percentage']/initial
        othersdf['before']=initial

        ratios_after_A8_B1.append(othersdf)

ratios=pd.concat(ratios_after_A8_B1)
ratios.to_csv(f'{results_output}ratios_proportions_above{threshold}_change.csv')


melted=pd.melt(ratios, id_vars=['Treatment', 'Proportion', 'Experiment'], value_vars=['ratio'], var_name=['stuff'])
order=['JB1-A8-110-5nM-', 'JB1-A8-110-0.5uM-', 'JB1-A8-110-2uM-','JB1-A8-SOD-0.5uM-']
fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(5,6))
ax = sns.barplot(data=melted, x='Treatment', y='value', palette='Greens', order=order,
                 capsize=0.3, errcolor='0.4', linewidth=2, edgecolor='0.5',  alpha=0.5)
#sns.scatterplot(data=melted, x='Treatment', y='value',  color='k', ax=ax)
# Exp_abbrev = melted['Experiment'].str.split('t').str[-1].tolist()
# for i, txt in enumerate(Exp_abbrev):
#     ax.annotate(txt, (melted.Treatment[i], melted.value[i]))
plt.title(f'proportion of mols >{threshold} after adding 110 or SOD')
plt.xticks(rotation=90)
plt.ylabel('prop after addition/prop w A8+B1 only (i.e. ratio increase)')
plt.tight_layout()

plt.savefig(f'{results_output}{threshold}_ratios_proportion_change.png')
plt.savefig(f'{results_output}{threshold}_ratios_proportion_change.svg')


ratios=pd.read_csv(f'{results_output}ratios_proportions_above{threshold}_change.csv')
sub = ['JB1-A8-110-0.5uM-']
subset=ratios[ratios['Treatment'].isin(sub)]
subset

melted=pd.melt(subset, id_vars=['Treatment', 'Proportion', 'Experiment'], value_vars=['percentage','before'], var_name=['after'])

#this plots the JB1 and A8 as 'before', and the +110 as 'after' rather than having SOD1 and 110 as relative increase compared to JB1 & A8, because we only have one replicate of the SOD1 FLOWED IN 
fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(5,6))
ax = sns.barplot(data=melted, x='Treatment', y='value', palette='Greens', hue='after', hue_order=['before', 'percentage'], capsize=0.3, errcolor='0.4', linewidth=2, edgecolor='0.5',  alpha=0.5)




#the below is plotting previous way of analysis (i.e. just the subset of treatments, ratio of change relative to another treatment)
fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(5,6))
ax = sns.barplot(data=melted, x='Treatment', y='value', palette='Greens', order=sub,
                 capsize=0.3, errcolor='0.4', linewidth=2, edgecolor='0.5',  alpha=0.5)
#sns.scatterplot(data=melted, x='Treatment', y='value',  color='k', ax=ax)
# Exp_abbrev = melted['Experiment'].str.split('t').str[-1].tolist()
# for i, txt in enumerate(Exp_abbrev):
#     ax.annotate(txt, (melted.Treatment[i], melted.value[i]))
plt.title(f'proportion of mols >{threshold} after adding 110 or SOD')
plt.xticks(rotation=90)
plt.ylabel('prop after addition/prop w A8+B1 only (i.e. ratio increase)')
plt.tight_layout()

plt.savefig(f'{results_output}{threshold}subset_ratios_proportion_change.png')
plt.savefig(f'{results_output}{threshold}subset_ratios_proportion_change.svg')

#------------------------------------------------
#now we want to do this for NOT proportions, but relative intensities

df = pd.read_csv(f'{input_folder}3_mols_per_foci/normalising/num_subs_not_rel_intens/normalised_intensities.csv')
#can also do this for NOT the proportions
df = df.drop([col for col in df.columns.tolist()
             if 'Unnamed: 0' in col], axis=1)
df = df[df['Treatment'] != 'JB1-A8-110-1nM-']


#df = df[df['Experiment'] != 'Experiment 90-2']
df.loc[df["Treatment"] == "JB1-110-0.5uM-", "Treatment"] = 'JB1-A8-110-0.5uM-'
#need to replace the b1 & A8 from exp 98 to say exp 100, as these were essentiall the same experiment and were always meant to be used for comparison
df.loc[(df["Experiment_number"] == 'Experiment 98-1') &
       (df["Treatment"] == 'JB1-A8-'), 'Experiment_number'] = 'Experiment 100-1'
#dropping sod from exp 98-1 because heath told me to
result = df.drop(df[(df['Treatment'] == 'JB1-A8-SOD-0.5uM-')
                 & (df['Experiment_number'] == 'Experiment 98-1')].index)


x = result.sort_values('Experiment_number')



medians=[]
for Experiment, d in x.groupby('Experiment_number'):
    for treatment, e in d.groupby('Treatment'):
        med=np.median(e['rel_intens'])
        medians.append([Experiment, treatment, med])
    

medians=pd.DataFrame(medians, columns=['Experiment', 'Treatment', 'Median_rel_intens'])


comparison_name = 'JB1-A8-'
ratios_after_A8_B1 = []
for experiment, d in medians.groupby('Experiment'):
    experiment, d
    z = d.Treatment.tolist()
    if comparison_name in z:
        initial = d['Median_rel_intens'][d['Treatment'] == comparison_name].tolist()[
            0]
        othersdf = d[d['Treatment'] != comparison_name]
        othersdf['ratio'] = othersdf['Median_rel_intens']/initial
        ratios_after_A8_B1.append(othersdf)

ratios_rel_intens= pd.concat(ratios_after_A8_B1)
ratios_rel_intens.to_csv(f'{results_output}ratios_rel_intens_increase_afterA8B1.csv')


ratios_rel_intens=pd.read_csv('results/3_mols_per_foci/normalising/num_subs_not_rel_intens/ratios_rel_intens_increase_afterA8B1.csv')
sub = ['A8-ATP-','JB1-A8-','JB1-A8-110-0.5uM-','JB1-A8-SOD-0.5uM-']
subset=ratios_rel_intens[ratios_rel_intens['Treatment'].isin(sub)]
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
ax = sns.barplot(data=subset, x='Treatment', y='ratio', palette='Reds',
                 capsize=0.3, errcolor='0.4', linewidth=2, edgecolor='0.5',  alpha=0.8)
sns.scatterplot(data=subset, x='Treatment', y='ratio', color='k', ax=ax)
plt.title(f'relative intensity change after adding 110 or SOD')
plt.xticks(rotation=90)
plt.ylabel('relative intensity after addition/rel intens b4 addition (ratio)')
plt.tight_layout()
plt.savefig(f'{results_output}subset_ratios_median_change.png')
plt.savefig(f'{results_output}subset_ratios_median_change.pdf')





#-----------just going to do the same thing but average relative intensity
test=pd.read_csv('results/3_mols_per_foci/normalising/normalised_intensities.csv')
x=[]
for experiment, df in test.groupby('Experiment_number'):
    for treatment, df1 in df.groupby('Treatment'):
        mean_rel_inten=np.mean(df1['rel_intens'])
        avg_relative_in=[treatment, experiment, mean_rel_inten]
        x.append(avg_relative_in)
x=pd.DataFrame(x, columns=['Treatment', 'Experiment', 'avg_relative_in'])
x.loc[(x["Experiment"] == 'Experiment 98-1') & (x["Treatment"] == 'JB1-A8-'), 'Experiment'] = 'Experiment 100-1'
sub = ['JB1-A8-','JB1-A8-110-0.5uM-']
subset=x[x['Treatment'].isin(sub)]
exps=['Experiment 100-1', 'Experiment 90-1']
subset=subset[subset['Experiment'].isin(exps)]

melted=pd.melt(subset, id_vars=['Treatment', 'Experiment'], value_vars=['avg_relative_in'], var_name=['after'])


fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(5,6))
ax = sns.barplot(data=melted, x='Treatment', y='value', palette='Greens', hue='Treatment',  capsize=0.3, errcolor='0.4', linewidth=2, edgecolor='0.5',  alpha=0.5)

#not significant because of the variation between experiments (p=0.7)
before=melted['value'][melted['Treatment']=='JB1-A8-']
after=melted['value'][melted['Treatment']=='JB1-A8-110-0.5uM-']
ttest_ind(before, after)

#okay now I want to see if I can account for the big difference between experiments by accounting for the number of binding sites?