"""This script normalises the initial intensity of each HSPA8 molecule colocalised with a fibril, to the max in that treatment. Then it plots this normalised fluorescence intensity as a KDE histogram plot. This corresponds to thesis Figure 5.4, Panel B
"""
from enum import unique
import os, re
import threading
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
import math
from scipy.stats import ks_2samp


def norm_intensity(df, groupby, norm_col):
    """now I want to normalise on each day, just the relative intensity not the number of subunits. so need to loop over for each experiment, and normalise that fluorescence intensity to between 0 and 1 (divide by the max)

    Args:
        df (df): dataframe with molecule sizes from py4bleaching
        groupby (str): the df column that you want to normalise by (for thesis, this is grouping by the day of the imaging i.e. Experiment number)
        norm_col (str): the column to normalise (intensity or molsize)

    Returns:
        df: the dataframe with normalised column
    """
    #
    normed_store=[]
    for experiment, w_df in df.groupby(groupby):
        experiment
        w_df
        maxo=w_df[norm_col].max()
        w_df['rel_intens']=w_df['norm_col']/maxo
        normed_store.append(w_df)
    nomalised_intensity=pd.concat(normed_store)
    return nomalised_intensity

def plot_kde(df, x, hue , palette):
    fig, ax = plt.subplots(figsize=(7,5))
    plt.rcParams.update({'font.size': 12})
    ax=sns.kdeplot(data=df, x=x, hue=hue, fill=True, alpha=0.3, linewidth= 0.5, common_norm=False, palette=palette)
    plt.xlabel('Normalised HSPA8 intensity')

if __name__ == "__main__":

    input_folder='data/Figures/Figure_4/Panel_B/'

    df = pd.read_csv(f'{input_folder}0_mol_per_foci.csv')
    df = df.drop([col for col in df.columns.tolist() if 'Unnamed: 0' in col], axis=1)
    nomalised_intensity=norm_intensity(df=df, groupby='Experiment_number', norm_col='Intensity')
    nomalised_intensity.to_csv(f'{input_folder}1_normed_mol_per_foci.csv')

    plot_kde(df=nomalised_intensity, x='rel_intens', hue='Treatment',palette='viridis')





    #okay so now want to find, as I have in previous script, how many molecules relative intensity is above a threshold. WIll set this based on the A8-ATP alone. so , above 0.2. 
bigs = []
Experiments = []

threshold = 0.3
for treatment, dif in nomalised_intensity.groupby('Treatment'):
    treatment
    dif
    for exp, df1 in dif.groupby('Experiment_number'):
        exp
        df1
        num_molecules1 = len(df1)
        print(exp)
        print(num_molecules1)
        bigs_df1 = df1[df1['rel_intens'] > threshold]
        num_bigs1 = len(bigs_df1)
        proportion_bigs1 = num_bigs1/num_molecules1
        print(proportion_bigs1)
        info = [treatment, exp, proportion_bigs1]
        Experiments.append(info)


experiment_split = pd.DataFrame(
    Experiments, columns=['Treatment', 'Experiment', 'Proportion'])
experiment_split['percentage'] = experiment_split['Proportion']*100

experiment_split.to_csv(f'{results_output}norm_relative_intensity_above_{threshold}.csv')

# now plot as we did in the subunits 'raw' script
results_output='results/3_mols_per_foci/normalising/num_subs_not_rel_intens/'
threshold=0.3
experiment_split=pd.read_csv(f'{results_output}norm_relative_intensity_above_{threshold}.csv')
melted=pd.melt(experiment_split, id_vars=['Treatment', 'Proportion', 'Experiment'], value_vars=['percentage'], var_name=['stuff'])

melted
#dropping replicates I don't want
#melted = melted[melted['Treatment'] != 'JB1-A8-110-1nM-']
melted = melted[melted['Experiment'] != 'Experiment 84-2']
melted.loc[(melted["Experiment"] == 'Experiment 90-2') & (melted["Treatment"] == 'JB1-A8-110-1nM-'), 'Treatment'] = 'JB1-A8-110-5nM-'
melted = melted.drop(melted[(melted['Treatment'] == 'JB1-A8-110-1nM-') & (melted['Experiment'] == 'Experiment 81-1')].index)
melted = melted.drop(melted[(melted['Treatment'] == 'JB1-A8-SOD-0.5uM-') & (melted['Experiment'] == 'Experiment 98-1')].index)
sub = ['A8-ATP-','JB1-A8-','JB1-A8-110-0.5uM-','JB1-A8-SOD-0.5uM-']

subset=melted[melted['Treatment'].isin(sub)]

Fig,ax=plt.subplots(figsize=(10,10))
ax=sns.barplot(data=subset, x='Treatment', y='value', palette='RdYlGn', saturation=1, capsize=0.3, errcolor='0.4',linewidth=2,edgecolor='0.5',alpha=0.5)
sns.scatterplot(data=subset, x='Treatment', y='value', legend=False, color='k', ax=ax)
# Exp_abbrev = melted['Experiment'].str.split('t').str[-1].tolist()
# for i, txt in enumerate(Exp_abbrev):
#     ax.annotate(txt, (melted.Treatment[i], melted.value[i]))
ax.set_ylabel('HSPA8 above threshold (%)')
ax.set_xlabel('Treatment')
plt.ylim(0,100)
plt.xticks(rotation=90)
plt.title(f'Proportion of HSPA8 molecules per foci > {threshold}')

plt.savefig(f'{results_output}subset_percentage_above_{threshold}-.svg')
plt.savefig(f'{results_output}subset_percentage_above_{threshold}-.pdf')
plt.savefig(f'{results_output}subset_percentage_above_{threshold}-.png')

subset.to_csv(f'{results_output}subset_percent_above_{threshold}_for_plot_and_stats.csv')



#now I want to find the fibril density in each experiment, for each treatment
foci_fibril_data=pd.read_csv(f'{input_folder}2_num_foci/all_foci_per_length.csv')
foci_fibril_data=foci_fibril_data.drop([col for col in foci_fibril_data.columns.tolist() if 'Unnamed: 0' in col], axis=1)
#in the foci data it only has exp 86-1 which isn't actually true. labelling this for this so it maps correctly, to 86-1


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


fibrils_together = pd.DataFrame(fibrils, columns=['Treatment', 'Experiment', 'num_fibs'])
fibrils_together.loc[fibrils_together["Treatment"] =="JB1-110-0.5uM-", "Treatment"] = 'JB1-A8-110-0.5uM-'
fibrils_together.loc[fibrils_together["Treatment"] == "JB1-110-2uM-", "Treatment"] = 'JB1-A8-110-2uM-'
#now I am making a temporary column which is the two columns smooshed together so I can map the number of fibrils to these things separately
fibrils_together['temp']=fibrils_together['Treatment'] + fibrils_together['Experiment']
experiment_split['temp'] = experiment_split['Treatment'] + experiment_split['Experiment']


d=dict(zip(fibrils_together['temp'], fibrils_together['num_fibs']))
experiment_split['num_fibs'] = experiment_split['temp'].map(d)

#now drop the temporary column
fibrils_together=fibrils_together.drop([col for col in fibrils_together.columns.tolist() if 'temp' in col], axis=1)
experiment_split=experiment_split.drop([col for col in experiment_split.columns.tolist() if 'temp' in col], axis=1)
experiment_split.to_csv(f'{results_output}normalised_intensity_num_fibrils.csv')

#plot this to see any correlation?
#plotted separately
for treatment, df in experiment_split.groupby('Treatment'):
    Fig, ax = plt.subplots()
    ax = sns.scatterplot(x='num_fibs', y='percentage', data=df,
                         hue='Treatment', legend='brief', palette='magma')
    plt.ylim(0, 100)

Fig, ax = plt.subplots()
ax = sns.scatterplot(x='num_fibs', y='percentage', data=experiment_split,
                     hue='Experiment', legend='brief', palette='viridis')
plt.ylim(0, 20)
plt.legend(bbox_to_anchor=(1.25, 1), borderaxespad=0)


#plotting on the same graph
Fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 6))
ax1.set_title('by_experiment')
sns.scatterplot(x='num_fibs', y='percentage', data=experiment_split,
                hue='Experiment', legend='brief', palette='Blues', ax=ax1)
ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.0),
           ncol=4, fancybox=True, shadow=True)
ax1.set_ylim(0, 100)

ax2.set_title('by_treatment')
sns.scatterplot(x='num_fibs', y='percentage', data=experiment_split,
                hue='Treatment', legend='brief', palette='Greens', ax=ax2)
ax2.legend(loc='upper center', bbox_to_anchor=(0.5, 1.0),
           ncol=2, fancybox=True, shadow=True)
ax2.set_ylim(0, 100)
plt.savefig('results/3_mols_per_foci/proportion_NORMALISED_A8_big8_vs_fibs.png')





