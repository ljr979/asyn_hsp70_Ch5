"""This script finds the proportion of molecules per experiment (i.e., replicate) that have a normalised intensity above a given threshold, and pots them as a barplot with error and scatter for replicates. (Panel C, thesis Figure 5.4)

"""
from enum import unique
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math

def prop_over_thresh(df, threshold, col, group):
    """finds all data in the column of interest that is above the defined threshold in that column. Finds the percentage of 

    Args:
        df (df): dataframe with columns of interest
        threshold (float): value above which we are filtering the data
        col (str): the column to look through for values above the threshold
        group(str): what column to group on (looks through this group and finds values) above threshold for each group. i.e., replicates (experiment number looks for proportion abov threshold for each experiment)

    Returns:
        experiment_split: a dataframe with the proportion above the threshold, foreach group
    """
    Exp_split_above = []
    for treatment, dif in df.groupby('Treatment'):
        treatment
        dif
        for exp, df1 in dif.groupby(group):
            exp
            df1
            num_molecules1 = len(df1)
            print(exp)
            print(num_molecules1)
            bigs_df1 = df1[df1[col] > threshold]
            num_bigs1 = len(bigs_df1)
            proportion_bigs1 = num_bigs1/num_molecules1
            print(proportion_bigs1)
            info = [treatment, exp, proportion_bigs1]
            Exp_split_above.append(info)
            
    experiment_split = pd.DataFrame(Exp_split_above, columns=['Treatment', 'Experiment', 'Proportion'])
    experiment_split['percentage'] = experiment_split['Proportion']*100
    return experiment_split


def plot_prop_over_threshold(df, thresh, ylim, palette):
    """plot as a barplot the % of molecules that are over the defined threshold

    Args:
        df (df): df containing data
        thresh (float): threshold
        ylim (int): int to get rid of white space to 100 or show to 100%
        palette (str): colour palette
    """
    Fig,ax = plt.subplots(figsize=(5,5))
    ax = sns.barplot(data=df, x='Treatment', y='value', palette=palette, saturation=1, capsize=0.3, errcolor='0.4',linewidth=2, edgecolor='0',alpha=0.8)
    sns.scatterplot(data=df, x='Treatment', y='value', legend=False, hue='Treatment', ax=ax)
    ax.set_ylabel('HSPA8 above threshold (%)')
    ax.set_xlabel('Treatment')
    plt.ylim(0,ylim)
    plt.xticks(rotation=90)
    plt.title(f'Proportion of HSPA8 molecules per foci > {thresh}')

if __name__ == "__main__":
    
    input_folder = 'data/Figures/Figure_4/Panel_B/'
    output_folder = 'data/Figures/Figure_4/Panel_C/'

    df = pd.read_csv(f'{input_folder}1_normed_mol_per_foci.csv')
    df = df.drop([col for col in df.columns.tolist() if 'Unnamed: 0' in col], axis=1)

    thresh = 0.3
    experiment_split = prop_over_thresh(df, threshold=thresh, col='rel_intens', group='Experiment_number')

    experiment_split.to_csv(f'{output_folder}1_prop_above_threshold_{thresh}.csv')

    melted = pd.melt(experiment_split, id_vars=['Treatment', 'Proportion', 'Experiment'], value_vars=['percentage'], var_name=['stuff'])

    order = ['A8-ATP-','JB1-A8-','JB1-A8-110-0.5uM-','JB1-A8-SOD-0.5uM-']

    plot_prop_over_threshold(df=melted, thresh=thresh, ylim=50, palette='RdYlGn')



