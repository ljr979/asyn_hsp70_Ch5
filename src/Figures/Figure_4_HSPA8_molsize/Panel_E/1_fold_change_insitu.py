"""This script compares the proportion of 'medium' and 'large' molecules that are present in the DNAJB1 + HSPA8 (only) treatment, with the proportion after ADDING hsp110 or SOD1. The threshold is 0.3 used here as this is the centre of the m2 (medium) gaussian distrbution. This ONLY compares directly between the same experiment. that is, only for experiments where b1+a8 are incubated together, and then 110 or sod are added. this accounts for the RELATIVE CHANGE that is DIRECTLY due to 110 or SOD1 and aimed to account for some of the experimental variability between flow cells and experiments. Corresponds to Thesis Figure 5.4, Panel E. Only have two replicates for this data, unfortunately. 

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
from scipy.stats import ttest_ind

def compare_increase(x, comparison_name):
    """_summary_

    Args:
        x (_type_): _description_
        comparison_name (_type_): _description_

    Returns:
        _type_: _description_
    """
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
    return ratios

def plot_prop_increase(df, order, palette):
    """_summary_

    Args:
        df (_type_): _description_
        order (_type_): _description_
        palette (_type_): _description_
    """
    Fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(5,6))
    plt.rcParams.update({'font.size': 12})
    ax = sns.barplot(data=df, x='Treatment', y='value', palette=palette, order=order,capsize=0.3, errcolor='0.4', linewidth=2, edgecolor='0',  alpha=0.8)
    sns.scatterplot(data=melted, x='Treatment', y='value',  color='k', ax=ax)
    plt.title(f'Prop med/large mols + 110 vs SOD')
    plt.xticks(rotation=90)
    plt.ylabel('Fold change med/large HSPA8 mols')
    plt.xlabel(' ')
    plt.tight_layout()

if __name__ == "__main__":
    
    input_folder = 'data/Figures/Figure_4/Panel_E/'
    #this file is just containing the data from panel C whereby the proportion of molecules above a certain threshold has been determined for each experiment, and then it's been sorted by experiment. Has also been filtered for those that were inthe same flow cell (i.e. ensuring that we only measure DIRECT change in HSPA8 intensity by experiment in same flow cell & same day)
    df = pd.read_csv(f'{input_folder}1_prop_over_thresh_by_exp.csv')

    #can also do this for NOT the proportions
    df = df.drop([col for col in df.columns.tolist() if 'Unnamed: 0' in col], axis=1)
    x = df.sort_values('Experiment')

    comparison_name = 'JB1-A8-'
    #loop through and only the experiments which actually directly compared the increase from JB1-A8 will be compared (from the same flow cell, where the additional treatment was directly after the 'comparison name' treatment)
    ratios = compare_increase(x, comparison_name)
    
    melted = pd.melt(ratios, id_vars=['Treatment', 'Proportion', 'Experiment'], value_vars=['ratio'], var_name=['stuff'])

    order = [ 'JB1-A8-110-0.5uM-','JB1-A8-SOD-0.5uM-']

    plot_prop_increase(melted, order, palette='YlGn')

