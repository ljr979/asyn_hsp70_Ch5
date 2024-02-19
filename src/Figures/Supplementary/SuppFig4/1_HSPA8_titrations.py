
from enum import unique
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from loguru import logger
from statistics import mean as mean
from scipy.stats import sem

def plot_titration(collated_info, summary, colour, xlabel, protein):
    font = {
    'family' : 'arial',
    'weight' : 'normal',
    'size'   : 12
    }
    plt.rc('font', **font)
    plt.rcParams['svg.fonttype']= 'none'

    fig, ax = plt.subplots(figsize=(5,3))
    ax = plt.errorbar(x=summary['timepoint'], y=summary['mean'], yerr=summary['sem'], c=colour, capsize=4, linewidth=2, alpha=0.5)
    ax = plt.scatter(x=collated_info['conc_num'], y=collated_info['percent_fibs_colocalised'], c=colour, alpha=0.5, linewidths=2)
    plt.ylabel('Fibrils colocalised (%)')
    plt.xlabel (xlabel)
    plt.ylim(0,100)
    plt.title(f'Percentage of fibrils bound by {protein}')
    plt.show()

def summaries(test):
    #find the average between the replicates, and SEM
    summary_dict={}
    for timepoint, df in test.groupby('conc_num'):
        timepoint
        df
        m=mean(df['value'])
        s=sem(df['value'])
        ms=[m, s]
        summary_dict.update({timepoint:ms})
    summary_df=pd.DataFrame.from_dict(summary_dict, orient='index', columns=['mean', 'sem']).reset_index().rename(columns={'index':'timepoint'})
    return summary_df

if __name__ == "__main__":

    protein = 'HSPA8'
    input_path = f'data/Figures/Supplementary/Fig_4/'
    collated_info = pd.read_csv(f'{input_path}{protein}_titrations_collated.csv')
   
    #now I just want to plot the same thing, but I want to do this as a lineplot, so need to make the x axis numerical


    collated_info_dict={
        '0.3nM':0.3, 
        '1nM': 1, 
        '3nM': 3, 
        '10nM': 10, 
        '30nM': 30
    }


    collated_info['conc_num']=collated_info['concentration'].map(collated_info_dict)

    for_plotting = pd.melt(collated_info, id_vars=[
                        'conc_num', 'Experiment_number'], value_vars='percent_fibs_colocalised')
    summary = summaries(test=for_plotting)

    plot_titration(collated_info, summary, colour='darkviolet', xlabel=f'{protein} concentration', protein=protein)
