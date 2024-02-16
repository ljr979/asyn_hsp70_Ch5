
"""Script to normalise and plot the molecules/foci at each region (end or middle). This is corresponding to the plots in thesis Figure 5.4 Panel F


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

def calculate_percent_summary(molecule_counts, output_folder):
    """For each treatment, finds the total number of foci colocalised with a fibril (each row is a foci). Then finds the percentage of foci that are bound at the end vs. the middle. 

    Args:
        molecule_counts (_type_): _description_
        output_folder (_type_): _description_

    Returns:
        _type_: _description_
    """
    coloc_dict = []
    for treatment, df in molecule_counts.groupby('Treatment'):
        LIST = []
        total_hsp_fibs = len(df)
        total_ends = len(df[df['end_or_middle']=='END'])
        total_mids = len(df[df['end_or_middle']=='MIDDLE'])

        percent_end = total_ends/total_hsp_fibs*100

        percent_mid = total_mids/total_hsp_fibs*100

        LIST = treatment, percent_end, percent_mid
        coloc_dict.append(LIST)

    
    coloc_dict = pd.DataFrame(coloc_dict, columns=['treatment', 'percentage_end', 'percentage_middle'])
    coloc_dict.to_csv(f'{output_folder}/Summary_percentage_end_middle.csv')
    return coloc_dict

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
        w_df['rel_intens']=w_df[norm_col]/maxo
        normed_store.append(w_df)
    nomalised_intensity=pd.concat(normed_store)
    return nomalised_intensity

def plot_boxes_count(df, output_folder, dicto, order_of_experiment, save=False):
    """_summary_

    Args:
        df (_type_): _description_
        output_folder (_type_): _description_
        dicto (_type_): _description_
        order_of_experiment (_type_): _description_
        save (bool, optional): _description_. Defaults to False.
    """
    df['treatment_to_plot'] = df['Treatment'].map(dicto)
    order_of_experiment = order_of_experiment

    Fig, ax = plt.subplots(figsize=(5,5))
    ax = sns.boxplot(
        x='treatment_to_plot',
        y='rel_intens', 
        data=df, 
        hue='end_or_middle', 
        palette='Greens', 
        order=order_of_experiment,
        )
    ax.legend(loc='upper left', ncol=2)
    ax.set_ylim(-0.1,0.8)
    plt.ylabel(f'Relative HSPA8 intensity')
    plt.xlabel(' ')
    plt.title(f'relative brightness/foci @ end or middle ')
    plt.xticks(rotation=45)
    plt.tight_layout()
    if save==True:
        plt.savefig(f'{output_folder}/box_end_middle_molecules_per_foci.png')
        plt.savefig(f'{output_folder}/box_end_middle_molecules_per_foci.svg')
    plt.show()

def experiment_summary(df, col_interest, group):
    """calcualte the ratio between the molecule size at the end vs. the middle, for each experiment and treatment (using the mean)

    Args:
        df (_type_): _description_
        col_interest (_type_): _description_
        group (_type_): _description_

    Returns:
        _type_: _description_
    """
    per_exp_mean=[]
    for experiment, df in df.groupby(group):
        #calculate the ratio between middle and end (subunit/brightness relative to each other)
        for treatment, df1 in df.groupby('Treatment'):
            mean_end=np.mean(df1[df1['end_or_middle']=='END'][col_interest])
            mean_mid=np.mean(df1[df1['end_or_middle']=='MIDDLE'][col_interest])
            ratio=mean_mid/mean_end
            per_exp_mean.append([experiment, treatment, mean_end, mean_mid, ratio])
    per_exp_mean=pd.DataFrame(per_exp_mean, columns=['experiment', 'treatment', 'mean_end', 'mean_mid', 'ratio_mid_end'])
    return per_exp_mean

if __name__ == "__main__":
        
    input_folder = f'data/Figures/Figure_4/Panel_F/'

    dicto={
        'A8-ATP-':'+ATP',
        'JB1-A8-':'JB1+A8',
        'JB1-A8-110-0.5uM-':'+110 0.5uM', 
        'JB1-A8-SOD-0.5uM-':'+SOD 0.5uM'
        }
        
    order_of_experiment= [
    '+ATP',
    'JB1+A8',
    '+110 0.5uM',
    '+SOD 0.5uM']

    molecule_counts=pd.read_csv(f'{input_folder}1_all_molsize_end_middle.csv')
    #-----------------------------------------
    #now I just want to normalise the intensity as I had previously done for molsize experiment, and check if the trend is the same.
    nomalised_intensity = norm_intensity(df=molecule_counts, groupby='Experiment_number', norm_col='last_step_mol_count')

    nomalised_intensity.to_csv(f'{input_folder}2_norm_molsize_end_middle.csv')

    plot_boxes_count(nomalised_intensity, input_folder, dicto, order_of_experiment)

    coloc_dict=calculate_percent_summary(nomalised_intensity, output_folder=input_folder)

    exp_summary=experiment_summary(nomalised_intensity, col_interest='rel_intens', group='Experiment_number')



