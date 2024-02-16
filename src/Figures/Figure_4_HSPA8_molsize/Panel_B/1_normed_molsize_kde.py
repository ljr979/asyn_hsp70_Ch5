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


