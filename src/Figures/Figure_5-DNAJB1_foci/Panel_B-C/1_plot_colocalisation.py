"""generates colocalisation plots from thesis figure 5.5 panels B-D
"""
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

def plot_coloc(df, palette, title, ylabel):
    """_summary_

    Args:
        df (_type_): _description_
        palette (_type_): _description_
        title (_type_): _description_
        ylabel (_type_): _description_
    """
    Fig, ax = plt.subplots(figsize=(5,5))
    ax = sns.barplot(x='Treatment', y='percent_colocalised', data=df, palette=palette, hue='proteins_colocalised', alpha=0.9, edgecolor='black', linewidth=3)
    sns.stripplot(x='Treatment', y='percent_colocalised', data=df, palette=palette, hue='proteins_colocalised',dodge=True, alpha=0.5, size=16, edgecolor='Black', linewidth=2, ax=ax)
    plt.legend(bbox_to_anchor=(1.1, 0), loc='upper left', ncol=1)
    ax.set_xlabel(' ')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    plt.xticks(rotation=90)
    ax.set_ylim(0,100)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    input_folder = 'data/Figures/Figure_5/Panel_B/'
    filtered = pd.read_csv(f'{input_folder}0_collated_colocalisation.csv')
    melted = pd.melt(filtered, id_vars=['Experiment_number', 'Treatment', 'proteins_colocalised'], value_vars=['proteins_colocalised','percent_colocalised'], var_name=['stuff'])

    #now we want to plot these separately. so one plot with chaps and fibrils, and one with A8 + B1
    #FIRST subset for fibril and chaps
    filt = ['JB1_fibril','HSPA8_fibril']
    fibs_chaps = filtered[filtered['proteins_colocalised'].isin(filt)]

    plot_coloc(df=fibs_chaps, palette='Reds', title='Percentage of fibrils bound by HSPA8/DNAJB1', ylabel='Fibrils bound (%)')


    #now for just A8 and b1 together
    filt = ['DNAJB1_HSPA8']
    CHAPS_chaps = filtered[filtered['proteins_colocalised'].isin(filt)]
    plot_coloc(df=CHAPS_chaps, palette='rocket_r', title='Percentage of HSPA8 bound by DNAJB1', ylabel='HSPA8 foci bound by DNAJB1 (%)')
