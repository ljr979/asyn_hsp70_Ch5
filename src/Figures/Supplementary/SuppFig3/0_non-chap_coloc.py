"""_summary_
"""
from enum import unique
import os
import re
import threading
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
import math

def plot_coloc(df, palette, title, ylabel, hue, ylim):
    """_summary_

    Args:
        df (_type_): _description_
        palette (_type_): _description_
        title (_type_): _description_
        ylabel (_type_): _description_
    """
    Fig, ax = plt.subplots(figsize=(5,5))
    ax = sns.barplot(x='Treatment', y='percent_colocalised', data=df, palette=palette, hue=hue, alpha=0.9, edgecolor='black', linewidth=3)
    sns.stripplot(x='Treatment', y='percent_colocalised', data=df, palette=palette, hue=hue,dodge=True, alpha=0.5, size=16, edgecolor='Black', linewidth=2, ax=ax)
    plt.legend(bbox_to_anchor=(1.1, 0), loc='upper left', ncol=1)
    ax.set_xlabel(' ')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    plt.xticks(rotation=90)
    ax.set_ylim(0,ylim)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
        
    input_folder='data/Figures/Supplementary/Fig_3/'
    
    files=[x for x in os.listdir(input_folder) if 'colocalisation' in x]

    f=[]
    for n in files:
        x=pd.read_csv(f'{input_folder}{n}')
        x['Treatment']=n.split('_')[0]
        f.append(x)

    f=pd.concat(f)

    f.rename(columns={'percent_colocalisation': 'percent_colocalised'}, inplace=True)
    
    plot_coloc(df=f, palette='BuPu', title='non-specific colocalisation', ylabel='fibrils colocalised (%)', hue='Treatment', ylim=20)