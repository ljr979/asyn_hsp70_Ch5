
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

def plotting (df, order, ylim=2, palette='RdYlGn'):
    fig, ax = plt.subplots(figsize=(5,5))
    plt.rcParams.update({'font.size': 12})
    ax = sns.boxplot(x='Treatment', y='log10_foci*100', data=df,
                    order=order, palette=palette)
    plt.title(f'Treatment')
    plt.ylabel('log10(# of foci/pixel)*100')
    plt.ylim(0,ylim)
    plt.xlabel(f'Treatment')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    input_folder = 'data/Figures/Figure_3/Panel_C/'

    fpp = pd.read_csv(f'{input_folder}foci_pix_all.csv')

    #log transform foci per pixel data, and multiply by 100 for plotting as many are small numbers
    fpp['log10_foci*100']=np.log10(fpp['foci_per_pixel']*100)
    fpp.to_csv(f'{input_folder}foci_per_pixel_log.csv')
    #order for plotting
    order=['A8-ATP-','JB1-A8-','JB1-A8-110-0.5uM-','JB1-A8-SOD-0.5uM-']
    plotting (df=fpp, order=order, ylim=2, palette='RdYlGn')