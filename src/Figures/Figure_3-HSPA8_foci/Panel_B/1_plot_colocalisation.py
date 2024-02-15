
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

def plotting(df, order):
    fig, ax = plt.subplots(figsize=(10,10))
    ax = sns.barplot(x='Treatment', y='value', data=df,
                    order=order, palette='RdYlGn', alpha=0.7, capsize=0.3,edgecolor='black')
    sns.stripplot(x='Treatment', y='value', data=df, order=order,
                palette='RdYlGn', dodge=True, size=8, edgecolor='grey', linewidth=1, ax=ax)
    ax.set_xlabel('HSPA8 treatment')
    ax.set_ylabel(f'Fibrils bound by HSPA8 (%)')
    ax.set_title('Percentage of fibrils bound by HSPA8')
    plt.xticks(rotation=90)
    plt.ylim(0, 100)
    plt.show()

if __name__ == "__main__":
    #input folder
    input_folder='data/Figures/Figure_3/Panel_B/'
    #read in colcalisation file (chance colocalisation is subtracted already)
    coloc_collated=pd.read_csv(f'{input_folder}colocalisation.csv')
    coloc_collated.drop([col for col in coloc_collated.columns.tolist() if 'Unnamed: 0' in col], axis=1, inplace=True)
    #melt for plotting
    melted=pd.melt(coloc_collated, id_vars=['replicate','Treatment'], value_vars=['percent_colocalised'], var_name=['coloc'])

    #order for plotting
    order=['A8-ATP-', 'JB1-A8-', 'JB1-A8-110-0.5uM-', 'JB1-A8-SOD-0.5uM-']

    plotting(df=melted, order=order, )