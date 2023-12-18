
from enum import unique
import os, re
import threading
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
import math

#output_folder='data/0_colocalisation/'
#defining a dictionary which tells us the experiment number in an easy word as the key, and the PATH to that repo as the value matching the key to call on later

results_output = f'results/2_num_foci/'


    



normalised_foci=pd.read_csv(f'{results_output}norm_per_experiment_foci_per_pixel.csv')
#I want to look at any difference (or lack thereof) in the number of fibrils vs. the number of colocalised spots on the fibril.
fibrils_count=[]
for experiment, df in normalised_foci.groupby('Experiment'):
       for treatment, df1 in df.groupby('treatment'):
        num_fibs=len(df1)
        med_foci_100=np.median(df1['foci_per_pixel'])*100
        data=[experiment, treatment, num_fibs, med_foci_100]
        fibrils_count.append(data)

        print(f'{experiment, treatment, num_fibs, med_foci_100}')

fibrils_count=pd.DataFrame(fibrils_count, columns=['Experiment', 'Treatment', 'fibril_count', 'median_foci_pixelx100'])


conditions_list=['A8-noATP-','A8-ATP-','JB1-A8-','JB1-A8-110-5nM-', 'JB1-A8-110-0.5uM-','JB1-A8-110-2uM-','JB1-A8-SOD-0.5uM-']

select=fibrils_count[fibrils_count['Treatment'].isin(conditions_list)]

for treatment, df in select.groupby('Treatment'):
    fig,ax=plt.subplots(figsize=(10,10))
    ax=sns.scatterplot(data=df, x='fibril_count', y='median_foci_pixelx100', hue='Experiment')
    plt.ylim(0,20)
    plt.title(f'{treatment}')


sns.hexbin
t=sns.light_palette("seagreen", as_cmap=True)
sns.set_palette('BuGn_r', 10)
fig, ax= plt.subplots()
ax=sns.kdeplot(select['fibril_count'], select['median_foci_pixelx100'], levels=np.arange(0, 1, 0.1),  fill=False)
plt.savefig(f'{results_output}fibril_count_vs_foci_per_pix.png')
plt.savefig(f'{results_output}fibril_count_vs_foci_per_pix.svg')

