
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


results_output=f'results/3_mols_per_foci/normalising/'
nomalised_intensity=pd.read_csv(f'{results_output}normalised_intensities.csv')

equal=[]
for treatment, df in nomalised_intensity.groupby('Treatment'):
    random=df.sample(n=1400, random_state=1400)
    equal.append(random)
equal=pd.concat(equal)

ks_2samp(equal['rel_intens'][equal['Treatment'] == 'JB1-A8-110-0.5uM-'], equal['rel_intens'][equal['Treatment'] == 'JB1-A8-'])
#OR TRY kruskal wallis?
#visualise as a histogram?
sub = ['A8-ATP-','JB1-A8-','JB1-A8-110-5nM-','JB1-A8-110-0.5uM-','JB1-A8-110-2uM-','JB1-A8-SOD-0.5uM-']
subset=nomalised_intensity[nomalised_intensity['Treatment'].isin(sub)]
fig, ax= plt.subplots()
ax=sns.kdeplot(data=subset, x='rel_intens', hue='Treatment', fill=True, alpha=0.3, linewidth= 0.5, common_norm=False, palette='viridis')
plt.savefig(f'{results_output}kde_normed_nofilter_all_dsns.png')
plt.savefig(f'{results_output}kde_normed_nofilter_all_dsns.svg')




###### FOR RIDGELINE PLOTS
###################################################
###################################################
labels = {data_name:label for data_name, (label, data_path) in data_paths.items()}

order = sub
### sets order for histogram
font = {'weight' : 'normal', 'size'   : 12 }
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
plt.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'
####### Sets order for ridgeline plot
data_paths_ridgeline = {
    "1":'A8-ATP-',
    '2':'JB1+A8-',
    "3":'JB1+A8+110-5nM-',
    "4":'JB1+A8+110-0.5uM-',
    "5":'JB1+A8+110-2uM-',
    "6":"JB1+A8+SOD1-2uM-",
}
n_colors = len(data_paths_ridgeline)
pal = sns.color_palette(palette='BuPu', n_colors=n_colors)
plt.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'




g = sns.FacetGrid(subset, row='Treatment', hue='Treatment', aspect=10, height=3, palette='Blues')
# then we add the densities kdeplots for each condition
g.map(sns.kdeplot, 'rel_intens',
      bw_adjust=1, clip_on=False,
      fill=True, alpha=0.6, linewidth=10)
# here we add a white line that represents the contour of each kdeplot
g.map(sns.kdeplot, 'rel_intens',
      bw_adjust=1, clip_on=False,
      color="white", linewidth=6)
# here we add a horizontal line for each plot
g.map(plt.axhline, y=0,
      lw=2, clip_on=False)
# we loop over the FacetGrid figure axes (g.axes.flat) and add the condition as text with the right color
# notice how ax.lines[-1].get_color() enables you to access the last line's color in each matplotlib.Axes
for i, ax in enumerate(g.axes.flat):
    ax.text(0, .5, data_paths_ridgeline[f'{i+1}'],
            fontweight='bold', fontsize=100,
            color=ax.lines[-1].get_color())
    ax.set_facecolor((0, 0, 0, 0))  ### removes background so KDE plots are not obstructed when stacked
# we use matplotlib.Figure.subplots_adjust() function to get the subplots to overlap
g.fig.subplots_adjust(hspace=-0.7)
# eventually we remove axes titles, yticks and spines
g.set_titles("")
g.set(yticks=([]))
g.despine(bottom=True, left=True)
plt.setp(ax.get_xticklabels(), fontsize=20, fontweight='bold')
plt.xlabel('relative intensity', fontweight='bold', fontsize=100)
plt.xlim(0, 0.7)
g.savefig(f'{output_folder}/Histogram-ridgeline.svg', dpi = 600)
plt.show()