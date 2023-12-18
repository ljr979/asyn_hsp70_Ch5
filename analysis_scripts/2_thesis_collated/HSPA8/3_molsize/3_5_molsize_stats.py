
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
from scipy import stats
import scikit_posthocs as sp

import statsmodels.stats.multicomp as mc

inputu=('results/3_mols_per_foci/normalising/')
#make this the same as what it was in the previous script for plotting as this is what it's saved as and we will look for that

all_ratios_df=pd.read_csv(f'{inputu}num_subs_not_rel_intens/percent_above_0.3_for_plot_and_stats.csv')



#this part tells us whether there is a difference in some of the means of these data sets.

stats.f_oneway(all_ratios_df['value'][all_ratios_df['Treatment'] == 'JB1-A8-SOD-0.5uM-'],
               all_ratios_df['value'][all_ratios_df['Treatment']
                                           == 'JB1-A8-110-2uM-'],
               all_ratios_df['value'][all_ratios_df['Treatment']
                                           == 'JB1-A8-110-5nM-'],
               all_ratios_df['value'][all_ratios_df['Treatment']
                                           == 'JB1-A8-110-0.5uM-'],
               all_ratios_df['value'][all_ratios_df['Treatment']
                                           == 'JB1-A8-'],
               all_ratios_df['value'][all_ratios_df['Treatment']
                                           == 'A8-ATP-'],
               all_ratios_df['value'][all_ratios_df['Treatment'] == 'A8-noATP-'])


#now we have a dataframe where every row is a timepoint, and compares every other timepoint to itself. each column the one being compared moves across. To get a simpler output, we can just say p_values < 0.05, and any entry that is statistically significant will be 'True'.
comp = mc.MultiComparison(all_ratios_df['value'], all_ratios_df['Treatment'])
post_hoc_res = comp.tukeyhsd()
summary=pd.DataFrame(post_hoc_res.summary())
summary.to_csv(f'{inputu}anova_tukeyhsd_summary.csv')


#I can't test the difference between the m3 population because it doesn't have error!
#will try with the RATIOS data
gauss_pop=pd.read_csv(f'results/3_mols_per_foci/normalising/num_subs_not_rel_intens/ratios_proportions_above0.3_change.csv')
melted=pd.melt(gauss_pop, id_vars=['Treatment', 'Proportion', 'percentage', 'Unnamed: 0','Experiment'], value_vars=['ratio'], var_name=['stuff'])
stats.f_oneway(melted['value'][melted['Treatment'] == 'JB1-A8-SOD-0.5uM-'],
melted['value'][melted['Treatment']== 'JB1-A8-110-2uM-'],
melted['value'][melted['Treatment']== 'JB1-A8-110-5nM-'],
melted['value'][melted['Treatment']== 'JB1-A8-110-0.5uM-'])


#now we have a dataframe where every row is a timepoint, and compares every other timepoint to itself. each column the one being compared moves across. To get a simpler output, we can just say p_values < 0.05, and any entry that is statistically significant will be 'True'.
comp = mc.MultiComparison(melted['value'], melted['Treatment'])
post_hoc_res = comp.tukeyhsd()
summary=pd.DataFrame(post_hoc_res.summary())
summary.to_csv(f'{inputu}anova_tukeyhsd_summary.csv')
