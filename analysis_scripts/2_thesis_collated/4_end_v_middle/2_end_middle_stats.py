
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
from scipy import stats
import scikit_posthocs as sp

import statsmodels.stats.multicomp as mc

inputu=('results/4_end_middle/step_two/')
#make this the same as what it was in the previous script for plotting as this is what it's saved as and we will look for that

all_ratios_df=pd.read_csv(f'{inputu}HSPA8thesis_df_for_plotting_nozero.csv')

def oneway_anova_treatments(df, col1, col2):

    x=stats.f_oneway(df[f'{col1}'][df[f'{col2}'] == 'JB1-A8-SOD-0.5uM-'],
                df[f'{col1}'][df[f'{col2}']
                                            == 'JB1-A8-110-2uM-'],
                df[f'{col1}'][df[f'{col2}']
                                            == 'JB1-A8-110-5nM-'],
                df[f'{col1}'][df[f'{col2}']
                                            == 'JB1-A8-110-0.5uM-'],
                df[f'{col1}'][df[f'{col2}']
                                            == 'JB1-A8-'],
                df[f'{col1}'][df[f'{col2}']
                                            == 'A8-ATP-'],
                df[f'{col1}'][df[f'{col2}'] == 'A8-noATP-'])
    return x
    
    
#NOW TRY AN ANOVA? FOR THE DIF OF THE MEANS
result=all_ratios_df[all_ratios_df['Treatment']!='JB1-A8-110-1nM-']
result=result[result['Experiment_number']!='Experiment 90-2']
#dropping replicates I don't want
# result =all_ratios_df.drop(all_ratios_df[(all_ratios_df['Treatment']=='JB1-A8-SOD-0.5uM-') & (all_ratios_df['Experiment_number']=='Experiment 98-1')].index)

ends=result[result['end_or_middle']=='chaps_per_end_length']
mids=result[result['end_or_middle']=='chaps_per_middle_length']

#first, want to compare the mean of all the treatments in ends, then in middles

#this part tells us whether there is a difference in some of the means of these data sets.

x=oneway_anova_treatments(df=ends, col1='Chaperones_per_length', col2='Treatment')


#now we have a dataframe where every row is a timepoint, and compares every other timepoint to itself. each column the one being compared moves across. To get a simpler output, we can just say p_values < 0.05, and any entry that is statistically significant will be 'True'.
comp = mc.MultiComparison(ends['Chaperones_per_length'], ends['Treatment'])
post_hoc_res = comp.tukeyhsd()
summary=pd.DataFrame(post_hoc_res.summary())
summary.to_csv(f'{inputu}ends_anova_tukeyhsd_summary.csv')


x=oneway_anova_treatments(df=mids, col1='Chaperones_per_length', col2='Treatment')


#now we have a dataframe where every row is a timepoint, and compares every other timepoint to itself. each column the one being compared moves across. To get a simpler output, we can just say p_values < 0.05, and any entry that is statistically significant will be 'True'.
comp = mc.MultiComparison(mids['Chaperones_per_length'], mids['Treatment'])
post_hoc_res = comp.tukeyhsd()
summary=pd.DataFrame(post_hoc_res.summary())
summary.to_csv(f'{inputu}mids_anova_tukeyhsd_summary.csv')

