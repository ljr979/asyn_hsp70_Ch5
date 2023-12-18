
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
from scipy import stats
import scikit_posthocs as sp

import statsmodels.stats.multicomp as mc

inputu=('data/0_colocalisation/')
#make this the same as what it was in the previous script for plotting as this is what it's saved as and we will look for that

all_ratios_df=pd.read_csv(f'{inputu}filter_colocalisation.csv')

#NOW TRY AN ANOVA? FOR THE DIF OF THE MEANS
all_ratios_df=all_ratios_df[all_ratios_df['Treatment']!='JB1-A8-110-1nM-']
all_ratios_df=all_ratios_df[all_ratios_df['Experiment_number']!='Experiment 90-2']
#dropping replicates I don't want
result =all_ratios_df.drop(all_ratios_df[(all_ratios_df['Treatment']=='JB1-A8-SOD-0.5uM-') & (all_ratios_df['Experiment_number']=='Experiment 98-1')].index)

#this part tells us whether there is a difference in some of the means of these data sets.

stats.f_oneway(result['percent_colocalised'][result['Treatment'] == 'JB1-A8-SOD-0.5uM-'],
               result['percent_colocalised'][result['Treatment']
                                           == 'JB1-A8-110-2uM-'],
               #result['percent_colocalised'][result['Treatment']
                                          # == 'JB1-A8-110-5nM-'],
               result['percent_colocalised'][result['Treatment']
                                           == 'JB1-A8-110-0.5uM-'],
               result['percent_colocalised'][result['Treatment']
                                           == 'JB1-A8-'],
               result['percent_colocalised'][result['Treatment']
                                           == 'A8-ATP-'],
               result['percent_colocalised'][result['Treatment'] == 'A8-noATP-'])


#now we have a dataframe where every row is a timepoint, and compares every other timepoint to itself. each column the one being compared moves across. To get a simpler output, we can just say p_values < 0.05, and any entry that is statistically significant will be 'True'.
comp = mc.MultiComparison(result['percent_colocalised'], result['Treatment'])
post_hoc_res = comp.tukeyhsd()
summary=pd.DataFrame(post_hoc_res.summary())
summary.to_csv(f'{inputu}anova_tukeyhsd_summary.csv')