
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
from scipy import stats
import scikit_posthocs as sp
import statsmodels.stats.multicomp as mc

input_folder = ('data/Figures/Figure_4/Panel_C/')
#make this the same as what it was in the previous script for plotting as this is what it's saved as and we will look for that
all_ratios_df = pd.read_csv(f'{input_folder}1_melt-prop_above_threshold_0.3.csv')

#this part tells us whether there is a difference in some of the means of these data sets.
stats.f_oneway(all_ratios_df['value'][all_ratios_df['Treatment'] == 'JB1-A8-SOD-0.5uM-'],
               all_ratios_df['value'][all_ratios_df['Treatment']
                                           == 'JB1-A8-110-0.5uM-'],
               all_ratios_df['value'][all_ratios_df['Treatment']
                                           == 'JB1-A8-'],
               all_ratios_df['value'][all_ratios_df['Treatment']
                                           == 'A8-ATP-'],
                )

#each column the one being compared moves across. To get a simpler output, we can just say p_values < 0.05, and any entry that is statistically significant will be 'True'.
comp = mc.MultiComparison(all_ratios_df['value'], all_ratios_df['Treatment'])
post_hoc_res = comp.tukeyhsd()
summary = pd.DataFrame(post_hoc_res.summary())
summary.to_csv(f'{input_folder}anova_tukeyhsd_summary.csv')
