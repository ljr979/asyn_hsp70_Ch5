
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
from scipy import stats
import scikit_posthocs as sp

import statsmodels.stats.multicomp as mc

inputu=('data/DNAJB1/1_gather_filter_foci/')
#make this the same as what it was in the previous script for plotting as this is what it's saved as and we will look for that

all_ratios_df=pd.read_csv(f'{inputu}filtered_foci_data.csv')


#split this up and turn each timepont into an array, so that we can compare between them easily
a8atp=all_ratios_df[all_ratios_df['treatment']=='JB1-alone-']
a8atp=np.array(a8atp['foci_per_pixel'])

b1a8=all_ratios_df[all_ratios_df['treatment']=='JB1-A8-']
b1a8=np.array(b1a8['foci_per_pixel'])

B1A8110mid=all_ratios_df[all_ratios_df['treatment']=='JB1-A8-110-0.5uM-']
B1A8110mid=np.array(B1A8110mid['foci_per_pixel'])



comparisons_list=['a8atp', 'b1a8',

'B1A8110mid',
]
#this part tells us whether there is a difference in some of the means of these data sets.
result = stats.kruskal(a8atp, b1a8,

B1A8110mid,
)
print(result)



#if p<0.05 for these comparisons, we need to perform Dunn's test for multiple comparisons on the same datasets
data = [
        all_ratios_df[all_ratios_df['treatment'] == 'JB1-alone-']['foci_per_pixel'], 
        all_ratios_df[all_ratios_df['treatment'] == 'JB1-A8-']['foci_per_pixel'],
       
        all_ratios_df[all_ratios_df['treatment'] == 'JB1-A8-110-0.5uM-']['foci_per_pixel'],
        
        ]
p_values = sp.posthoc_dunn(data, p_adjust='holm')

#now we have a dataframe where every row is a timepoint, and compares every other timepoint to itself. each column the one being compared moves across. To get a simpler output, we can just say p_values < 0.05, and any entry that is statistically significant will be 'True'.


signif=p_values < 0.05
signif.columns=comparisons_list

names_rows={
    1:'a8atp',
    2:'b1a8',
    3: 'B1A8110mid',

}
signif.rename(index=names_rows)
signif.to_csv(f'{inputu}stats_kruskal_dunns.csv')


#this part tells us whether there is a difference in some of the means of these data sets.

stats.f_oneway(
               all_ratios_df['foci_per_pixel'][all_ratios_df['treatment']
                                           == 'JB1-A8-110-0.5uM-'],
               all_ratios_df['foci_per_pixel'][all_ratios_df['treatment']
                                           == 'JB1-A8-'],
               all_ratios_df['foci_per_pixel'][all_ratios_df['treatment']
                                           == 'JB1-alone-'],
               )


#now we have a dataframe where every row is a timepoint, and compares every other timepoint to itself. each column the one being compared moves across. To get a simpler output, we can just say p_values < 0.05, and any entry that is statistically significant will be 'True'.
comp = mc.MultiComparison(all_ratios_df['foci_per_pixel'], all_ratios_df['treatment'])
post_hoc_res = comp.tukeyhsd()
summary=pd.DataFrame(post_hoc_res.summary())
summary.to_csv(f'{inputu}anova_tukeyhsd_summary.csv')