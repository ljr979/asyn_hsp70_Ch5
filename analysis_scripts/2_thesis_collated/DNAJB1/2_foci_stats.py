
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
from scipy import stats
import scikit_posthocs as sp

import statsmodels.stats.multicomp as mc

inputu=('results/2_num_foci/')
#make this the same as what it was in the previous script for plotting as this is what it's saved as and we will look for that

all_ratios_df=pd.read_csv(f'{inputu}filtered_normalised_by_exp.csv')


#split this up and turn each timepont into an array, so that we can compare between them easily
a8noatp=all_ratios_df[all_ratios_df['treatment']=='A8-noATP-']
a8noatp=np.array(a8noatp['normed'])

A8ATP=all_ratios_df[all_ratios_df['treatment']=='A8-ATP-']
A8ATP=np.array(A8ATP['normed'])

b1a8=all_ratios_df[all_ratios_df['treatment']=='JB1-A8-']
b1a8=np.array(b1a8['normed'])

B1A8110low=all_ratios_df[all_ratios_df['treatment']=='JB1-A8-110-5nM-']
B1A8110low=np.array(B1A8110low['normed'])

B1A8110mid=all_ratios_df[all_ratios_df['treatment']=='JB1-A8-110-0.5uM-']
B1A8110mid=np.array(B1A8110mid['normed'])

B1A8110high=all_ratios_df[all_ratios_df['treatment']=='JB1-A8-110-2uM-']
B1A8110high=np.array(B1A8110high['normed'])

B1A8SOD=all_ratios_df[all_ratios_df['treatment']=='JB1-A8-SOD-0.5uM-']
B1A8SOD=np.array(B1A8SOD['normed'])


comparisons_list=['a8noatp','A8ATP', 'b1a8',
'B1A8110low',
'B1A8110mid',
'B1A8110high',
'B1A8SOD']
#this part tells us whether there is a difference in some of the means of these data sets.
result = stats.kruskal(a8noatp,A8ATP, b1a8,
B1A8110low,
B1A8110mid,
B1A8110high,
B1A8SOD)
print(result)



#if p<0.05 for these comparisons, we need to perform Dunn's test for multiple comparisons on the same datasets
data = [all_ratios_df[all_ratios_df['treatment'] == 'A8-noATP-']['normed'], 
        all_ratios_df[all_ratios_df['treatment'] == 'A8-ATP-']['normed'], 
        all_ratios_df[all_ratios_df['treatment'] == 'JB1-A8-']['normed'],
        all_ratios_df[all_ratios_df['treatment'] == 'JB1-A8-110-5nM-']['normed'],
        all_ratios_df[all_ratios_df['treatment'] == 'JB1-A8-110-0.5uM-']['normed'],
        all_ratios_df[all_ratios_df['treatment'] == 'JB1-A8-110-2uM-']['normed'],
        all_ratios_df[all_ratios_df['treatment'] == 'JB1-A8-SOD-0.5uM-']['normed']]
p_values = sp.posthoc_dunn(data, p_adjust='holm')

#now we have a dataframe where every row is a timepoint, and compares every other timepoint to itself. each column the one being compared moves across. To get a simpler output, we can just say p_values < 0.05, and any entry that is statistically significant will be 'True'.


signif=p_values < 0.05
signif.columns=comparisons_list

names_rows={
    1:'a8noatp',
    2:'A8ATP',
    3: 'b1a8',
4:'B1A8110low',
5:'B1A8110mid',
6:'B1A8110high',
7:'B1A8SOD'
}
signif.rename(index=names_rows)
signif.to_csv(f'{inputu}stats_kruskal_dunns.csv')


#NOW TRY AN ANOVA? FOR THE DIF OF THE MEANS
all_ratios_df=all_ratios_df[all_ratios_df['treatment']!='JB1-A8-110-1nM-']
all_ratios_df=all_ratios_df[all_ratios_df['Experiment']!='Experiment 90-2']
#dropping replicates I don't want
result =all_ratios_df.drop(all_ratios_df[(all_ratios_df['treatment']=='JB1-A8-SOD-0.5uM-') & (all_ratios_df['Experiment']=='Experiment 98-1')].index)

#this part tells us whether there is a difference in some of the means of these data sets.

stats.f_oneway(result['normed'][result['treatment'] == 'JB1-A8-SOD-0.5uM-'],
               result['normed'][result['treatment']
                                           == 'JB1-A8-110-2uM-'],
               result['normed'][result['treatment']
                                           == 'JB1-A8-110-5nM-'],
               result['normed'][result['treatment']
                                           == 'JB1-A8-110-0.5uM-'],
               result['normed'][result['treatment']
                                           == 'JB1-A8-'],
               result['normed'][result['treatment']
                                           == 'A8-ATP-'],
               result['normed'][result['treatment'] == 'A8-noATP-'])


#now we have a dataframe where every row is a timepoint, and compares every other timepoint to itself. each column the one being compared moves across. To get a simpler output, we can just say p_values < 0.05, and any entry that is statistically significant will be 'True'.
comp = mc.MultiComparison(result['normed'], result['treatment'])
post_hoc_res = comp.tukeyhsd()
summary=pd.DataFrame(post_hoc_res.summary())
summary.to_csv(f'{inputu}anova_tukeyhsd_summary.csv')