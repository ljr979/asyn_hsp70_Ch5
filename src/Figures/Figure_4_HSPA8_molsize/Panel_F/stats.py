"""statistical tests for panel F, Figure 5.4
 kruskal wallis on the distributions of normalised (relative) molecules SIZES i.e. relative number of subunits. stats on the a) distribution of molecules at end v middle for each treatment (kruskal wallis and dunns, but also b)the means (anova)'
"""
from enum import unique
import os, re
import threading
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
import math
from scipy import stats
import scikit_posthocs as sp
import statsmodels.stats.multicomp as mc
# 
def kruskal_treatments(df):
    """_summary_

    Args:
        df (_type_): _description_

    Returns:
        _type_: _description_
    """
    A8ATP=df[df['Treatment']=='A8-ATP-']
    A8ATP=np.array(A8ATP['rel_size'])

    b1a8=df[df['Treatment']=='JB1-A8-']
    b1a8=np.array(b1a8['rel_size'])

    B1A8110mid=df[df['Treatment']=='JB1-A8-110-0.5uM-']
    B1A8110mid=np.array(B1A8110mid['rel_size'])

    B1A8SOD=df[df['Treatment']=='JB1-A8-SOD-0.5uM-']
    B1A8SOD=np.array(B1A8SOD['rel_size'])

    #this part tells us whether there is a difference in some of the means of these data sets.
    result = stats.kruskal(A8ATP, b1a8, B1A8110mid, B1A8SOD)
    print(result)
    return result

def dunns_treatments(df, comparisons_list):
    """_summary_

    Args:
        df (_type_): _description_
        comparisons_list (_type_): _description_

    Returns:
        _type_: _description_
    """
    data = [
        df[df['Treatment'] == 'A8-ATP-']['rel_intens'], 
        df[df['Treatment'] == 'JB1-A8-']['rel_intens'],
        df[df['Treatment'] == 'JB1-A8-110-0.5uM-']['rel_intens'],
        df[df['Treatment'] == 'JB1-A8-SOD-0.5uM-']['rel_intens']]
    p_values = sp.posthoc_dunn(data, p_adjust='holm')
    signif=p_values < 0.05
    signif.columns=comparisons_list

    names_rows={
        1:'A8ATP',
        2: 'b1a8',
        3:'B1A8110mid',
        4:'B1A8SOD'
    }
    signif.rename(index=names_rows, inplace=True)
    return signif

def kruskal_endvmid(df):
    """_summary_

    Args:
        df (_type_): _description_

    Returns:
        _type_: _description_
    """
    end=df[df['end_or_middle']=='END']
    end=np.array(end['rel_intens'])
    
    middle=df[df['end_or_middle']=='MIDDLE']
    middle=np.array(middle['rel_intens'])
    #this part tells us whether there is a difference in some of the means of these data sets.
    result = stats.kruskal(end, middle)
    print(result)
    return result

def dunns_endvmid(df):
    """_summary_

    Args:
        df (_type_): _description_

    Returns:
        _type_: _description_
    """
    data = [df[df['end_or_middle'] == 'END']['rel_intens'], 
        df[df['end_or_middle'] == 'MIDDLE']['rel_intens']]

    p_values = sp.posthoc_dunn(data, p_adjust='holm')

    signif=p_values < 0.05
        

    return signif

def anova_treatments(df):
    """_summary_

    Args:
        df (_type_): _description_

    Returns:
        _type_: _description_
    """
    result=df
    #this part tells us whether there is a difference in some of the means of these data sets.

    stats.f_oneway(result['rel_intens'][result['Treatment'] == 'JB1-A8-SOD-0.5uM-'],
                result['rel_intens'][result['Treatment']
                                            == 'JB1-A8-110-0.5uM-'],
                result['rel_intens'][result['Treatment']
                                            == 'JB1-A8-'],
                result['rel_intens'][result['Treatment']
                                            == 'A8-ATP-'],
)


    #now we have a dataframe where every row is a timepoint, and compares every other timepoint to itself. each column the one being compared moves across. To get a simpler output, we can just say p_values < 0.05, and any entry that is statistically significant will be 'True'.
    comp = mc.MultiComparison(result['rel_intens'], result['Treatment'])
    post_hoc_res = comp.tukeyhsd()
    summary=pd.DataFrame(post_hoc_res.summary())
    return summary

def anova_endmid(df, treatment):
    """_summary_

    Args:
        df (_type_): _description_
        treatment (_type_): _description_
    """
    result=df
    #this part tells us whether there is a difference in some of the means of these data sets.
    x=stats.f_oneway(result['rel_intens'][result['end_or_middle'] == 'END'],
                result['rel_intens'][result['end_or_middle']
                                            == 'MIDDLE'])
    y=pd.DataFrame(x).T.rename(columns={1:'pval'})
    y['Treatment']= treatment
   # if y['pval'].values[0] < 0.05:


    return y

if __name__ == "__main__":

   
    input_folder = 'data/Figures/Figure_4/Panel_F/'
    #read in normalised molsize end middle data
    nomalised_intensity = pd.read_csv(f'{input_folder}2_norm_molsize_end_middle.csv')

    comparisons_list=['A8ATP', 'b1a8','B1A8110mid','B1A8SOD']
    
    #subset all middle bound molecules
    mids=nomalised_intensity[nomalised_intensity['end_or_middle']=='MIDDLE']

    #this test will tell whether BETWEEN TREATMENTS, the distribution of molecule sizes at the MIDDLE of the fibril is different
    result=kruskal_treatments(df=mids)
    signif=dunns_treatments(df=mids, comparisons_list=comparisons_list)
    signif.to_csv(f'{input_folder}2_stats_kruskal_dunns_mids.csv')
    #subset end bound mols
    ends=nomalised_intensity[nomalised_intensity['end_or_middle']=='END']
    result=kruskal_treatments(df=ends)
    signif=dunns_treatments(df=ends, comparisons_list=comparisons_list)
    signif.to_csv(f'{input_folder}2_stats_kruskal_dunns_ends.csv')


    #now need to loop through each treatment and compare end and middle
    comparisons_list=['end', 'middle']
    kruskal_store=[]
    Dunns_store=[]
    for treatment, df in nomalised_intensity.groupby('Treatment'):
        treatment
        df
        result=kruskal_endvmid(df)
        t=pd.DataFrame(result).T.rename(columns={0:'stat', 1:'pval'})
        names_rows={0:f'{treatment}'}
        t.rename(index=names_rows, inplace=True)
        kruskal_store.append(t)
        if t['pval'].values[0] < 0.05:
            p_values=dunns_endvmid(df=df)
            p_values.columns=['end','middle']
            p_values['Treatment'] = treatment
            names_rows={
                1:'end',
                2:'middle',
            }
            p_values.rename(index=names_rows, inplace=True)

            Dunns_store.append(p_values)

    kruskal_end_mid=pd.concat(kruskal_store)
    dunns_end_mid=pd.concat(Dunns_store)

    kruskal_end_mid.to_csv(f'{input_folder}kruskal_wallis_end_v_middle.csv')
    dunns_end_mid.to_csv(f'{input_folder}dunns_posthoc_end_v_middle.csv')
    #----------------------------------------------------------------------------
    #now do an anova on the middle v end. do this for ends df, and mids df. then repeat for each treatment, end v mid
    summary_ends=anova_treatments(df=ends)
    summary_mids=anova_treatments(df=mids)

    anova_store=[]
    posthoc=[]
    for treatment, df in nomalised_intensity.groupby('Treatment'):
        treatment
        df
        anova=anova_endmid(df, treatment)
        anova_store.append(anova)
        if anova['pval'].values[0] <0.05:
            comp = mc.MultiComparison(df['rel_size'], df['end_or_middle'])
            post_hoc_res = comp.tukeyhsd()
            summary=pd.DataFrame(post_hoc_res.summary())
            summary['Treatment']=treatment
            posthoc.append(summary)
    anova=pd.concat(anova_store)
    tukeys=pd.concat(posthoc)

    summary_ends.to_csv(f'{input_folder}molsize_treatments_end_anova.csv')
    summary_mids.to_csv(f'{input_folder}molsize_treatments_mid_anova.csv')