"""_summary_

"""
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
from scipy import stats
import scikit_posthocs as sp
import statsmodels.stats.multicomp as mc

def summaries(test): 
    """Find the average between the replicates, and SEM

    Args:
        test (df): dataframe with the colocalisation data

    Returns:
        df: dataframe with the summaries in it
    """
    
    agg_func_math = {
    'log10_foci*100': [ 'mean',  'sem']
    }
    summary_df = test.groupby(['Experiment'], as_index=False).agg(agg_func_math).round(2).reset_index()
    summary_df.columns = ['_'.join(col).rstrip('_') for col in summary_df.columns.values]
    return summary_df

def kruskal_treatments(df, value_col,tlist=[]):
    """_summary_

    Args:
        df (_type_): _description_

    Returns:
        _type_: _description_
    """
    one=df[df['Treatment']==tlist[0]]
    one=np.array(one[value_col])

    two=df[df['Treatment']==tlist[1]]
    two=np.array(two[value_col])

    three=df[df['Treatment']==tlist[2]]
    three=np.array(three[value_col])

    four=df[df['Treatment']==tlist[3]]
    four=np.array(four[value_col])

    #this part tells us whether there is a difference in some of the means of these data sets.
    result = stats.kruskal(one, two, three, four)
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

def anova(value_col, df, treatments_list):
    """_summary_

    Args:
        protein (_type_): _description_
        df (_type_): _description_
        treatments_list (_type_): _description_
    """
    if len(treatments_list)==4:
        x= stats.f_oneway(
                df[value_col][df['Treatment']==treatments_list[0]],
                df[value_col][df['Treatment']==treatments_list[1]],
                df[value_col][df['Treatment']==treatments_list[2]],
                df[value_col][df['Treatment']==treatments_list[3]]
                        )

        print(x)

if __name__ == "__main__":

    input_folder = 'data/Figures/Figure_5/Panel_D/'

    all_ratios_df = pd.read_csv(f'{input_folder}1_log_100_foci_pix.csv')
    all_ratios_df = all_ratios_df[[
                                'hsp_count',
                                'lengths (pixel)',
                                'protein',
                                'foci_per_pixel',
                                'replicate',
                                'Experiment',
                                'Treatment',
                                'log10_foci*100']]


    treatments_list = all_ratios_df['Treatment'].unique().tolist()

    anova(value_col='log10_foci*100', df=all_ratios_df, treatments_list=treatments_list)
    comp = mc.MultiComparison(all_ratios_df['log10_foci*100'], all_ratios_df['Treatment'])
    post_hoc_res = comp.tukeyhsd()
    summary = pd.DataFrame(post_hoc_res.summary())
    summary.to_csv(f'{input_folder}anova_tukeyhsd_summary.csv')
