"""perform statistics (anova) on colocalisation data in panels B-C (DNAJB1 colocalised w fibrils, A8 w fibrils, and A8 + JB1)
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


def anova_loop(protein, df, treatments_list):
    """_summary_

    Args:
        protein (_type_): _description_
        df (_type_): _description_
        treatments_list (_type_): _description_
    """
    if len(treatments_list)==4:
        x= stats.f_oneway(
                df[df['proteins_colocalised']==protein]['percent_colocalised'][df['Treatment']==treatments_list[0]],
                df[df['proteins_colocalised']==protein]['percent_colocalised'][df['Treatment']==treatments_list[1]],
                df[df['proteins_colocalised']==protein]['percent_colocalised'][df['Treatment']==treatments_list[2]],
                df[df['proteins_colocalised']==protein]['percent_colocalised'][df['Treatment']==treatments_list[3]]
                        )
        print(x)
        

    if len(treatments_list)==3:
        y= stats.f_oneway(
            df[df['proteins_colocalised']==protein]['percent_colocalised'][df['Treatment']==treatments_list[0]],
            df[df['proteins_colocalised']==protein]['percent_colocalised'][df['Treatment']==treatments_list[1]],
            df[df['proteins_colocalised']==protein]['percent_colocalised'][df['Treatment']==treatments_list[2]]
                    )
        print(y)


if __name__ == "__main__":
    inputu=('data/Figures/Figure_5/Panel_B/')

    all_ratios_df = pd.read_csv(f'{inputu}/0_collated_colocalisation.csv')

    for protein, df in all_ratios_df.groupby('proteins_colocalised'):
        treatments_list = df['Treatment'].unique().tolist()
        print(protein)
        print(f'testing for difference between treatments for {protein}')
        anova_loop(protein, df, treatments_list)
            


