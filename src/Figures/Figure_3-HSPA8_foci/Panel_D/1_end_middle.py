import os, re
import threading
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
#this script aims to take the data obtained by running the end vs middle R script and plotting it (that is, looking at the binding of chaperones to the end 3 pixels vs. the middle of the fibril)

def wrangling_data_for_plotting(input_folder, filename, radius):
    df= pd.read_csv(f'{input_folder}{filename}')
    df['end_length_pixels']=radius
    df=df[df['Length']>10]


    middle_end_lengths_df=[]

    for fibril, df in df.groupby('ID_length_conc_image'):
        fibril
        end_length_pix=(radius)

        number_of_chaps_END=len(df[df["Where"]=='END'])
        number_of_chaps_MIDDLE=len(df[df["Where"]=='MIDDLE'])
        middle_length_pix=(df['Length'].values.tolist()[0])-(end_length_pix)
        chaps_per_end_length=number_of_chaps_END / end_length_pix
        if number_of_chaps_MIDDLE > 0:
            chaps_per_middle_length=number_of_chaps_MIDDLE/middle_length_pix
        if number_of_chaps_MIDDLE == 0:
            chaps_per_middle_length=0

        df['end_length_pix']=end_length_pix
        df['middle_length_pix']=middle_length_pix
        df['number_of_chaps_END']=number_of_chaps_END
        df['number_of_chaps_MIDDLE']=number_of_chaps_MIDDLE
        df['chaps_per_end_length']=chaps_per_end_length
        df['chaps_per_middle_length']=chaps_per_middle_length
        
        middle_end_lengths_df.append(df)
    middle_end_lengths_df=pd.concat(middle_end_lengths_df)
    middle_end_lengths_df.to_csv(f'{input_folder}1_foci_per_length_pix.csv')
    return df, middle_end_lengths_df

def plotting_chaps_per_pix_end_middle(data, palette, protein, legendtitle, plot_title, order_of_experiment):
    fig, ax= plt.subplots (figsize=(8, 8))
    g=sns.barplot(data=data, x="Treatment", y="Chaperones_per_length", hue="end_or_middle", 
    order=order_of_experiment, 
    palette=palette, alpha=.6, edgecolor="grey")
    plt.ylabel(f'{protein} per unit length (pixel)')
    plt.legend(title=f'{legendtitle}', labels=['End', 'Middle'], loc='upper right')
    plt.title(f'{plot_title}')
    plt.tight_layout()
    plt.ylim(0,2)
    plt.xticks(rotation=45)
    plt.show()

def melt_df (df):
    df_melt=pd.melt(df, id_vars=['Unnamed: 0',
        'Pos.',
        'X',
        'Y',
        'Length',
        'distance',
        'Point_X',
        'Point_Y',
        'protein',
        'image_number',
        'ID_length_conc_image',
        'new_ID_hspX_hspY_num',
        'EndX1',
        'EndX2',
        'EndY1',
        'EndY2',
        'dist_1',
        'dist_2',
        'Where',
        'Experiment_number',
        'replicate',
        'Treatment',
        'end_length_pixels',
        'end_length_pix',
        'middle_length_pix',
        'number_of_chaps_END',
        'number_of_chaps_MIDDLE',], value_vars=['chaps_per_end_length',
        'chaps_per_middle_length'], var_name="end_or_middle", value_name='Chaperones_per_length')
    return df_melt


if __name__ == "__main__":
    input_folder = 'data/Figures/Figure_3/Panel_D/'

    protein='HSPA8'
    radius = 2
    data, middle_end_lengths_df = wrangling_data_for_plotting(input_folder=input_folder, filename=f'0_all_end_middle_HSPA8.csv', radius=2)

    palette = 'Blues'
    sns.set_palette(sns.color_palette(palette))
    legendtitle ='Location of molecule'
    plot_title = f'{protein} foci per pixel fibril length'
    df_melt=melt_df(df=middle_end_lengths_df)
    data=df_melt

    #adding in this line makes it so that we don't include the pixels where there are NO chaperones bound, so that we are comparing only those pixels with chaperones bound
    d=data[data['Chaperones_per_length']>0]
    d['log_chaps_per_length*100']=np.log10(d['Chaperones_per_length']*100)

    d.to_csv(f'{input_folder}2_foci_per_length_plot.csv')

    #order to plot
    order_of_experiment= [ 'A8-ATP-','JB1-A8-', 'JB1-A8-110-0.5uM-','JB1-A8-SOD-0.5uM-']

    plotting_chaps_per_pix_end_middle(d, palette, protein, legendtitle, plot_title, order_of_experiment)

