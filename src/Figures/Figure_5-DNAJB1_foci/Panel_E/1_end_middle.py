import os, re
import threading
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

def wrangling_data_for_plotting(input_folder, filename, radius, output_folder):

    df = pd.read_csv(f'{input_folder}{filename}')
    df['end_length_pixels'] = radius
    df = df[df['Length']>10]

    middle_end_lengths_df=[]

    for fibril, df1 in df.groupby('ID_length_conc_image'):
        fibril
        end_length_pix=(radius)

        number_of_chaps_END=len(df1[df1["Where"]=='END'])
        number_of_chaps_MIDDLE=len(df1[df1["Where"]=='MIDDLE'])
        middle_length_pix=(df1['Length'].values.tolist()[0])-(end_length_pix)
        chaps_per_end_length=number_of_chaps_END / end_length_pix
        if number_of_chaps_MIDDLE > 0:
            chaps_per_middle_length=number_of_chaps_MIDDLE/middle_length_pix
        if number_of_chaps_MIDDLE == 0:
            chaps_per_middle_length=0


        df1['end_length_pix']=end_length_pix
        df1['middle_length_pix']=middle_length_pix
        df1['number_of_chaps_END']=number_of_chaps_END
        df1['number_of_chaps_MIDDLE']=number_of_chaps_MIDDLE
        df1['chaps_per_end_length']=chaps_per_end_length
        df1['chaps_per_middle_length']=chaps_per_middle_length
        
        middle_end_lengths_df.append(df1)
    middle_end_lengths_df=pd.concat(middle_end_lengths_df)
    middle_end_lengths_df.to_csv(f'{output_folder}1_foci_per_length_.csv')
    return middle_end_lengths_df

def plotting_chaps_per_pix_end_middle(data, palette, protein, yval, legendtitle, plot_title, order_of_experiment):
    fig, ax= plt.subplots (figsize=(5, 5))
    g=sns.barplot(data=data, x="Treatment", y=yval, hue="end_or_middle", 
    order=order_of_experiment, 
    palette=palette, alpha=.9, edgecolor="black", linewidth=3)
    plt.ylabel(f'{protein} per unit length (pixel)')
    plt.legend(title=f'{legendtitle}', labels=['End', 'Middle'], loc='upper right')
    plt.title(f'{plot_title}')
    plt.tight_layout()
    plt.ylim(0,2.5)
    plt.xticks(rotation=45)

def melt(df):
    df_melt=pd.melt(df, id_vars=['Pos.', 'X', 'Y', 'Length',
       'distance', 'Point_X', 'Point_Y', 'concentration', 'protein',
       'new_ID_hspX_hspY_num', 'EndX2',
       'EndY1', 'EndY2', 'dist_1', 'dist_2', 'Where', 'end_length_pixels',
       'end_length_pix', 'middle_length_pix', 'number_of_chaps_END',
       'number_of_chaps_MIDDLE', 'Treatment', 'Experiment_number'], value_vars=['chaps_per_end_length',
       'chaps_per_middle_length'], var_name="end_or_middle", value_name='Chaperones_per_length')
    return df_melt

if __name__ == "__main__":
    input_folder='data/Figures/Figure_5/Panel_E/'

    protein='DNAJB1'

    radius=2

    middle_end_lengths_df = wrangling_data_for_plotting(input_folder=input_folder, filename='0_filtered_end_vs_middle.csv', radius=radius, output_folder=input_folder)

    df_melt=melt(df=middle_end_lengths_df)
    #adding in this line makes it so that we don't include the pixels where there are NO chaperones bound, so that we are comparing only those pixels with chaperones bound
    d=df_melt[df_melt['Chaperones_per_length']>0]
    #now I want to make this more consistent with the foci per pixel figure, and the HSPA8 figure. so I am going to log10 the chaps / length *100 (i.e. per nm) for plotting (now it is density / nm per region)
    d['log_chaps_per_length*100']=np.log10(d['Chaperones_per_length']*100)
    #save!
    d.to_csv(f'{input_folder}2_end_midle_log100_plotting.csv')


    #-----------plotting section
    order_of_experiment= ['JB1-alone-', 'JB1-A8-','JB1-A8-110-0.5uM-', 'JB1-A8-SOD-0.5uM-']

    plotting_chaps_per_pix_end_middle(data=d, palette='rocket_r', protein='DNAJB1', yval='log_chaps_per_length*100', legendtitle='Location of molecule', plot_title=f'{protein} foci per pixel fibril length', order_of_experiment=order_of_experiment)

