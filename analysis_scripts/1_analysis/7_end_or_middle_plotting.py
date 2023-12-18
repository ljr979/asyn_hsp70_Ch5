import os, re
import threading
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger


#this script aims to take the data obtained by running the end vs middle R script and plotting it (that is, looking at the binding of chaperones to the end 3 pixels vs. the middle of the fibril)

input_folder='python_results/end_vs_middle/step_one/'
output_folder='python_results/end_vs_middle/step_two/'

    
if not os.path.exists(output_folder):
        os.makedirs(output_folder)

protein='HSPA8'

def wrangling_data_for_plotting(input_folder, radius, output_folder, protein):
    files=[filename for filename in os.listdir(input_folder) if 'end_vs_middle' in filename]
    
    # end_length=(2*radius)*0.099
    # middle_length=(fibril_length*0.099)-end_length
    data=[]
    for filename in files:
        df= pd.read_csv(f'{input_folder}{filename}')
        df['end_length_pixels']=radius
        #df['concentration']=filename
        #protein=filename.split('_')[1]
        #df['protein']=protein
        df=df[df['Length']>8]
        # mean_chaps_per_end=df['Chaps_per_end_length'].mean()
        # mean_chaps_per_middle=data['Chaps_per_middle_length'].mean()
        #summary.append(mean_chaps_per_end, mean_chaps_per_middle)
        
        data.append(df)

    data=pd.concat(data)

    middle_end_lengths_df=[]

    for fibril, df in data.groupby('ID_length_conc_image'):
        fibril
        end_length_um=(radius*2)*0.160

        number_of_chaps_END=len(df[df["Where"]=='END'])
        number_of_chaps_MIDDLE=len(df[df["Where"]=='MIDDLE'])
        middle_length_um=(df['Length'].values.tolist()[0]*0.160)-(end_length_um)
        chaps_per_end_length=number_of_chaps_END / end_length_um
        if number_of_chaps_MIDDLE > 0:
            chaps_per_middle_length=number_of_chaps_MIDDLE/middle_length_um
        if number_of_chaps_MIDDLE == 0:
            chaps_per_middle_length=0


        df['end_length_um']=end_length_um
        df['middle_length_um']=middle_length_um
        df['number_of_chaps_END']=number_of_chaps_END
        df['number_of_chaps_MIDDLE']=number_of_chaps_MIDDLE
        df['chaps_per_end_length']=chaps_per_end_length
        df['chaps_per_middle_length']=chaps_per_middle_length
        
        middle_end_lengths_df.append(df)
    middle_end_lengths_df=pd.concat(middle_end_lengths_df)
    middle_end_lengths_df.to_csv(f'{output_folder}middle_end_lengths_df.csv')
    return data, middle_end_lengths_df



def plotting_chaps_per_um_end_middle(data, palette, protein, legendtitle, plot_title, order_of_experiment):
    
    g=sns.barplot(data=data, x="concentration", y="Chaperones_per_length", hue="end_or_middle", 
    order=order_of_experiment, 
    palette=palette, alpha=.6, edgecolor="grey")
    #g.despine(left=True)
    sns.stripplot(
    x="concentration", 
    y="Chaperones_per_length", 
    hue="end_or_middle", 
    order=order_of_experiment, 
    data=data, dodge=True, alpha=0.6, ax=g
    )
    #g.set_axis_labels(f"Concentration of {protein} titrated onto fibrils", f"# of {protein} bound per unit fibril length (um)")
    #plt.ylim(0,7)
    plt.legend(title=f'{legendtitle}', labels=['End', 'Middle'], loc='upper right')
    plt.title(f'{plot_title}')
    plt.tight_layout()
    plt.savefig(f'{output_folder}End_vs_middle_{protein}_nozeros.png')
    #plt.savefig(f'{output_folder}End_vs_middle_{protein}.png')
    plt.savefig(f'{output_folder}End_vs_middle_{protein}_nozeros.svg')
    #plt.savefig(f'{output_folder}End_vs_middle_{protein}.svg')
    



proteins=[protein for protein in os.listdir(input_folder)]
radius=2
data, middle_end_lengths_df = wrangling_data_for_plotting(input_folder=input_folder, radius=radius, output_folder=output_folder, protein=protein)



df_melt=pd.melt(middle_end_lengths_df, id_vars=['Unnamed: 0', 'Frame', 'Contour_ID', 'Pos.', 'X', 'Y', 'Length',
       'distance', 'Point_X', 'Point_Y', 'concentration', 'protein',
       'new_ID_hspX_hspY_num', 'EndX1, EndX2, EndY1, EndY2', 'EndX1', 'EndX2',
       'EndY1', 'EndY2', 'dist_1', 'dist_2', 'Where', 'end_length_pixels',
       'end_length_um', 'middle_length_um', 'number_of_chaps_END',
       'number_of_chaps_MIDDLE', ], value_vars=['chaps_per_end_length',
       'chaps_per_middle_length'], var_name="end_or_middle", value_name='Chaperones_per_length')



# palette="Greens"
# palette_5_col=['#ef946c', '#c4a77d', '#70877f', '#454372','#2f2963']
two_colours=['#d20f46', '#d39ca2']
sns.set_palette(sns.color_palette(two_colours))
palette=two_colours

legendtitle='Location of molecule'
plot_title=f'{protein} foci per um fibril length'
data=df_melt

#adding in this line makes it so that we don't include the pixels where there are NO chaperones bound, so that we are comparing only those pixels with chaperones bound
d=data[data['Chaperones_per_length']>0]
data=d


order_of_experiment = ['+hsp110-0.5uM', '+SOD1-0.5uM']


plotting_label_dict = {'A8-B1-110': '+hsp110-0.5uM',
                       'A8-B1-SOD1': '+SOD1-0.5uM'

}
data['concentration']=data['concentration'].map(plotting_label_dict)
data.to_csv(f'{output_folder}{protein}df_for_plotting_nozero.csv')
plotting_chaps_per_um_end_middle(data, palette, protein, legendtitle, plot_title, order_of_experiment)



