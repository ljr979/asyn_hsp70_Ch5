import os, re
import threading
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger


#this script aims to take the data obtained by running the end vs middle R script and plotting it (that is, looking at the binding of chaperones to the end 3 pixels vs. the middle of the fibril)

input_folder='data/DNAJB1/2_gather_filter_endmid/'
output_folder='results/DNAJB1/3_end_middle/step_two/'

    
if not os.path.exists(output_folder):
        os.makedirs(output_folder)

protein='DNAJB1'

def wrangling_data_for_plotting(input_folder, radius, output_folder, protein):
    files=[filename for filename in os.listdir(input_folder) if 'end_vs_middle.csv' in filename]
    
    # end_length=(2*radius)*0.099
    # middle_length=(fibril_length*0.099)-end_length
    data=[]
    for filename in files:
        df= pd.read_csv(f'{input_folder}{filename}')
        df['end_length_pixels']=radius
        #df['concentration']=filename
        #protein=filename.split('_')[1]
        #df['protein']=protein
        df=df[df['Length']>10]
        # mean_chaps_per_end=df['Chaps_per_end_length'].mean()
        # mean_chaps_per_middle=data['Chaps_per_middle_length'].mean()
        #summary.append(mean_chaps_per_end, mean_chaps_per_middle)
        
        data.append(df)

    data=pd.concat(data)

    middle_end_lengths_df=[]

    for fibril, df in data.groupby('ID_length_conc_image'):
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
    middle_end_lengths_df.to_csv(f'{output_folder}middle_end_lengths_df.csv')
    return data, middle_end_lengths_df



def plotting_chaps_per_pix_end_middle(data, palette, protein, legendtitle, plot_title, order_of_experiment):
    fig, ax= plt.subplots (figsize=(8, 8))
    g=sns.barplot(data=data, x="Treatment", y="Chaperones_per_length", hue="end_or_middle", 
    order=order_of_experiment, 
    palette=palette, alpha=.6, edgecolor="grey")
    #g.despine(left=True)
    #sns.stripplot(
    # x="Treatment", 
    # y="Chaperones_per_length", 
    # hue="end_or_middle", 
    # order=order_of_experiment, 
    # data=data, dodge=True, alpha=0.6, ax=g
    #)
    #g.set_axis_labels(f"Concentration of {protein} titrated onto fibrils", f"# of {protein} bound per unit fibril length (um)")
    #plt.ylim(0,7)
    plt.ylabel(f'{protein} per unit length (pixel)')
    plt.legend(title=f'{legendtitle}', labels=['End', 'Middle'], loc='upper right')
    plt.title(f'{plot_title}')
    plt.tight_layout()
    plt.ylim(0,2)
    plt.xticks(rotation=45)
    plt.savefig(f'{output_folder}thesis_End_vs_middle_{protein}_nozeros.png')
    #plt.savefig(f'{output_folder}End_vs_middle_{protein}.png')
    plt.savefig(f'{output_folder}thesis_End_vs_middle_{protein}_nozeros.svg')
    #plt.savefig(f'{output_folder}End_vs_middle_{protein}.svg')
    



#proteins=[protein for protein in os.listdir(input_folder)]
radius=2
data, middle_end_lengths_df = wrangling_data_for_plotting(input_folder=input_folder, radius=radius, output_folder=output_folder, protein=protein)



df_melt=pd.melt(middle_end_lengths_df, id_vars=['Pos.', 'X', 'Y', 'Length',
       'distance', 'Point_X', 'Point_Y', 'concentration', 'protein',
       'new_ID_hspX_hspY_num', 'EndX2',
       'EndY1', 'EndY2', 'dist_1', 'dist_2', 'Where', 'end_length_pixels',
       'end_length_pix', 'middle_length_pix', 'number_of_chaps_END',
       'number_of_chaps_MIDDLE', 'Treatment', 'Experiment_number'], value_vars=['chaps_per_end_length',
       'chaps_per_middle_length'], var_name="end_or_middle", value_name='Chaperones_per_length')



# palette="Greens"
# palette_5_col=['#ef946c', '#c4a77d', '#70877f', '#454372','#2f2963']
two_colours=['#d20f46', '#d39ca2']
palette='Blues'
sns.set_palette(sns.color_palette(palette))


legendtitle='Location of molecule'
plot_title=f'{protein} foci per pixel fibril length'
data=df_melt

#adding in this line makes it so that we don't include the pixels where there are NO chaperones bound, so that we are comparing only those pixels with chaperones bound
d=data[data['Chaperones_per_length']>0]
data=d




order_of_experiment= ['JB1-alone-', 'JB1-A8-','JB1-A8-110-0.5uM-']


data.to_csv(f'{output_folder}{protein}df_for_plotting_nozero.csv')
plotting_chaps_per_pix_end_middle(data, palette, protein, legendtitle, plot_title, order_of_experiment)


#for thesis figure
order_of_experiment= [ 'A8-ATP-','JB1-A8-','JB1-A8-110-5nM-', 'JB1-A8-110-0.5uM-','JB1-A8-110-2uM-','JB1-A8-SOD-0.5uM-']


data.to_csv(f'{output_folder}{protein}thesis_df_for_plotting_nozero.csv')
plotting_chaps_per_pix_end_middle(data, palette, protein, legendtitle, plot_title, order_of_experiment)


#everything below here does not seem to work


palette = {'A8-noATP-': '#f46d43','A8-ATP-': '#f46d43', 'JB1-A8-': '#fee999', 'JB1-A8-110-1nM-': '#e3f399', 'JB1-A8-110-5nM-': '#e3f399', 'JB1-A8-110-0.5uM-': '#9dd569','JB1-A8-110-2uM-': '#9dd569', 'JB1-A8-SOD-0.5uM-':'#d39ca2'}
#palette = 'Oranges'
protein = 'HSPA8'
legendtitle = 'Location'
plot_title = f'Location of {protein} on fibril'
change_labels = {'chaps_per_end_length': 'End',
                 'chaps_per_middle_length': 'Middle'}
data['end_or_middle'] = data['end_or_middle'].map(change_labels)
fig, ax = plt.subplots()
ax=sns.catplot(
    data=data, x="Treatment", y="Chaperones_per_length", col="end_or_middle",
    kind="bar", aspect=.6, palette=palette)
plt.xticks(rotation=90)
plt.ylabel('# HSPA8 foci per region')
plt.ylim(0, 4)

plt.title(f'{plot_title}')
#plt.tight_layout()
plt.savefig(f'{output_folder}{treat}_End_vs_middle_{protein}.svg')
plt.savefig(f'{output_folder}{treat}_End_vs_middle_{protein}.png')
