"""this script aims to take the data obtained by running the end vs middle script and plotting it (that is, looking at the binding of chaperones to the end 3 pixels vs. the middle of the fibril)

"""
import os, re
import threading
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

def wrangling_data(input_folder, radius, output_folder, zoom_factor):
    """gathers all end v middle files and finds only those fibrils that are longer than 10 pixels (long enought to actually detemine an end from a middle). Then, converts lengths to um, and counts the number of chaperones bound in each region PER fibril. 

    Args:
        input_folder (str): location of end middle files
        radius (int): size in pixels of an end region
        output_folder (str): where to save this data
        zoom_factor (float): the number of um per pixel

    Returns:
        _type_: _description_
    """
   #gather files you just saved from the previous end v middle script
    files = [filename for filename in os.listdir(input_folder) if 'end_vs_middle' in filename]
    
    data = []
    for filename in files:
        df = pd.read_csv(f'{input_folder}{filename}')
        df['end_length_pixels'] = radius
        #filter for fibrils which are greater than 10 pixels, as otherwise the end region is 'too much' of the fibril
        df = df[df['Length']>10]
        data.append(df)
        #collate all files
    data = pd.concat(data)

    middle_end_lengths_df = []
    for fibril, df in data.groupby('ID_length_conc_image'):
        #now we are looping over each fibril as a df, these were just gathered and concatinated( during the loop, each fibril has multiple pixels, each row is one of these pixels)
        fibril
        #convert the end length to um
        end_length_um = (radius*2)*zoom_factor
        #count the number of molecules that ahve been assigned 'end' on this fibril
        number_of_chaps_END = len(df[df["Where"]=='END'])
        #and middle
        number_of_chaps_MIDDLE = len(df[df["Where"]=='MIDDLE'])
        #how many um length is the end region?
        middle_length_um = (df['Length'].values.tolist()[0]*0.160)-(end_length_um)
        #now calculate the number of chaperones in each region, per um fibril length
        chaps_per_end_length = number_of_chaps_END / end_length_um
        #this if loop accounts for those molecules that didn't have any bound in the middle region
        if number_of_chaps_MIDDLE > 0:
            chaps_per_middle_length = number_of_chaps_MIDDLE/middle_length_um
        if number_of_chaps_MIDDLE == 0:
            chaps_per_middle_length = 0

        #SAVE ALL THIS INFO in new columns of the dataframe
        df['end_length_um'] = end_length_um
        df['middle_length_um'] = middle_length_um
        df['number_of_chaps_END'] = number_of_chaps_END
        df['number_of_chaps_MIDDLE'] = number_of_chaps_MIDDLE
        df['chaps_per_end_length'] = chaps_per_end_length
        df['chaps_per_middle_length'] = chaps_per_middle_length
        
        middle_end_lengths_df.append(df)
    middle_end_lengths_df = pd.concat(middle_end_lengths_df)
    #save 
    middle_end_lengths_df.to_csv(f'{output_folder}middle_end_lengths_df.csv')
    return data, middle_end_lengths_df

def plotting_chaps_per_um_end_middle(data, palette, protein, legendtitle, plot_title, order_of_experiment):
    """PLOTS this data as a barplot with overlaid datapoints

    Args:
        data (df): dataframe with the number of hcaperones bound per unit length (not including pixels less than zero)
        palette (str): colour palette to use
        protein (str): the protein we are looking at
        legendtitle (str): what to call the legend
        plot_title (str): name of the plot
        order_of_experiment (list): order to plot 
    """
    g = sns.barplot(data=data, x="Concentration", y="Chaperones_per_length", hue="end_or_middle", 
    order=order_of_experiment, 
    palette=palette, alpha=.6, edgecolor="grey")
    sns.stripplot(
    x="concentration", 
    y="Chaperones_per_length", 
    hue="end_or_middle", 
    order=order_of_experiment, 
    data=data, dodge=True, alpha=0.6, ax=g
    )
    plt.legend(title=f'{legendtitle}', labels=['End', 'Middle'], loc='upper right')
    plt.title(f'{plot_title}')
    plt.tight_layout()
    plt.show()
    
if __name__ == "__main__":

    input_folder = 'data/Analysis_workflow/2_example_python_output/5_end_vs_middle_init/step_one/'
    output_folder = 'data/Analysis_workflow/2_example_python_output/6_end_or_middle_plot/step_two/'

    if not os.path.exists(output_folder):
            os.makedirs(output_folder)

    protein = 'HSPA8'
    #this is important! needs to be adjusted according to whether the microscope was 'zoomed' or not. this is the size of each pizel in um (not nm!) we multiply each region by this to get the size in um
    zoom_factor = 0.160
    proteins = [protein for protein in os.listdir(input_folder)]
    #distance that classifies an end
    radius = 2
    data, middle_end_lengths_df = wrangling_data(input_folder=input_folder, radius=radius, output_folder=output_folder, zoom_factor=zoom_factor)

    #melt dataframe
    df_melt = pd.melt(middle_end_lengths_df, id_vars=['Unnamed: 0', 'Frame', 'Contour_ID', 'Pos.', 'X', 'Y', 'Length',
        'distance', 'Point_X', 'Point_Y', 'concentration', 'protein',
        'new_ID_hspX_hspY_num', 'EndX1, EndX2, EndY1, EndY2', 'EndX1', 'EndX2',
        'EndY1', 'EndY2', 'dist_1', 'dist_2', 'Where', 'end_length_pixels',
        'end_length_um', 'middle_length_um', 'number_of_chaps_END',
        'number_of_chaps_MIDDLE', ], value_vars=['chaps_per_end_length',
        'chaps_per_middle_length'], var_name="end_or_middle", value_name='Chaperones_per_length')


    #define hex codes for colours we want
    two_colours = ['#d20f46', '#d39ca2']
    #tell seaborn that is the palette to use
    sns.set_palette(sns.color_palette(two_colours))

    #adding in this line makes it so that we don't include the pixels where there are NO chaperones bound, so that we are comparing only those pixels with chaperones bound
    d = df_melt[df_melt['Chaperones_per_length']>0]
    #dictionary to make plot grouping names easy/attractive
    plotting_label_dict = {
      't1-t2-t3': 't1-t2-t3',
    }
    #map these new names onto the df
    d['concentration'] = d['concentration'].map(plotting_label_dict)

    #save the df that you are plotting from
    d.to_csv(f'{output_folder}{protein}_df_for_plotting.csv')
    #plot!
    plotting_chaps_per_um_end_middle(data=d, palette=two_colours, protein=protein, legendtitle='Location of molecule', plot_title=f'{protein} foci per um fibril length', order_of_experiment=['t1-t2-t3'])



