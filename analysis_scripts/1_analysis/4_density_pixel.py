import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

#Script to find the density of fluorescently labelled molecules bound per individual fibril (as molecule / pixel)


input_top='imagejresults/Colocalisation/'
output_top='python_results/Colocalisation/'

#find sub-treatments within colocalisation folder
treatments= [folder for folder in os.listdir(input_top)]


for exp in treatments: 

    input_folder=f'{input_top}{exp}/'
    #make new folder to save this data
    output_folder=f'{output_top}{exp}/3_foci_per_pixel/'
    #define which experiment to analyse (i.e., hspa8 and fibrils only , or multiple colocalisations)
    experiments=['fibrils-647']

    alls=[]
    #loop over the experiment folders
    for c in experiments:
        collated=[]
        folder = f'{input_folder}{c}/'
        #gather data on treatment
        concentration = folder.split('/')[-3].split('_')[-1]
        #gather data of experiment
        experiment_number= folder.split('/')[-3].split('_')[0]
        #make output folder
        output = f'{output_folder}{c}/'
        
        if not os.path.exists(output):
                os.makedirs(output)
        #gather files
        colocalised_files = [filename for filename in os.listdir(f'{folder}') if '.csv' in filename] 
        
        colocalised_fibrils = []
        for filepath in colocalised_files: 
            #loop over files
            data = pd.read_csv(f'{folder}{filepath}')
            #data = data.drop(['Class'], axis=1)
            protein=filepath.split('_')[0]
            #make new columns with identifiers about the experiment/treatment
            data['concentration']=concentration
            data['protein']=protein
            data['image_number']=filepath
            colocalised_fibrils.append(data)
        #collate all of them
        colocalised_fibrils=pd.concat(colocalised_fibrils)
    
        collated.append(colocalised_fibrils)
        collated=pd.concat(collated)

        collated.rename(columns={'Contour ID':'Contour_ID'}, inplace=True)

        #now we have a big dataframe with only the colocalisation data and we can get the number of foci on each fibril from this
        collated2=[]
        #iterate over each 'fibril'. Every fibril has a row assigned to each individual pixel of that fibril, so we are going over each pixel, for each fibril and giving it a unique name with the length of that fibril, the image it came from, the treatment and the 'contour id' assigned at imagej
        for Contour_ID, row in collated.iterrows():
            Contour_ID
            row
            row['ID_length_conc_image']=str(row['Contour_ID'])+'_'+str(row['Length'])+'_'+str(row['concentration'])+'_'+str(row['image_number'])
            row=pd.DataFrame(row).T
            collated2.append(row)

        collated2=pd.concat(collated2)
        #then we filter this for only those rows that are colocalised with a molecule
        coloc_collated2=collated2[collated2['distance']!=-1]
        #save this file to work on from here
        coloc_collated2.to_csv(f'{output_folder}{c}/fibril_collated-colocal-only_data.csv')
        #save the non-filtered for colocalisation version too.
        collated2.to_csv(f'{output_folder}{c}/HSPA8_fibril_collated-colocal_data.csv')
        
        #define two dictionaries
        proteins = {}
        lengths = {}


        #then make a new dataframe where I count the number of colocalised SPOTS on each contour, and then add in the matching fibril length to this counted number of foci
        #counts the number of individual fibril pixels
        counts = pd.DataFrame(coloc_collated2['ID_length_conc_image'].value_counts()).reset_index().rename(columns={'index':'ID_length_conc_image', 'ID_length_conc_image':'hsp_count'})

        #now get all of the unique names
        contourIDs = [name for name in counts['ID_length_conc_image']]

        #now get the lengths of each of these
        for contourID in contourIDs:
            length = coloc_collated2.loc[coloc_collated2['ID_length_conc_image'] == contourID, 'Length'].iloc[0]
            #assign the length of the fibril to the dictionary we made earlier, and the key is the contour ID name
            lengths.update({contourID:length})

        #now get the protein that matches that colocalised spot
        for contourID in contourIDs:
            protein = coloc_collated2.loc[coloc_collated2['ID_length_conc_image'] == contourID, 'protein'].iloc[0]
            #assign this to the contour ID in the dictionary
            proteins.update({contourID:protein})



        #now make a new dictionary where the keys are all values in the column with the unique, colocalised fibril names, and the values are the concentrations from that experiment
        concentrations_dict = dict(zip(coloc_collated2.ID_length_conc_image, coloc_collated2.concentration))
        #now make a column in counts df, which is the length of the fibril each of those hsps was coloc with
        counts['lengths (pixel)'] = counts['ID_length_conc_image'].map(lengths)
        #now map on the concentration matching each hsp onto that df
        counts['concentration'] = counts['ID_length_conc_image'].map(concentrations_dict)
        # same for the protein
        counts['protein']= counts['ID_length_conc_image'].map(proteins)

        # with this new dataframe, divide the count of foci, by the length in pixels for each contour
        counts['foci_per_pixel'] = counts['hsp_count']/counts['lengths (pixel)']

        alls.append(counts)

    alls=pd.concat(alls)

    #find the median of this
    median_foci_per_pixel = counts['foci_per_pixel'].median()
    alls.to_csv(f'{output_folder}foci_per_length_unit.csv')
    #plot as a violinplot for preliminary data vis.
    ax = sns.violinplot(x='concentration', y='log_foci', data=alls, scale='width', palette='Reds', hue='protein')
    plt.title(f'{protein}foci per length unit fibril (pixel)')
    plt.ylabel(f'{protein}exp(foci per pixel)')
    plt.ylim(1.0,1.8)
    plt.savefig(f'{output_folder}HSPA8_foci_per_length_unit_colocalised.png')
    plt.savefig(f'{output_folder}HSPA8_foci_per_length_unit_colocalised.svg')
    plt.show()



#---------------------------------PROPER PLOTTING SECTION
#for plotting multiple treatments on the same graph from the % colocalisation dataframes saved before
input_folder='python_results/Colocalisation/'
output_folder='python_results/Colocalisation/'
experiment_type=[name for name in os.listdir(input_folder) ]
experiment_type=[name for name in experiment_type if not 'png' in name ]
experiment_type=[name for name in experiment_type if not '.csv' in name ]
#GATHER all of the output files you generated for each treatment above
concat=[]
for name in experiment_type:
    filepath=f'{input_folder}{name}/3_foci_per_pixel/foci_per_length_unit.csv'
    test=pd.read_csv(filepath)
    test['treatment']=name
    concat.append(test)
#smoosh them together and you have a big df to plot from (which has treatments, concentrations, proteins to plot by)
concat=pd.concat(concat)

#make a dictionary to shorten names for plotting
dicto = {
        'Experiment100-FC1-1_HSPA8_A8-B1-110': '+110-0.5uM',
        'Experiment100-FC2-1_HSPA8_A8-B1-SOD1': '+SOD-2h'
        }
#define the order you want to plot (by its new name)
order_of_experiment = ['+110-0.5uM','+SOD-2h']
#map this on as a new col
concat['treatment_to_plot']=concat['treatment'].map(dicto)



#now plot as a boxplot the density of foci/pixel
ax=sns.boxplot(data=concat, x='treatment_to_plot', y='foci_per_pixel', palette='Reds', order=order_of_experiment)
ax.set_ylim(0,0.5)
ax.set_xlabel('Protein colocalised with fibrils/each other')
ax.set_ylabel('number of foci per pixel')
ax.set_title('Number foci per pixel colocalised by A8/JB1')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f'{output_folder}FOCI_per_pixel.png')
plt.savefig(f'{output_folder}FOCI_per_pixel.svg')
plt.show()



