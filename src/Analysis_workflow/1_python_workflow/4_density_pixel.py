"""Script to find the density of fluorescently labelled molecules bound per individual fibril (as molecule / pixel)

"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def organise_unique_fibs(input_folder, output_folder, c):
    """find fibrils, and assign each fibril a unique name. this is assigned to each pixel of the ridge which is colocalised with a chaperone molecule

    Args:
        input_folder (str): location of colocalisation files
        output_folder (str): location to save the new fibril colocalisation data
        c (_type_): this is looped over within a for loop, it should be the type of colocalisation we are looking at (i.e. fibrils + 647, or fibrils + 488). should be the name of the folder that the different colocalisation files are within

    """
    collated = []
    folder = f'{input_folder}{c}/'
            #gather data on treatment
    concentration = folder.split('/')[-3].split('_')[-1]
            #gather data of experiment
    experiment_number = folder.split('/')[-3].split('_')[0]
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
        protein = filepath.split('_')[0]
                #make new columns with identifiers about the experiment/treatment
        data['concentration'] = concentration
        data['protein'] = protein
        data['image_number'] = filepath
        colocalised_fibrils.append(data)
            #collate all of them
    colocalised_fibrils = pd.concat(colocalised_fibrils)
        
    collated.append(colocalised_fibrils)
    collated = pd.concat(collated)

    collated.rename(columns={'Contour ID':'Contour_ID'}, inplace=True)

            #now we have a big dataframe with only the colocalisation data and we can get the number of foci on each fibril from this
    collated2 = []
            #iterate over each 'fibril'. Every fibril has a row assigned to each individual pixel of that fibril, so we are going over each pixel, for each fibril and giving it a unique name with the length of that fibril, the image it came from, the treatment and the 'contour id' assigned at imagej
    for Contour_ID, row in collated.iterrows():
        Contour_ID
        row
        row['ID_length_conc_image'] = str(row['Contour_ID'])+'_'+str(row['Length'])+'_'+str(row['concentration'])+'_'+str(row['image_number'])
        row = pd.DataFrame(row).T
        collated2.append(row)

    collated2 = pd.concat(collated2)
            #t hen we filter this for only those rows that are colocalised with a molecule
    coloc_collated2 = collated2[collated2['distance']!=-1]
            #save this file to work on from here
    coloc_collated2.to_csv(f'{output_folder}{c}/fibril_collated-colocal-only_data.csv')
            #save the non-filtered for colocalisation version too.
    collated2.to_csv(f'{output_folder}{c}/HSPA8_fibril_collated-colocal_data.csv')
    return coloc_collated2

def count_coloc_spots(alls, coloc_collated2, proteins, lengths):
    """counts the actual number of colocalised spots per unique ridge name

    Args:
        alls (list): list to append the results to
        coloc_collated2 (df): df containing the unique ridges and their colocalisation data
        proteins (dict): dictionary to map the protein type to the colocalsied spot for each unique ID, to map back later. This is useful when there are multiple types of proteins whose colocalisation we are analysing
        lengths (dict): same as proteins, but the length

    """
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
    return counts,protein

def concatinate(df, files_list):
    """concatinate files in a list

    Args:
        df (list): empty list to put files into
        files_list (list): list of files to read in and concat

    Returns:
        _type_: _description_
    """
    for fc in files_list:
        path = fc.replace('//', '/')
        path = path.replace('\\', '/')
        test = pd.read_csv(fc)
        cond = path.split('/')[5]
        test['treatment'] = cond
        df.append(test)
    df = pd.concat(df)
    return df

def concat_density_files(input_folder, experiment_type):
    """find and concatinate the files that have the foci per length unit saved in them (For each fibril) i.e. each datapoint is one fibril

    Args:
        input_folder (str): folder path above the fibrils-647 or fibrils-488 folders
        experiment_type (list): list containing the folders to loop over (this will be the fibrils-647 etc.)


    """
    concat = []
    for name in experiment_type:
        filepath = f'{input_folder}{name}/3_foci_per_pixel/foci_per_length_unit.csv'
        test = pd.read_csv(filepath)
        test['treatment'] = name
        concat.append(test)
    #smoosh them together and you have a big df to plot from (which has treatments, concentrations, proteins to plot by)
    concat = pd.concat(concat)
    return concat

def boxplot(df, order_of_experiment, ylim):
    ax = sns.boxplot(data=df, x='treatment', y='foci_per_pixel', palette='Reds', order=order_of_experiment)
    ax.set_ylim(0,ylim)
    ax.set_xlabel('Protein colocalised with fibrils/each other')
    ax.set_ylabel('number of foci per pixel')
    ax.set_title('Number foci per pixel colocalised by A8/JB1')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":

    input_top = 'data/Analysis_workflow/2_example_python_output/0_collect_data/Colocalisation/'
    output_top = 'data/Analysis_workflow/2_example_python_output/4_density_pixel/'

    #find sub-treatments within colocalisation folder
    treatments = [folder for folder in os.listdir(input_top)]


    for exp in treatments: 

        input_folder = f'{input_top}{exp}/'
        #make new folder to save this data
        output_folder = f'{output_top}{exp}/3_foci_per_pixel/'
        #define which experiment to analyse (i.e., hspa8 and fibrils only , or multiple colocalisations)
        experiments = ['fibrils-647']
 
        alls = []
        #loop over the experiment folders
        for c in experiments:
            coloc_collated2 = organise_unique_fibs(input_folder, output_folder, c)
            #define two dictionaries
            proteins = {}
            lengths = {}
            #then make a new dataframe where I count the number of colocalised SPOTS on each contour, and then add in the matching fibril length to this counted number of foci
            #counts the number of individual fibril pixels
            counts, protein = count_coloc_spots(alls, coloc_collated2, proteins, lengths)

        alls = pd.concat(alls)

        #find the median of this
        median_foci_per_pixel = counts['foci_per_pixel'].median()
        alls.to_csv(f'{output_folder}foci_per_length_unit.csv')

    #---------------------------------PROPER PLOTTING SECTION
    #for plotting multiple treatments on the same graph from the % colocalisation dataframes saved before
    experiment_type = [name for name in os.listdir(output_top) ]
    experiment_type = [name for name in experiment_type if not 'png' in name ]
    experiment_type = [name for name in experiment_type if not '.csv' in name ]
    #GATHER all of the output files you generated for each treatment above
    concat = concat_density_files(output_top, experiment_type)

    #make a dictionary to shorten names for plotting
    dicto = {
            't1-t2-t3': '1-2-3',

            }

    #map this on as a new col
    concat['treatment'] = concat['concentration'].map(dicto)

    #now plot as a boxplot the density of foci/pixel
    boxplot(df=concat, order_of_experiment=['1-2-3'], ylim=0.5)



