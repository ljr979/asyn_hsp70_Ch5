"""This script aims to take the data obtained by running the end vs middle script and plot it (that is, looking at the binding of chaperones to the end pixels vs. the middle of the fibril)
This version adds loop to go over multiple proteins rather than multiple concentrations and adds in the size based on py4bleaching analysis (i.e., number of subunits per molecule, per location)

"""
import os, re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#NOTE: this script will only work if you ran the 0_add_coords.py script in the beginning, before running any py4bleaching analysis.

def get_files(input_folder, filters):
    """colect molecule counts files

    Args:
        input_folder (_type_): _description_
        filters (str): any files you want to ignore that are also within the folder (e.g. if there are two proteins and you want to include only one)

    Returns:
        _type_: _description_
    """
    input_folder_counts = f'{input_folder}Trajectories/coloc/'
    treatments = [treatment for treatment in os.listdir(input_folder_counts)]
    #path to molecule counts files
    collated_molecule_counts_paths = []
    for treatment in treatments: 
        filepaths = [[f'{root}/{name}' for name in files if 'molecule_counts.csv' in name]for root, dirs, files in os.walk(f'{input_folder_counts}/')]
        filepaths = [item for sublist in filepaths for item in sublist]
        collated_molecule_counts_paths.append(filepaths)

    #make sure they're only read in once and no data double up
    path = []
    for x in collated_molecule_counts_paths:
        if x not in path:
            path.append(x)
    collated_molecule_counts_paths = path 
    #flatten list
    filepaths = [item for sublist in collated_molecule_counts_paths for item in sublist]
   
    #read in the the files with the counts
    collated_molecule_counts = []    
    for filepath in filepaths:
        df = pd.read_csv(filepath)
        collated_molecule_counts.append(df)
    all_molecule_counts = pd.concat(collated_molecule_counts)

    #find the ends data you generated in step one
    input_folder_ends = f'{output_folder}step_one/'
    #filter only for the protein you want to look at here
    end_vs_middle_files = [filename for filename in os.listdir(input_folder_ends) if f'{filters}' not in filename]

    #collate if there are multiple, smoosh together
    collated_end_middle = []
    for filename in end_vs_middle_files:
        df = pd.read_csv(f'{input_folder_ends}{filename}')
        collated_end_middle.append(df)
    collated_end_middle = pd.concat(collated_end_middle)

    return treatments, all_molecule_counts, collated_end_middle


def end_or_middle_all(end_vs_middle):
    """loops over end v middle dataframe and groups on the unique IDS. split these names into separate columns and gives a dictionary which assigns the unique coordinate/name, as end or middle, to map onto some other dataframe later.

    Args:
        end_vs_middle (df): df containing end v middle data

    Returns:
        dict: dictionary with the ends, middles and matching unique names
    """
    end_or_middle = []
    for row, df in end_vs_middle.groupby('new_ID_hspX_hspY_num'):
        row
        #for each of the fibrils that have been analysed for their end and middle colocalisations (in the previous script), we want to split their unique names up into separate columns to map the trajectories to later
        df[['Contour_ID','coordsX','coordsY', 'number']] = df['new_ID_hspX_hspY_num'].str.split('_', expand = True)
        df['Contour_ID'] = df['Contour_ID'].astype(float)
        cols = ['Contour_ID','coordsX','coordsY']
        #we are also popping off the enumerated part so that we can match this more easily to the  coords
        df['new_ID_hspX_hspY'] = df[cols].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
        end_or_middle.append(df)
    end_or_middle = pd.concat(end_or_middle)
    #now make a dictionary which usess the identifier of the molecule, and the allocation to end or middle, as the value.
    ends_dict = dict(zip(end_or_middle.new_ID_hspX_hspY, end_or_middle.Where))
    return ends_dict


def molecule_counts_format(molecule_counts, ends_dict):
    """format the molecule counts dataframe so each 'molecule count' has also got the end or middle location

    Args:
        molecule_counts (df): df with the subunit count for each colocalised mol
        ends_dict (dict): dictoinary with the unique name (coords and fibril contour id) matching the end v middle location

    Returns:
        df: dfs of ends, mids, and together
    """
    cols = ['Contour_ID','coordsX','coordsY']
    #now we want to join together for the molecule counts dataframe (ie. the py4bleaching output) the columns that are identifiers and their locations, so that they match exactly to the ones from the fibrils (end v middle) output.
    molecule_counts['new_ID_hspX_hspY'] = molecule_counts[cols].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
    #now map them together so that the end and middle comes into the molecule counts dataframe 
    molecule_counts['end_or_middle'] = molecule_counts['new_ID_hspX_hspY'].map(ends_dict)
    #drop unnecessary columns
    molecule_counts = drop_stuff(molecule_counts=molecule_counts)
    #filter / separate end and middle molecules
    Ends = molecule_counts[molecule_counts['end_or_middle']=='END']
    Mids = molecule_counts[molecule_counts['end_or_middle']=='MIDDLE']
    return molecule_counts, Ends, Mids


def drop_stuff(molecule_counts):
    """drop these columns

    Args:
        molecule_counts (df): df which needs columsn dropped

    """
    molecule_counts.drop([col for col in molecule_counts.columns.tolist() if ' ' in col], axis=1, inplace=True)
    molecule_counts.drop([col for col in molecule_counts.columns.tolist() if 'all_small_mol_count' in col], axis=1, inplace=True)
    molecule_counts.drop([col for col in molecule_counts.columns.tolist() if 'single_step_mol_count' in col], axis=1, inplace=True)
    molecule_counts.drop([col for col in molecule_counts.columns.tolist() if 'max_fluorescence' in col], axis=1, inplace=True)
    molecule_counts.drop([col for col in molecule_counts.columns.tolist() if 'colocalisation' in col], axis=1, inplace=True)
    return molecule_counts


def plot_box_count(molecule_counts, dicto, order_of_experiment, ylim):
    """plots subunit count at each region as a boxplot

    Args:
        molecule_counts (df): df with mol count and location
        output_folder (str): where to save
        dicto (dict): dictionary to map the variable 1 to the label you want to use for plotting
        order_of_experiment (list): order to plot
    """
    molecule_counts['treatment_to_plot']=molecule_counts['variable1'].map(dicto)


    ax = sns.boxplot(
        x='treatment_to_plot',
        y= 'last_step_mol_count', 
        data=molecule_counts, 
        hue='end_or_middle', 
        palette='Greens', 
        order=order_of_experiment)
    
    ax.legend(loc='upper center', ncol=2)
    ax.set_ylim(0,ylim)
    plt.ylabel(f'# of molecules/foci')
    plt.xlabel('Experimental condition')
    plt.title(f'# of molecules/foci @ end or middle ')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    input_folder = 'data/Analysis_workflow/2_example_python_output/1_py4bleaching/'
    output_folder = f'data/Analysis_workflow/2_example_python_output/5_end_vs_middle_init/'

    current_output = f'data/Analysis_workflow/2_example_python_output/7_subunits_end_middle/step_3/'

    if not os.path.exists(current_output):
            os.makedirs(current_output)


    #read in the 'molecule_counts' file from this treatment.
    treatments, molecule_counts, end_vs_middle = get_files(input_folder=f'{input_folder}/', filters='jb1')
    molecule_counts.drop([col for col in molecule_counts.columns.tolist() if 'Unnamed: 0' in col], axis=1, inplace=True)

    #figure out where each molecule is
    ends_dict = end_or_middle_all(end_vs_middle=end_vs_middle)

    #map on the location of the molecule to the trajectory it matches (molecule size)
    molecule_counts, Ends, Mids = molecule_counts_format(molecule_counts, ends_dict)
    #save!
    molecule_counts.to_csv(f'{current_output}/end_vs_middle_collated.csv')

    #assign nice names for plotting
    dicto = {
        
        'a8-b1-flowout': 't1-t2-t3',
        'a8-b1-sod1': 't3-t4-t5'
      

        }

    #plot as violinplots
    plot_box_count(molecule_counts=molecule_counts, output_folder=current_output, dicto=dicto, order_of_experiment=['t1-t2-t3', 't3-t4-t5' ])




