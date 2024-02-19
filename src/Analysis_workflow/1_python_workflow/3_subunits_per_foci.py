"""This script gathers, organizes and plots the data from py4bleaching 

"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def gather_coloc_stoich(stoich_files):
    """gathers all of the molecule counts files, makes sure the paths make sense. gets some metadata (variable2, which is concentration, so we can split between treatments later/ identify things more easily), and assigns this variable to a column. 

    Args:
        stoich_files (list): list of molecule_counts.csv files

    Returns:
        list: all files added to list for concatination later
    """
    molecule_counts=[]
    for filepath in stoich_files:
             
        filepath=filepath.replace('\\', '/')
        data= pd.read_csv(filepath)
        if 'non-coloc' not in filepath:
             data['colocalisation']='coloc'
        if 'datatype' in data.columns.tolist():
             data.drop([col for col in data.columns.tolist() if 'datatype' in col], axis=1, inplace=True)
        variable2s = data.variable2
        concentrations=[]
        for item in variable2s:
            concentration=item.split('-')[0]
            concentrations.append(concentration)
        data['concentration']=concentrations

        molecule_counts.append(data)
    return molecule_counts

def plot_violins( df, order_of_experiment, x, y, palette):
    """plots the number of subunits per colocalised foci. As a violinplot.

    Args:
        output_folder (str): wher eto save
        df (df): dataframe with concatinated molecule count
        order_of_experiment (list): list of treatments in order for plotting
        x (str): the variable (column name) to plot on the x axis
        y (str): volumn variable to plot on the y axis
        palette (str): colour palette for plotting
        ylim (int): max on y axis
    """
    ax = sns.boxplot(x=x,
    y=y, 
    data=df,  
    order=order_of_experiment,
    palette=palette)
    plt.xticks(rotation=45)
    ax.set_ylabel('# of subunits')
    ax.set_xlabel('Treatment')
    #ax.set_ylim(ylim)
    plt.tight_layout()
    plt.title(f'# of HSPA8 molecules per foci')
    plt.show()

if __name__ == "__main__":
    
    input_folder='data/Analysis_workflow/2_example_python_output/1_py4bleaching/Trajectories/'
    output_folder='data/Analysis_workflow/2_example_python_output/3_subunits_per_foci/'

    #find stoichiometry files output from py4bleaching
    stoich_files =[[f'{root}/{filename}' for filename in files if 'molecule_counts.csv' in filename] for root, dirs, files in os.walk(f'{input_folder}')]
    stoich_files=[item for sublist in stoich_files for item in sublist]
    stoich_files= [x for x in stoich_files if 'all_combined' not in x]

    if not os.path.exists(output_folder):
            os.makedirs(output_folder)

    molecule_counts=gather_coloc_stoich(stoich_files)
    #concatinate into a dataframe
    molecule_counts=pd.concat(molecule_counts)
    #account for all molecules that are < 1 due to variation in step size
    molecule_counts['all_small_mol_count']=molecule_counts['all_small_mol_count']+1
    #dictionary to convert long names to their treatment for easy plotting
    concentration_dict={
        'Experiment98-fc1-1': 'Enum-treat1',
        'Experiment98-fc1-2': 'Enum-treat2'

        }
    molecule_counts['concentration']=molecule_counts['Experiment_number'].map(concentration_dict)
    order_of_experiment = ['Enum-treat1', 'Enum-treat2']

    #drop columns you don't want
    drop_list=['last_step_small_mol_count','single_step_mol_count','max_fluorescence', 'Unnamed: 0', 'Unnamed: 0.1']
    molecule_counts.drop([col for col in molecule_counts.columns.tolist() if col in drop_list], axis=1, inplace=True)

    #plot molecule size as number of subunits per molecule, in a boxplot
    plot_violins(molecule_counts, order_of_experiment, x='concentration', y="all_small_mol_count", palette="YlOrBr")

