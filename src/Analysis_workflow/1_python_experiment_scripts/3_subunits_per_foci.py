import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#This script gathers, organizes and plots the data from py4bleaching 

input_folder='python_results/py4bleaching/'
output_folder='python_results/subunits_per_foci/'
#find stoichiometry files output from py4bleaching
stoich_files =[[f'{root}/{filename}' for filename in files if 'molecule_counts.csv' in filename] for root, dirs, files in os.walk(f'{input_folder}')]
stoich_files=[item for sublist in stoich_files for item in sublist ]

if not os.path.exists(output_folder):
        os.makedirs(output_folder)


def gather_coloc_stoich(stoich_files):
    """gathers all of the molecule counts files, makes sure the paths make sense. gets some metadata (variable2, which is concentration, so we can split between treatments later/ identify things more easily), and assigns this variable to a column. 

    Args:
        stoich_files (list): list of molecule_counts.csv files

    Returns:
        list: all files added to list for concatination later
    """
    molecule_counts=[]
    for filepath in stoich_files:
        if 'noncoloc' not in filepath:
            filepath=filepath.replace('\\', '/')
            data= pd.read_csv(filepath)
            variable2s = data.variable2
            concentrations=[]
            for item in variable2s:
                concentration=item.split('-')[0]
                concentrations.append(concentration)
            data['concentration']=concentrations
            molecule_counts.append(data)
    return molecule_counts


molecule_counts=gather_coloc_stoich(stoich_files)
#concatinate into a dataframe
molecule_counts=pd.concat(molecule_counts)
#define order you want it plotted


#dictionary to convert long names to their treatment for easy plotting
concentration_dict={
    
    'Experiment100b-fc1-1': 'A8+JB1+110-5nM',
    'Experiment100b-fc1-2': 'A8+JB1+110-0.5uM',
    'Experiment100b-fc2-1': 'A8+JB1+SOD1'

    }
molecule_counts['concentration']=molecule_counts['Experiment_number'].map(concentration_dict)

#drop columns you don't want
molecule_counts.drop([col for col in molecule_counts.columns.tolist() if ' ' in col], axis=1, inplace=True)
molecule_counts.drop([col for col in molecule_counts.columns.tolist() if 'last_step_small_mol_count' in col], axis=1, inplace=True)
molecule_counts.drop([col for col in molecule_counts.columns.tolist() if 'single_step_mol_count' in col], axis=1, inplace=True)
molecule_counts.drop([col for col in molecule_counts.columns.tolist() if 'max_fluorescence' in col], axis=1, inplace=True)

#define the order to plot
order_of_experiment = ['A8+JB1+110-5nM', 'A8+JB1+110-0.5uM', 'A8+JB1+SOD1']
#plot molecule size as number of subunits per molecule, in a boxplot
ax = sns.boxplot(x="concentration",
 y="all_small_mol_count", 
 data=molecule_counts,  
 order=order_of_experiment,
  #scale='width',
   palette='YlOrBr')
plt.xticks(rotation=45)
ax.set_ylabel('# of subunits')
ax.set_xlabel('Concentration (nM)')
ax.set_ylim(0,30)
plt.tight_layout()
plt.title(f'# of HSPA8 molecules per foci')
plt.savefig(f'{output_folder}number_subunits_per_foci.png')
plt.show()

