import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
input_folder='python_results/py4bleaching/'
output_folder='python_results/subunits_per_foci/'

stoich_files =[[f'{root}/{filename}' for filename in files if 'molecule_counts.csv' in filename] for root, dirs, files in os.walk(f'{input_folder}')]
stoich_files=[item for sublist in stoich_files for item in sublist ]

if not os.path.exists(output_folder):
        os.makedirs(output_folder)

def gather_coloc_stoich(stoich_files):
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
molecule_counts=pd.concat(molecule_counts)

order_of_experiment = ['A8+JB1+110-5nM', 'A8+JB1+110-0.5uM', 'A8+JB1+SOD1']


concentration_dict={
    
    'Experiment100b-fc1-1': 'A8+JB1+110-5nM',
    'Experiment100b-fc1-2': 'A8+JB1+110-0.5uM',
    'Experiment100b-fc2-1': 'A8+JB1+SOD1'

    }
molecule_counts['concentration']=molecule_counts['Experiment_number'].map(concentration_dict)


molecule_counts.drop([col for col in molecule_counts.columns.tolist() if ' ' in col], axis=1, inplace=True)
molecule_counts.drop([col for col in molecule_counts.columns.tolist() if 'last_step_small_mol_count' in col], axis=1, inplace=True)
molecule_counts.drop([col for col in molecule_counts.columns.tolist() if 'single_step_mol_count' in col], axis=1, inplace=True)
molecule_counts.drop([col for col in molecule_counts.columns.tolist() if 'max_fluorescence' in col], axis=1, inplace=True)

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



#test using the log of mol counts
molecule_counts['log_count']=np.log10(molecule_counts['all_small_mol_count'])

ax = sns.violinplot(x="concentration",
 y="log_count", 
 data=molecule_counts,  
 order=order_of_experiment,
  scale='width',
   palette='YlOrBr')
plt.xticks(rotation=45)
ax.set_ylabel('log10(# of subunits)')
ax.set_xlabel('Treatment')
#ax.set_ylim(0,5)
plt.tight_layout()
plt.title(f'log10 # of HSPA8 molecules per foci')
plt.savefig(f'{output_folder}log10_number_subunits_per_foci.png')
plt.show()