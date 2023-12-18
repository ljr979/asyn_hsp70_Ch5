import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

#this is the path that contains the trajectories that are just a couple of frames (initial intensity)
input_folder='python_results/initial_fluo_data/'
output_folder='python_results/FC1/subunits_per_foci/avg-t0-size-nobleach/'
#this is the path to the folder containing output from py4bleaching
pythons='python_results/py4bleaching/FC1/coloc/Trajectories/'



    

#look for end_vs_middle.csv files
if not os.path.exists(output_folder):
        os.makedirs(output_folder)


#this first section grabs the trajectories that are not full bleaching traces, and makes the average of those values the 'initial intensity'. These are for the SAME fields of view that were included in the 1h timepoint in this flow cell. it then takes the step sizes made by py4bleaching using the trajectories gathered later in the experiment, and gets the 'initial' molecule size (at time 0 before adding hsp110) it does this for noncolocalised too.
yikes_=[[f'{root}/{name}' for name in files if 'colocal_traj' in name]for root, dirs, files in os.walk(f'{input_folder}/')]
yikes_=[item for sublist in yikes_ for item in sublist if 'jb1' not in item]


collated=[]
for item in yikes_:
        df=pd.read_csv(f'{item}')
        df=df.T.reset_index().rename(columns = {'index': 'molecule_number'}).tail(-1)
        #find average intensity of each spot
        df['Avg-intensity']=df.mean(axis=1)

        item=item.replace('\\', '/')

        c=item.split('/')[-2]

        df['colocalisation']=c
        collated.append(df)

collated=pd.concat(collated)

#filter to keep only columns in which the average intensity is NOT negative
collated = collated.loc[collated["Avg-intensity"] > 0 ]


#find step sizes

stepspath=[[f'{root}/{name}' for name in files if 'median_steps' in name]for root, dirs, files in os.walk(f'{pythons}/')]
stepspath=[item for sublist in stepspath for item in sublist][0]
steps=pd.read_csv(f'{stepspath}')
stepsize=steps[steps['step_type']=='last_step']
stepszie=stepsize['step_size'][1]


collated = collated.loc[collated["Avg-intensity"] >= stepszie ]

collated['approx_mol_size']=collated['Avg-intensity']/stepszie
collated.to_csv(f'{output_folder}avg-intensity-nobleach.csv')

ax = sns.violinplot(
    x='colocalisation',
    y= 'approx_mol_size', 
    data=collated, 
    scale='width', 
    #hue='end_or_middle', 
    palette='Greens', 
    #order=order_of_experiment,
    alpha=0.45)
ax.legend(loc='upper center', ncol=2)
ax.set_ylim(0,20)
# ax.annotate(f"median colocalised length = {median_length_colocal_fibrils}nm", xy=(0.1, 0.8), xycoords='figure fraction')
plt.ylabel(f'# of molecules/foci')
plt.xlabel('coloc vs. non-coloc')

plt.title(f'# of molecules/foci (from average initial intensity)')
#plt.xticks(rotation=45)
plt.tight_layout()
#plt.savefig(f'{output_folder}/avg-molecules_per_foci.png')
plt.show()



#this section includes the mols per foci from bleaching , and plots them together with the subunits plotted from above (average intensity of spots when first detected)
#here I have also not included the non-colocalised spots, as I hadn't analysed them yet for the other treatments (other than the averaged non bleached spots)

stoich_files =[[f'{root}/{filename}' for filename in files if 'molecule_counts.csv' in filename] for root, dirs, files in os.walk(f'{pythons}')]
stoich_files=[item for sublist in stoich_files for item in sublist ]
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
molecule_counts = molecule_counts.loc[molecule_counts["max_fluorescence"] >= stepszie ]
collated2=molecule_counts[['molecule_number', 'last_step_mol_count', 'variable1', 'variable2']]

collated=collated[['molecule_number', 'colocalisation', 'approx_mol_size']]

#collated=collated[collated['colocalisation']=='coloc']

#this section is making sure that the two dataframes have matching columns names before concatinating them (variable 1 should be the thing you want to plot separately i.e. different treatments. variable 2 is another, but I've been keeping it as the protein we are looking at)
collated=collated[['molecule_number', 'approx_mol_size']]
collated['variable2']='HSPA8'
collated['variable1'] = 'Experiment96-FC1-1_HSPA8_A8-B1-flowout'
collated=collated.rename(columns = {'approx_mol_size':'last_step_mol_count', 'colocalisation':'variable1'})
for_plots=pd.concat([collated, collated2])


#make sure the names of the variable 1 in the dataframe are good for plotting not long and annoying/impractical

plot_dict = {'Experiment96-FC1-1_HSPA8_A8-B1-flowout': 'A8-B1',
             'a8-b1-+sod1-2h': '+SOD-2h', 'a8-b1-+hsp110-0.5um': '+hsp110-0.5uM-2h'}
for_plots['variable_plot']=for_plots['variable1'].map(plot_dict)
#now decide the order of the experiment and match this strings to the new names you just mapped onto your df
order_of_experiment=['A8-B1', '+SOD-2h', '+hsp110-0.5uM-2h']

ax = sns.boxplot(
    x='variable_plot',
    y= 'last_step_mol_count', 
    data=for_plots, 
    #scale='width', 
    #hue='end_or_middle', 
    palette='Greens', 
    order=order_of_experiment,
    #alpha=0.45
    )
#ax.legend(loc='upper center', ncol=2)
ax.set_ylim(0,15)
# ax.annotate(f"median colocalised length = {median_length_colocal_fibrils}nm", xy=(0.1, 0.8), xycoords='figure fraction')
plt.ylabel(f'# of molecules/foci')
plt.xlabel('Condition')

plt.title(f'# of molecules/foci JB1 + A8 before and after + 0.5uM HSP110')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f'{output_folder}/molecules_per_foci_v2.png')
plt.show()
