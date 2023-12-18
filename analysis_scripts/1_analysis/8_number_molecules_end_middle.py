import os, re
import threading
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger


#this script aims to take the data obtained by running the end vs middle R script and plotting it (that is, looking at the binding of chaperones to the end 3 pixels vs. the middle of the fibril)
#this version adds loop to go over multiple proteins rather than multiple concentrations
input_folder='python_results/py4bleaching/'
output_folder=f'python_results/end_vs_middle/'
current_output=f'{output_folder}step_3/'
if not os.path.exists(current_output):
         os.makedirs(current_output)

#addinga little bit at the start here that combines all of the molecule counts from the two different rounds of py4bleaching I had to do on this experiment. for some reason they were wildly different so  I ran them separately. But, for this script to work I kind of need them to be togehter so I just read them in and concatinate then save again

counto =[[f'{root}/{filename}' for filename in files if 'molecule_counts.csv' in filename] for root, dirs, files in os.walk(f'{input_folder}')]
counto=[item for sublist in counto for item in sublist ]
combine=[]
for filename in counto:
    d=pd.read_csv(f'{filename}')
    combine.append(d)
combine=pd.concat(combine)

new_input=f'{input_folder}all_combined/Trajectories/'
if not os.path.exists(new_input):
         os.makedirs(new_input)
combine.to_csv(f'{new_input}molecule_counts.csv')


def get_files(input_folder):
    input_folder_counts=f'{input_folder}Trajectories/'
    treatments=[treatment for treatment in os.listdir(input_folder_counts)]


    collated_molecule_counts_paths=[]
    for treatment in treatments: 
        filepaths=[[f'{root}/{name}' for name in files if 'molecule_counts.csv' in name]for root, dirs, files in os.walk(f'{input_folder_counts}/')]
        filepaths=[item for sublist in filepaths for item in sublist]
        collated_molecule_counts_paths.append(filepaths)


    path=[]
    for x in collated_molecule_counts_paths:
        if x not in path:
            path.append(x)
    collated_molecule_counts_paths=path 
    filepaths=[item for sublist in collated_molecule_counts_paths for item in sublist]
   

    collated_molecule_counts=[]    
    for filepath in filepaths:
        df=pd.read_csv(filepath)
        
        collated_molecule_counts.append(df)
    all_molecule_counts=pd.concat(collated_molecule_counts)

    input_folder_ends=f'{output_folder}step_one/HSPA8/'
    end_vs_middle_files = [filename for filename in os.listdir(input_folder_ends) if 'jb1' not in filename]
    collated_end_middle=[]
    for filename in end_vs_middle_files:
        df=pd.read_csv(f'{input_folder_ends}{filename}')
        # treatment=filename.split('_')[0:3]
        # treatment='_'.join(treatment)
        # df['treatment']=treatment
        collated_end_middle.append(df)

    collated_end_middle=pd.concat(collated_end_middle)
    return treatments, all_molecule_counts, collated_end_middle


def drop_stuff(molecule_counts):
    molecule_counts.drop([col for col in molecule_counts.columns.tolist() if ' ' in col], axis=1, inplace=True)
    molecule_counts.drop([col for col in molecule_counts.columns.tolist() if 'all_small_mol_count' in col], axis=1, inplace=True)
    molecule_counts.drop([col for col in molecule_counts.columns.tolist() if 'single_step_mol_count' in col], axis=1, inplace=True)
    molecule_counts.drop([col for col in molecule_counts.columns.tolist() if 'max_fluorescence' in col], axis=1, inplace=True)
    molecule_counts.drop([col for col in molecule_counts.columns.tolist() if 'colocalisation' in col], axis=1, inplace=True)
    return molecule_counts


def plot_violin_count(molecule_counts, output_folder, dicto, order_of_experiment):
#plot distribution of fibril lengths cooc and not coloc
    #for treatment, df in molecule_counts.groupby('treatment'): 
        #fig, ax = plt.subplots()
    
    molecule_counts['treatment_to_plot']=molecule_counts['variable1'].map(dicto)
    order_of_experiment= order_of_experiment


    ax = sns.violinplot(
        x='treatment_to_plot',
        y= 'last_step_mol_count', 
        data=molecule_counts, 
        scale='width', 
        hue='end_or_middle', 
        palette='Greens', 
        order=order_of_experiment,
        alpha=0.45)
    ax.legend(loc='upper center', ncol=2)
    ax.set_ylim(0,25)
    # ax.annotate(f"median colocalised length = {median_length_colocal_fibrils}nm", xy=(0.1, 0.8), xycoords='figure fraction')
    plt.ylabel(f'# of molecules/foci')
    plt.xlabel('Experimental condition')
    #plt.show()
    
    plt.title(f'# of molecules/foci @ end or middle ')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f'{output_folder}/end_middle_molecules_per_foci.png')
    plt.savefig(f'{output_folder}/end_middle_molecules_per_foci.svg')
    plt.show()




def end_or_middle_all(end_vs_middle):

    end_or_middle=[]
    for row, df in end_vs_middle.groupby('new_ID_hspX_hspY_num'):
        row
        
        df[['Contour_ID','coordsX','coordsY', 'number']]=df['new_ID_hspX_hspY_num'].str.split('_', expand = True)
        df['Contour_ID']=df['Contour_ID'].astype(float)
        cols=['Contour_ID','coordsX','coordsY']
        df['new_ID_hspX_hspY'] = df[cols].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
        end_or_middle.append(df)
    end_or_middle=pd.concat(end_or_middle)

    ends_dict=dict(zip(end_or_middle.new_ID_hspX_hspY, end_or_middle.Where))
    return ends_dict


def molecule_counts_format(molecule_counts, ends_dict):
    #molecule_counts['Contour_ID']=molecule_counts['Contour_ID'].astype(int)
    cols=['Contour_ID','coordsX','coordsY']
    molecule_counts['new_ID_hspX_hspY'] = molecule_counts[cols].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
        
    molecule_counts['end_or_middle']=molecule_counts['new_ID_hspX_hspY'].map(ends_dict)

    #map this molecule count dictionary onto the end vs middle dataframe
    molecule_counts=drop_stuff(molecule_counts=molecule_counts)
    Ends=molecule_counts[molecule_counts['end_or_middle']=='END']
    Mids=molecule_counts[molecule_counts['end_or_middle']=='MIDDLE']
    return molecule_counts, Ends, Mids


def calculate_percent_summary(molecule_counts, output_folder):
    coloc_dict=[]
    #
    for treatment, df in molecule_counts.groupby('variable1'):
        LIST=[]
        total_hsp_fibs=len(df)
        total_ends=len(df[df['end_or_middle']=='END'])
        total_mids=len(df[df['end_or_middle']=='MIDDLE'])

        percent_end=total_ends/total_hsp_fibs*100

        percent_mid=total_mids/total_hsp_fibs*100

        LIST=treatment, percent_end, percent_mid
        coloc_dict.append(LIST)

    
    coloc_dict=pd.DataFrame(coloc_dict, columns=['treatment', 'percentage_end', 'percentage_middle'])
    coloc_dict.to_csv(f'{current_output}/Summary_percentage_end_middle.csv')
    return coloc_dict

    
treatments, all_molecule_counts, collated_end_middle=get_files(input_folder=f'{input_folder}all_combined/')

#first step read in end vs. middle data WITH the column which has the contour ID _ X COORD_YCOORD so that this name matches those associated with the molecule count from py4bleaching. should also have the END or MIDDLE definition
end_vs_middle=collated_end_middle

#read in the 'molecule_counts' file from this treatment.
molecule_counts=all_molecule_counts
molecule_counts.drop([col for col in molecule_counts.columns.tolist() if 'Unnamed: 0' in col], axis=1, inplace=True)
ends_dict=end_or_middle_all(end_vs_middle=end_vs_middle)

molecule_counts, Ends, Mids=molecule_counts_format(molecule_counts, ends_dict)
molecule_counts.to_csv(f'{current_output}/end_vs_middle_collated.csv')


dicto={
    
    'a8-b1-110': 'A8+JB1+110',
    'a8-b1-sod1': 'A8+JB1+SOD1'

    }

order_of_experiment = ['A8+JB1+110', 'A8+JB1+SOD1']

plot_violin_count(molecule_counts=molecule_counts, output_folder=output_folder, dicto=dicto, order_of_experiment=order_of_experiment)

coloc_dict=calculate_percent_summary(molecule_counts, output_folder)

#put the above data int oa dataframe and then plot it as a bar chart (this equates to the TOTAL NUMBER OF HSPA8 molecules that are colocalised with fibrils and whether they have a preference for end or middle)
conc_plot=pd.melt(coloc_dict, id_vars=['treatment' ], value_vars=['percentage_end','percentage_middle'], var_name="end or middle", value_name='percentage of total chaperones at each location')
f=sns.barplot(data=conc_plot, x='treatment', y='percentage of total chaperones at each location', hue='end or middle', palette='rocket')
plt.title('Percentage of colocalised HSPA8 @ end vs middle ')
plt.savefig(f'{current_output}/percent_summary_endmiddle.svg')
plt.savefig(f'{current_output}/percent_summary_endmiddle.png')
plt.show()



