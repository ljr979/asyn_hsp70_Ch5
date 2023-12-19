import os, re
import threading
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger


#this script aims to take the data obtained by running the end vs middle script and plot it (that is, looking at the binding of chaperones to the end pixels vs. the middle of the fibril)
#this version adds loop to go over multiple proteins rather than multiple concentrations and adds in the size based on py4bleaching analysis (i.e., number of subunits per molecule, per location)

#NOTE: this script will only work if you ran the 0_add_coords_save_new_copy.py script in the beginning, before running any py4bleaching analysis
 
input_folder='python_results/py4bleaching/'
output_folder=f'python_results/end_vs_middle/'
current_output=f'{output_folder}step_3/'
if not os.path.exists(current_output):
         os.makedirs(current_output)

#adding a little bit at the start here that combines all of the molecule counts from the py4bleaching I had did on this experiment. But, for this script to work I kind of need them to be togehter so I just read them in and concatinate then save again

counto =[[f'{root}/{filename}' for filename in files if 'molecule_counts.csv' in filename] for root, dirs, files in os.walk(f'{input_folder}')]
counto=[item for sublist in counto for item in sublist ]
combine=[]
#read in all the counts (in case imaging for the experiment had different conditions, and as a result you ahd to analyse photobleaching separately) and then combine them into one  df
for filename in counto:
    d=pd.read_csv(f'{filename}')
    combine.append(d)
combine=pd.concat(combine)

#now save this
new_input=f'{input_folder}all_combined/Trajectories/'
if not os.path.exists(new_input):
         os.makedirs(new_input)
combine.to_csv(f'{new_input}molecule_counts.csv')

#this is a string in the name of the files, which you want to filter OUT. so in this example, any files that have 'jb1' won't be included, can change when you then want to filter out HSPA8 and only look at jb1
filters='jb1'

def get_files(input_folder, filters):
    input_folder_counts=f'{input_folder}Trajectories/'
    treatments=[treatment for treatment in os.listdir(input_folder_counts)]
    #path to molecule counts files
    collated_molecule_counts_paths=[]
    for treatment in treatments: 
        filepaths=[[f'{root}/{name}' for name in files if 'molecule_counts.csv' in name]for root, dirs, files in os.walk(f'{input_folder_counts}/')]
        filepaths=[item for sublist in filepaths for item in sublist]
        collated_molecule_counts_paths.append(filepaths)

    #make sure they're only read in once and no data double up
    path=[]
    for x in collated_molecule_counts_paths:
        if x not in path:
            path.append(x)
    collated_molecule_counts_paths=path 
    #flatten list
    filepaths=[item for sublist in collated_molecule_counts_paths for item in sublist]
   
    #read in the the files with the counts
    collated_molecule_counts=[]    
    for filepath in filepaths:
        df=pd.read_csv(filepath)
        collated_molecule_counts.append(df)
    all_molecule_counts=pd.concat(collated_molecule_counts)

    #find the ends data you generated in step one
    input_folder_ends=f'{output_folder}step_one/HSPA8/'
    #filter only for the protein you want to look at here
    end_vs_middle_files = [filename for filename in os.listdir(input_folder_ends) if f'{filters}' not in filename]

    #collate if there are multiple, smoosh together
    collated_end_middle=[]
    for filename in end_vs_middle_files:
        df=pd.read_csv(f'{input_folder_ends}{filename}')
        collated_end_middle.append(df)
    collated_end_middle=pd.concat(collated_end_middle)

    return treatments, all_molecule_counts, collated_end_middle


def end_or_middle_all(end_vs_middle):
    end_or_middle=[]
    for row, df in end_vs_middle.groupby('new_ID_hspX_hspY_num'):
        row
        #for each of the fibrils that have been analysed for their end and middle colocalisations (in the previous script), we want to split their unique names up into separate columns to map the trajectories to later
        df[['Contour_ID','coordsX','coordsY', 'number']]=df['new_ID_hspX_hspY_num'].str.split('_', expand = True)
        df['Contour_ID']=df['Contour_ID'].astype(float)
        cols=['Contour_ID','coordsX','coordsY']
        #we are also popping off the enumerated part so that we can match this more easily to the  coords
        df['new_ID_hspX_hspY'] = df[cols].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
        end_or_middle.append(df)
    end_or_middle=pd.concat(end_or_middle)
    #now make a dictionary which usess the identifier of the molecule, and the allocation to end or middle, as the value.
    ends_dict=dict(zip(end_or_middle.new_ID_hspX_hspY, end_or_middle.Where))
    return ends_dict


def molecule_counts_format(molecule_counts, ends_dict):
    cols=['Contour_ID','coordsX','coordsY']
    #now we want to join together for the molecule counts dataframe (ie. the py4bleaching output) the columns that are identifiers and their locations, so that they match exactly to the ones from the fibrils (end v middle) output.
    molecule_counts['new_ID_hspX_hspY'] = molecule_counts[cols].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
    #now map them together so that the end and middle comes into the molecule counts dataframe 
    molecule_counts['end_or_middle']=molecule_counts['new_ID_hspX_hspY'].map(ends_dict)
    #drop unnecessary columns
    molecule_counts=drop_stuff(molecule_counts=molecule_counts)
    #filter / separate end and middle molecules
    Ends=molecule_counts[molecule_counts['end_or_middle']=='END']
    Mids=molecule_counts[molecule_counts['end_or_middle']=='MIDDLE']
    return molecule_counts, Ends, Mids


def drop_stuff(molecule_counts):
    molecule_counts.drop([col for col in molecule_counts.columns.tolist() if ' ' in col], axis=1, inplace=True)
    molecule_counts.drop([col for col in molecule_counts.columns.tolist() if 'all_small_mol_count' in col], axis=1, inplace=True)
    molecule_counts.drop([col for col in molecule_counts.columns.tolist() if 'single_step_mol_count' in col], axis=1, inplace=True)
    molecule_counts.drop([col for col in molecule_counts.columns.tolist() if 'max_fluorescence' in col], axis=1, inplace=True)
    molecule_counts.drop([col for col in molecule_counts.columns.tolist() if 'colocalisation' in col], axis=1, inplace=True)
    return molecule_counts


def plot_violin_count(molecule_counts, output_folder, dicto, order_of_experiment):
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
    plt.ylabel(f'# of molecules/foci')
    plt.xlabel('Experimental condition')
    plt.title(f'# of molecules/foci @ end or middle ')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f'{output_folder}/end_middle_molecules_per_foci.png')
    plt.savefig(f'{output_folder}/end_middle_molecules_per_foci.svg')
    plt.show()


def calculate_percent_summary(molecule_counts, output_folder):
    coloc_dict=[]
    #now make a summary 
    #loop through each treatment
    for treatment, df in molecule_counts.groupby('variable1'):
        LIST=[]
        #total number of hsps that are bound to a fibril
        total_hsp_fibs=len(df)
        #total num molecules at the end of the fibril theyre bound to
        total_ends=len(df[df['end_or_middle']=='END'])
        #total at the middle
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

#figure out where each molecule is
ends_dict=end_or_middle_all(end_vs_middle=end_vs_middle)

#map on the location of the molecule to the trajectory it matches (molecule size)
molecule_counts, Ends, Mids=molecule_counts_format(molecule_counts, ends_dict)
#save!
molecule_counts.to_csv(f'{current_output}/end_vs_middle_collated.csv')

#assign nice names for plotting
dicto={
    
    'a8-b1-110': 'A8+JB1+110',
    'a8-b1-sod1': 'A8+JB1+SOD1'

    }

#define the order of hte plotting
order_of_experiment = ['A8+JB1+110', 'A8+JB1+SOD1']

#plot as violinplots
plot_violin_count(molecule_counts=molecule_counts, output_folder=output_folder, dicto=dicto, order_of_experiment=order_of_experiment)

#make a summary of the percent of total coloc molecules at each location
coloc_dict=calculate_percent_summary(molecule_counts, output_folder)

#put the above data int oa dataframe and then plot it as a bar chart (this equates to the TOTAL NUMBER OF HSPA8 molecules that are colocalised with fibrils and whether they have a preference for end or middle)
conc_plot=pd.melt(coloc_dict, id_vars=['treatment' ], value_vars=['percentage_end','percentage_middle'], var_name="end or middle", value_name='percentage of total chaperones at each location')
f=sns.barplot(data=conc_plot, x='treatment', y='percentage of total chaperones at each location', hue='end or middle', palette='rocket')
plt.title('Percentage of colocalised HSPA8 @ end vs middle ')
plt.savefig(f'{current_output}/percent_summary_endmiddle.svg')
plt.savefig(f'{current_output}/percent_summary_endmiddle.png')
plt.show()



