"""This script purely looks at COLOCALISATION and LENGTH OF FIBRILS, no data on molecule size etc.

"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def get_experiment_info(input_folder, treatment):
    """sets parameters so that you can save the important info from the experiment to carry through

    Args:
        input_folder (str): colocalisation folder
        treatment (str): treatment

    Returns:
        strings: experiment information
    """
    folder = f'{input_folder}{treatment}/'
    concentration = treatment.split('/')[0].split('_')[-1]
    experiment_number= folder.split('/')[-3].split('_')[0]
    
    output = f'{output_folder}{treatment}/'

    if not os.path.exists(output):
            os.makedirs(output)
    return folder, concentration, experiment_number, output

def concatinate_coloc_data(files):
    """gathers all the fibril colocalisation files, fixes some columns, concatinates and returns all colocalised fibril output

    Args:
        files (list ): files that have fibril colocalisation data

    Returns:
        df: concatinated coloc data
    """
    colocalisation_output = []
    
    for filepath in files: 
        filepath
        #read in file
        data = pd.read_csv(f'{folder}{filepath}')
        if 'fibril' in filepath:
            cols=data.columns.tolist()
            #remove class column
            if 'Class' in cols:
                data = data.drop(['Class'], axis=1)
        #name protein colocalised with fibril, from experiment info
        protein=f'{filepath}'.split('_')[0]+'_'+f'{filepath}'.split('_')[1]
        #add column to say whether colocalised or non colocalised fibril
        data['colocalisation'] = data['distance'].apply(lambda x: 'True' if x > -1 else 'False')
        data['treatment']=treatment
        data['protein_colocalised']=protein

        colocalisation_output.append(data)



    return colocalisation_output

def concatinate_coloc_hsps(files):
    """concatinate the colocalisation chaperones files

    Args:
        files (list): list of hsp files

    Returns:
        df: concatinated hsp files
    """
    colocalisation_output = []
    
    for filepath in files: 
        data = pd.read_csv(f'{filepath}')

        filename=f'{filepath}'.split('/')[-1]
        protein_colocalised=f'{filename}'.split('_')[0]+'_'+f'{filename}'.split('_')[1]
        data['colocalisation'] = data['distance'].apply(lambda x: 'True' if x > -1 else 'False')
        data['treatment']=f'{filepath}'.split('/')[-2]
        data['protein_colocalised']=protein_colocalised
        colocalisation_output.append(data)

    colocalisation_output=pd.concat(colocalisation_output)
    
    return colocalisation_output

def find_coloc_fibs_filter(colocalisation_output_fibrils, median_seeds_length, protein):
    """function to find the length of fibrils and seeds, and to find the percentage of the fibrils that are colocalised with AF647 labelled molecules.

    Args:
        colocalisation_output_fibrils (df): dataframe with the colocalisation fibrils information
        median_seeds_length (int): length of seeds (nm)
        protein (str): the protein being analysed
    """
    #filter the dataset to find how many 'contours' are fibrils vs seeds (unique ID'S as there are duplicates in the coloc output df) 
    #filter bigger than average seed length (nm)  
    bigger_than_seeds = colocalisation_output_fibrils['real_length_nm']>median_seeds_length
    fibrils_only = colocalisation_output_fibrils[bigger_than_seeds]
    #split into fibrils and seeds
    seeds = colocalisation_output_fibrils['real_length_nm']<=median_seeds_length
    seeds_only = colocalisation_output_fibrils[seeds]

    #count the number of fibrils & seeds
    total_number_ridges_detected=len(colocalisation_output_fibrils['Contour_ID'].value_counts())
    total_number_ridges_fibrils=len(fibrils_only['Contour_ID'].value_counts())
    total_number_ridges_seeds=len(seeds_only['Contour_ID'].value_counts())

    #find the number of colocalised hsp points in total
    num_colocalised_chaperones = len(colocalisation_output_fibrils.loc[colocalisation_output_fibrils['distance'] > -1])

    #find all points on fibrils that are colocalised with chaperones
    colocalised_fibrils = fibrils_only[fibrils_only['colocalisation']=='True']
    
    #get list of unique contour ID names 
    coloc_contour_IDs_list=colocalised_fibrils['Contour_ID'].to_list()

    #new df with all points on fibrils that ARENT colocalised with chaperones
    noncolocal_fibrils = fibrils_only[fibrils_only['colocalisation']=='False']
    #now filter our the noncolocalised points, that are colocalised somewhere else along the ridge by only taking those that have contour ID's NOT in the list of contour IDs associated with colocalised fibrils
    noncolocal_fibrils=noncolocal_fibrils[~noncolocal_fibrils['Contour_ID'].isin(coloc_contour_IDs_list)]


    #find the percentage of the total number of fibrils that are colocalised with a chaperone (only those unique contour ID's, not double ups from multiple chaperones bound)
    percentage_fibrils_colocalised=len(colocalised_fibrils['Contour_ID'].value_counts())/total_number_ridges_fibrils*100

    #make a new dataframe which matches up unique contour IDs (individual fibrils) with their fibril length by making a dictionary then converting to a new dataframe
    colocalised_fibrils_lengths = dict(zip(colocalised_fibrils.Contour_ID, colocalised_fibrils.real_length_nm))

    colocalised_fibrils_lengths=pd.DataFrame(colocalised_fibrils_lengths.items(), columns=['Contour_ID', 'real_length_nm']).T
    colocalised_fibrils_lengths=colocalised_fibrils_lengths.T

    colocalised_fibrils_lengths['colocalisation']='colocalised'

    colocalised_fibrils_lengths['protein_colocalised']=protein

    noncolocalised_fibrils_lengths = dict(zip(noncolocal_fibrils.Contour_ID, noncolocal_fibrils.real_length_nm))

    noncolocalised_fibrils_lengths=pd.DataFrame(noncolocalised_fibrils_lengths.items(), columns=['Contour_ID', 'real_length_nm']).T.T
    noncolocalised_fibrils_lengths['colocalisation']='non_colocalised'
    noncolocalised_fibrils_lengths['protein_colocalised']=protein

    lengths_data=pd.concat([colocalised_fibrils_lengths,noncolocalised_fibrils_lengths])

    lengths_data['treatment']=f'{treatment}'
    lengths_data['protein_colocalised']=protein
    
    percentage_fibrils_colocalised
    #save the lengths files of colocalised vs. non-colocalised fibrils
    lengths_data.to_csv(f'{output}lengths_data.csv')
    lengths_collated.append(lengths_data)
    

    return(lengths_collated, percentage_fibrils_colocalised)

def calculations_hsps(colocalisation_output_hsps):
#find the number of colocalised hsp points in total
    """calculate the amount of hsp colocalisations with each other (this is for 3 colour experiments).

    Args:
        colocalisation_output_hsps (df): all colocalised 488 labelled chaperones 

    Returns:
        int: % of fibril colocalisation
    """
    num_JB1_on_HSPA8 = len(colocalisation_output_hsps.loc[colocalisation_output_hsps['distance'] > -1])

    #new df with all points on fibrils that ARENT colocalised with chaperones
    noncolocal_chaps= colocalisation_output_hsps[colocalisation_output_hsps['colocalisation']=='False']


    #find the percentage of the total number of HSPA8 (AF647 labelled molecule) that are colocalised with a 488 labelled molecule
    total_number_hspa8=len(colocalisation_output_hsps['Slice'])
    percentage_coloc_hsp=(num_JB1_on_HSPA8/total_number_hspa8)*100
    
    return percentage_coloc_hsp
#get mean and median fibril lengths between proteins and colocalisation states
def create_length_summary(lengths_collated):
    """create a summary of the length of fibrils that are colocalised with chaperones, what their colocalisation is with each chaperone/each other

    Args:
        lengths_collated (df): df with all of the lengths of the fibrils

    Returns:
        df: with summary data
    """
    summary_data_fibrils=[]
    for treatment, df in lengths_collated.groupby('treatment'):
        treatment
        df
        for coloc, df2 in df.groupby('colocalisation'):
            coloc
            df2
            for protein, df3 in df2.groupby('protein_colocalised'):
                protein
                df3
                median_fibril_length=df3['real_length_nm'].median()
                mean_fibril_length=df3['real_length_nm'].mean()
                identifier=f'{treatment}_{protein}_{coloc}'
                summary_data_fibrils.append([identifier,median_fibril_length,mean_fibril_length])

    summary_data_fibrils=pd.DataFrame(summary_data_fibrils, columns=['treatment_proteins_colocalisation', 'median fibril length', 'mean fibril length'])
    summary_data_fibrils.to_csv(f'{output_folder}summary_data_fibril_length.csv')

    return summary_data_fibrils

def find_colocalisation(concatinate_coloc_data, find_coloc_fibs_filter, percentage_coloc_fibrils, treatment, output, median_seeds_length, fibril_hsp_files, zoom_factor):
    """calculate the percentage of fibrils colocalised, and save.

    Args:
        concatinate_coloc_data (func): function to concatinate the colocalisation files
        find_coloc_fibs_filter (func): function to find colocalised fibrils in order to act on them and do calculations
        percentage_coloc_fibrils (df): dataframe with the colocalised fibrils
        treatment (str): treatment being looped over
        output (str): where to save the treatment specific data
        median_seeds_length (int): length of seeds above whic hare fibrils
        fibril_hsp_files (list): list of fibrils + hsp files with colocalisation data
        zoom_factor (int): size of pixel in nm

    Returns:
        df
    """
    colocalisation_output_fibrils = concatinate_coloc_data(files=fibril_hsp_files)
            
    colocalisation_output_fibrils = pd.concat(colocalisation_output_fibrils)
    colocalisation_output_fibrils.rename(columns={'Contour ID':'Contour_ID'}, inplace=True)
            
    #correct the lengths column to account for pixel size in nm (zoomed in or out)
    colocalisation_output_fibrils['real_length_nm'] = colocalisation_output_fibrils['Length']*zoom_factor
    colocalised_fibrils = colocalisation_output_fibrils[colocalisation_output_fibrils['colocalisation']=='True']
    colocalised_fibrils.to_csv(f'{output}colocalised_fibrils_only.csv')

    percentage_colocalised = {}
    for protein, df in colocalisation_output_fibrils.groupby('protein_colocalised'):
        protein = protein
        lengths_collated, percentage_fibrils_colocalised = find_coloc_fibs_filter(colocalisation_output_fibrils=df, median_seeds_length=median_seeds_length, protein=protein)
        #add into dictionary for each protein the percentage of total fibrils that are colocalised with a chaperone
        percentage_colocalised[protein]=percentage_fibrils_colocalised

    percentage_colocalised = pd.DataFrame(percentage_colocalised.items(), columns=['proteins_colocalised', 'percent_colocalised'])
    percentage_colocalised['treatment'] = treatment

    percentage_colocalised.to_csv(f'{output}{treatment}_percentage_colocalisation.csv')
    percentage_coloc_fibrils.append(percentage_colocalised)
    return lengths_collated, percentage_coloc_fibrils

def concatinate(df, files_list):
    """concatinate files in a list

    Args:
        df (list): empty list to put files into
        files_list (list): list of files to read in and concat

    Returns:
        _type_: _description_
    """
    for fc in files_list:
        path=fc.replace('//', '/')
        path=path.replace('\\', '/')
        test=pd.read_csv(fc)
        cond=path.split('/')[5]
        test['treatment']=cond
        df.append(test)
    df=pd.concat(df)
    return df

def plot_colocalisation(output_folder, df, order_of_experiment, ylim, hue_var, title):
    """plot barplot of each protein colocalsied w fibril or each other

    Args:
        output_folder (str): where to save the plot
        df (df): dataframe w colocalisation and error
        order_of_experiment (list): list of the order of treatments in experiment, for plotting
        ylim (int  ): limit for plot
        hue_var (str): the variable / column to use for hue
        title (str): name of plot
    """
    ax=sns.barplot(data=df, x='treatment_to_plot', y='percent_colocalised', hue=hue_var, palette='Purples', order=order_of_experiment, alpha=0.45, edgecolor='black')
    ax.set_ylim(0,ylim)
    ax.set_xlabel('Protein colocalised with fibrils/each other')
    ax.set_ylabel('Percentage (%) colocalised')
    ax.set_title(f'{title}')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f'{output_folder}percent_proteins_colocalised.png')
    plt.show()

def plot_violin_2protein(output_folder, length, order_of_experiment):
    """plot the length data as violinplot (coloc and non coloc fibrils). 

    Args:
        output_folder (str): where to save plots
        length (df): df with all lengths of fibrils from the experiment
        order_of_experiment (list): order to plot
    """
    for protein, df in length.groupby('protein_colocalised'): 
        protein
        df
        fig, ax = plt.subplots()
        ax = sns.violinplot(x='treatment_to_plot', y= 'real_length_nm', hue= 'colocalisation', data=df, scale='width', palette='Blues_r', order=order_of_experiment, alpha=0.25)
        ax.legend(loc='upper center', ncol=2)
        plt.ylabel('Fibril length (nm)')
        plt.xlabel('Experiment condition')
        protein = protein.split('_')[0]
        plt.title(f'Length of fibrils colocalised with {protein} (in nm)')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(f'{output_folder}/{protein}_fibril_length_distributions.png')
        plt.show()

def workflow(get_experiment_info, concatinate_coloc_data, find_coloc_fibs_filter, create_length_summary, find_colocalisation, input_folder, output_folder, treatments, percentage_coloc_fibrils, hsp_only_paths, zoom_factor, median_seeds_length):
    """organise functions to run on treatments and save legnth and colocalisation data for each.

    Args:
        get_experiment_info (func): gathers information from folders that coloc files are contained within
        concatinate_coloc_data (func): concat files that have colocalisation data
        find_coloc_fibs_filter (func): does the work= finds fibrils and colocalisation / length data.
        create_length_summary (func): summary (mean, sd) to check data makes sense
        find_colocalisation (func): other work funciton: finds and calculates the % colocalised proteins
        input_folder (str): location of experiment data (from 'colocalisation' data saved from ImageJ output)
        output_folder (str): where to save the python output on colocalisation, and the plots
        treatments (list): the treatments that colocalisation data are nested within
        percentage_coloc_fibrils (df): dataframe with colocalisation data
        hsp_only_paths (list): paths with the colocalisation between hsps only (not fibrils)
        zoom_factor (int): this is SPECIFIC TO MICROSCOPE : this is how many nm per pixel
        median_seeds_length (int): average length of a seed (above which fibrils are considered fibrils rather than seeds. )

    Returns:
        dummary data about length and colocalisation
    """
    for treatment in treatments:
        #make a giant for-loop OR define functions here that we can then act on using multiple different input folders. these input folders can then all be acted on, and a graph can plot each of the 'length' columns in each of these. p.s add in seeds so that we can filter out the seeds that are in this batch
        treatment
        folder, concentration, experiment_number, output = get_experiment_info(input_folder, treatment)
    
        folder = f'{folder}fibrils-647/'

        fibril_hsp_files = [filename for filename in os.listdir(f'{folder}/') if 'fibril' in filename]

        #if there is a file in the list of hsp files go ahead here
        if len(fibril_hsp_files)>0:
            lengths_collated, percentage_coloc_fibrils = find_colocalisation(concatinate_coloc_data, find_coloc_fibs_filter, percentage_coloc_fibrils, treatment, output, median_seeds_length, fibril_hsp_files, zoom_factor)


        hsp_only_colocalised_files = [filename for filename in os.listdir(folder) if 'fibril' not in filename]
        for filename in hsp_only_colocalised_files:
            path = f'{folder}{filename}'
            hsp_only_paths.append(path)

    lengths_collated = pd.concat(lengths_collated)
    lengths_collated.to_csv(f'{output_folder}lengths_collated.csv')
    percentage_coloc_fibrils=pd.concat(percentage_coloc_fibrils)
    percentage_coloc_fibrils.to_csv(f'{output_folder}fibrils_only_coloc_output.csv')
    summary_data_fibrils = create_length_summary(lengths_collated)

    return percentage_coloc_fibrils, summary_data_fibrils

if __name__ == "__main__":   

    input_folder = 'data/Analysis_workflow/2_example_python_output/0_collect_data/Colocalisation/'
    output_folder = 'data/Analysis_workflow/2_example_python_output/2_colocalisation-length/Colocalisation/'  

    treatments = [treatment for treatment in os.listdir(input_folder)]
    lengths_collated = []
    percentage_coloc_fibrils = []
    hsp_only_paths = []
    #this says how many nm per pixel
    zoom_factor = 160
    #define the length of seeds batch, in nm
    median_seeds_length = 200

    #loop over each treatment under the top folder
    percentage_coloc_fibrils, summary_data_fibrils = workflow(get_experiment_info, concatinate_coloc_data, find_coloc_fibs_filter, create_length_summary, find_colocalisation, input_folder, output_folder, treatments, percentage_coloc_fibrils, hsp_only_paths, zoom_factor, median_seeds_length)

    #---------------------------------------------------------------------------------
    #PLOTTING
    #for plotting multiple treatments on the same graph from the % colocalisation dataframes saved before
    experiment_type  = [name for name in os.listdir(input_folder) ]

    concat = []
    length = []
    #gather files we generated above
    yikes_coloc = [[f'{root}/{name}' for name in files if 'percentage_colocalisation.csv' in name]for root, dirs, files in os.walk(f'{output_folder}/')]
    yikes_coloc = [item for sublist in yikes_coloc for item in sublist]

    #gather files we generated above
    yikes_len = [[f'{root}/{name}' for name in files if 'lengths_data.csv' in name]for root, dirs, files in os.walk(f'{output_folder}')]
    yikes_len = [item for sublist in yikes_len for item in sublist]

    #concatinate colocalisation files
    concat = concatinate(df=concat, files_list=yikes_coloc)
    length = concatinate(df=length, files_list=yikes_len)

    #dictionary to shorten/simplify the names for plotting treatments (adjust accordingly)
    dicto = {'2_example_python_output': 'Treat_plot',
            }

    #map dictionary so we have a new col for plotting
    concat['treatment_to_plot']=concat['treatment'].map(dicto)
    length['treatment_to_plot']=length['treatment'].map(dicto)

    #define order of plotting with the 'treatment to plot' name
    order_of_experiment = ['Treat_plot']

    #plot barplot of colocalisation
    plot_colocalisation(output_folder=output_folder, df=concat, order_of_experiment=order_of_experiment, ylim=100, hue_var='proteins_colocalised', title='HSPA8 colocalised with fibrils or DNAJB1')

    #violinplot of fibril length colocalised and not
    #plot distribution of fibril lengths coloc and not coloc
    plot_violin_2protein(output_folder, length, order_of_experiment)

