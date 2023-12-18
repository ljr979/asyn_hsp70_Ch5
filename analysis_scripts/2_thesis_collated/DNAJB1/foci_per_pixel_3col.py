import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

input_top='imagejresults/FC3/Colocalisation/'
output_top='python_results/Colocalisation/FC3/'


treatments= [folder for folder in os.listdir(input_top)]


for exp in treatments: 
    input_folder=f'{input_top}{exp}/'
    output_folder=f'{output_top}{exp}/3_foci_per_pixel/'
    experiments=['fibrils-488','fibrils-647']

    alls=[]
    for c in experiments:
        collated=[]
        folder = f'{input_folder}{c}/'

        concentration = folder.split('/')[-3].split('_')[-1]

        experiment_number= folder.split('/')[-3].split('_')[0]
        
        output = f'{output_folder}{c}/'
        
        if not os.path.exists(output):
                os.makedirs(output)
        
        colocalised_files = [filename for filename in os.listdir(f'{folder}') if '.csv' in filename] 
        
        colocalised_fibrils = []
        for filepath in colocalised_files: 
            data = pd.read_csv(f'{folder}{filepath}')
            #data = data.drop(['Class'], axis=1)
            protein=filepath.split('_')[0]
            data['concentration']=concentration
            data['protein']=protein
            data['image_number']=filepath
            colocalised_fibrils.append(data)

        colocalised_fibrils=pd.concat(colocalised_fibrils)
        
        collated.append(colocalised_fibrils)
        collated=pd.concat(collated)

        collated.rename(columns={'Contour ID':'Contour_ID'}, inplace=True)
        #now we have a big dataframe with only the colocalisation data and we can get the number of foci on each fibril from this

        collated2=[]
        for Contour_ID, row in collated.iterrows():
            Contour_ID
            row
            row['ID_length_conc_image']=str(row['Contour_ID'])+'_'+str(row['Length'])+'_'+str(row['concentration'])+'_'+str(row['image_number'])
            row=pd.DataFrame(row).T
            collated2.append(row)

        collated2=pd.concat(collated2)
        coloc_collated2=collated2[collated2['distance']!=-1]
        coloc_collated2.to_csv(f'{output_folder}{c}/fibril_collated-colocal-only_data.csv')
        collated2.to_csv(f'{output_folder}{c}/HSPA8_fibril_collated-colocal_data.csv')

        #coloc_collated2.rename(columns={"ID_length_conc_image#": "ID_length_conc_image"}, inplace=True)
        #coloc_collated=collated[collated['distance']!=-1]

        proteins = {}
        lengths = {}

        # for conc, df in coloc_collated.groupby('concentration'):
        #     conc
        #     df
        #then make a new dataframe where I count the number of colocalised SPOTS on each contour, and then add in the matching fibril length to this counted number of foci
        counts = pd.DataFrame(coloc_collated2['ID_length_conc_image'].value_counts()).reset_index().rename(columns={'index':'ID_length_conc_image', 'ID_length_conc_image':'hsp_count'})
        contourIDs = [name for name in counts['ID_length_conc_image']]

        for contourID in contourIDs:
            length = coloc_collated2.loc[coloc_collated2['ID_length_conc_image'] == contourID, 'Length'].iloc[0]
            lengths.update({contourID:length})


        for contourID in contourIDs:
            protein = coloc_collated2.loc[coloc_collated2['ID_length_conc_image'] == contourID, 'protein'].iloc[0]
            proteins.update({contourID:protein})




        concentrations_dict = dict(zip(coloc_collated2.ID_length_conc_image, coloc_collated2.concentration))
        counts['lengths (pixel)'] = counts['ID_length_conc_image'].map(lengths)
        counts['concentration'] = counts['ID_length_conc_image'].map(concentrations_dict)
        counts['protein']= counts['ID_length_conc_image'].map(proteins)

        # with this new dataframe, divide the count of foci, by the length in pixels for each contour
        counts['foci_per_pixel'] = counts['hsp_count']/counts['lengths (pixel)']
        counts['log_foci']=np.exp(counts['foci_per_pixel'])
        alls.append(counts)

    alls=pd.concat(alls)
    #concentration_dict={'100pM':'0.1nM', '500pM':'0.5nM','1nM':'1nM', '5nM':'5nM', '10nM':'10nM', '50nM':'50nM', '100nM':'100nM'}
    #counts['concentration']=counts['concentration'].map(concentration_dict)


    #find the median of this?
    median_foci_per_pixel = counts['foci_per_pixel'].median()
    alls.to_csv(f'{output_folder}foci_per_length_unit.csv')
    ax = sns.violinplot(x='concentration', y='log_foci', data=alls, scale='width', palette='Reds', hue='protein')
    plt.title(f'{protein}foci per length unit fibril (pixel)')
    plt.ylabel(f'{protein}exp(foci per pixel)')
    plt.ylim(1.0,1.8)


    plt.savefig(f'{output_folder}HSPA8_foci_per_length_unit_colocalised.png')
    plt.show()




#for plotting multiple treatments on the same graph from the % colocalisation dataframes saved before
input_folder='python_results/Colocalisation/FC3/'
output_folder='python_results/Colocalisation/FC3/'
experiment_type=[name for name in os.listdir(input_folder) ]
experiment_type=[name for name in experiment_type if not 'png' in name ]
concat=[]
for name in experiment_type:

    filepath=f'{input_folder}{name}/3_foci_per_pixel/foci_per_length_unit.csv'
    test=pd.read_csv(filepath)
    test['treatment']=name
    concat.append(test)
    


concat=pd.concat(concat)

dicto={'A8-JB1-110-1nM-excess': 'A8+JB1+110 1nM(excess)',
       'A8-JB1-110-1nM-flowout':'A8+JB1+110 1nM(flowout)'
}

concat['treatment_to_plot']=concat['concentration'].map(dicto)

order_of_experiment= ['A8+JB1+110 1nM(excess)', 'A8+JB1+110 1nM(flowout)']

ax=sns.violinplot(data=concat, x='treatment_to_plot', y='log_foci', hue='protein', palette='Reds', order=order_of_experiment, alpha=0.45, edgecolor='black')
ax.set_ylim(1,1.8)
ax.set_xlabel('Protein colocalised with fibrils/each other')
ax.set_ylabel('exp(number of foci per pixel)')
ax.set_title('Number foci per pixel colocalised by A8/JB1')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f'{output_folder}FOCI_per_pixel_exponential.png')
plt.show()



