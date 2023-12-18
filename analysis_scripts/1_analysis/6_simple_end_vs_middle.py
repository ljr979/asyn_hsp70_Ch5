
from enum import unique
import os, re
import threading
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
import math

input_folder='python_results/Colocalisation/'
output_folder='python_results/end_vs_middle/step_one/'

treatments=[tm for tm in os.listdir(input_folder)]


if not os.path.exists(output_folder):
        os.makedirs(output_folder)

filepaths=[[f'{root}/{name}' for name in files if 'collated-colocal_data' in name]for root, dirs, files in os.walk(f'{input_folder}')]
#filepaths=[[f'{root}{name}' for name in files if 'fibrils-647HSPA8_fibril_collated-colocal_data.csv' in name]for root, dirs, files in os.walk(f'{input_folder}')]
filepaths=[item for sublist in filepaths for item in sublist if not 'coords_added' in item]

#these two lines FILTER for 488 OR 647 LABELLED CHAPERONE. so comment out whichever one you are looking at and change defined protein accordingly (need to do this in a more sophisticated way when I have time)

JB1_filepaths=[item for item in filepaths if 'fibrils-488' in item]
filepaths=[item for item in filepaths if 'fibrils-647' in item]

#filter for colocalised chaperone on fibril 
def filter_colocal_fix_glitch(fibrils):
    fibril_IDS_coloc_df=fibrils[fibrils['distance']>-1]

    #checks whether the second row is equal to -1 (not colocalised)
    if fibrils['distance'][1]==-1:
        #if it is, grabs that first row as an index and a series
        fkn_hell=fibrils.loc[1,:]
        #makes a new dataframe where this series is popped onto the end of our colocalised fibrils dataframe (which had been filtered by distance <1)
        annoying=fibril_IDS_coloc_df.append(fkn_hell)
        #makes sure this extra row goes at the start
        annoying=annoying.sort_index(ascending=True)
        #defines the hspcoordinates of the colocalised hsps (plus that annoying one row that will be (0,0)
        fibril_IDS_coloc_df=annoying

    return fibril_IDS_coloc_df, fibrils
        
#loop through colocalised dataframe and group by Contour ID and length, THEN we use this function to get the ENDPOINT (X,Y) of both ends of the fibrils, and map this onto the colocalised fibrils dataframe (now new_fibrils_df)

def find_fibril_ends(fibril_IDS_coloc_df, fibrils, IDS_list):
    new_fibrils_dict={}
    #find unique fibril names in this colocalised dataframe
    uniques=fibril_IDS_coloc_df['new_ID_hspX_hspY_num'].tolist()

    # IDS_list=list(fibril_IDS_coloc_df['Contour_ID'])   
    #loop over the fibrils (all points on a fibril) dataframe and group into individual fibrils (based on their contour ID, gives ALL points on the fibril)
    for IDs, df in fibrils.groupby('ID_length_conc_image'):
        IDs
        df
    # if IDs in IDS_list:
        #matched_coloc_df=fibril_IDS_coloc_df[fibril_IDS_coloc_df['new_ID_hspX_hspY']==df['new_ID_hspX_hspY']]
        #define the end point as the first and last row for x and y for each fibril
        EndX1=list(df['X'])[0]
        EndX2=list(df['X'])[-1]
        EndY1=list(df['Y'])[0]
        EndY2=list(df['Y'])[-1]
        #keep these in a list for each fibril
        dit_list=[EndX1, EndX2, EndY1, EndY2]

        #now loop over the column that contains their unique name assigned earlier with the coordinates in it, and if any row on the fibril contains a number that is also existing in the list of unique names, put all the ends data we just got out, into a dictionary with the unique fibril name as these are the data we want to look at later
        for x in df['new_ID_hspX_hspY_num']:
            x

            if x in uniques:
                
                newpair={x:dit_list}
                new_fibrils_dict.update(newpair)
            # else: 
            #     newpair={x:'None'}
            #     new_fibrils_dict.update(newpair)
    

    #map this dictionary onto our colocalised df, so that we now know for each colocalised spot on the fibril, where the start and end of that fibril is
    fibril_IDS_coloc_df['EndX1, EndX2, EndY1, EndY2']=fibril_IDS_coloc_df['new_ID_hspX_hspY_num'].map(new_fibrils_dict)

    new_fibrils_df=pd.concat([fibril_IDS_coloc_df, fibril_IDS_coloc_df['EndX1, EndX2, EndY1, EndY2'].apply(pd.Series)], axis = 1)
    new_fibrils_df.rename(columns={0:'EndX1',1:'EndX2',2:'EndY1',3:'EndY2', 'X2':'Point_X', 'Y2':'Point_Y'}, inplace=True)
    return new_fibrils_df, fibril_IDS_coloc_df


def end_or_middle(new_fibrils_df, radius):
    updated_distance_fibs=[]
    for row , df in new_fibrils_df.groupby('ID_length_conc_image'):
        row
        df
        EndX1=df['EndX1'].values[0]
        EndX2=df['EndX2'].values[0]
        EndY1=df['EndY1'].values[0]
        EndY2=df['EndY2'].values[0]
        Point_X=df['Point_X'].values[0]
        Point_Y=df['Point_Y'].values[0]
        dist_1=math.sqrt((math.pow((EndX1-Point_X), 2) + math.pow((EndY1-Point_Y), 2)))
        dist_2=math.sqrt((math.pow((EndX2-Point_X), 2) + math.pow((EndY2-Point_Y), 2)))
        df['dist_1']=dist_1
        df['dist_2']=dist_2
        if dist_1 < radius or dist_2 < radius:
            df['Where']='END'
        else:
            df['Where']='MIDDLE'

        updated_distance_fibs.append(df)
    updated_distance_fibs=pd.concat(updated_distance_fibs)
    return updated_distance_fibs

def workflow(filepath, output_folder, Experiment_number):

#first loop reads in all coloc data and defines stuff about the experiment
    # new_filename=f'{flow}_'+ f'{ATP}_' + '-'.join(proteins)
    fibrils=pd.read_csv(f'{filepath}')
    fibrils['new_ID_hspX_hspY_num'] = fibrils["Contour_ID"].astype(str) + "_"+ fibrils["X2"].astype(str)+ "_"+ fibrils["Y2"].astype(str)
    fibrils['new_ID_hspX_hspY_num']=[f'{new_ID_hspX_hspY}_{x}' for x, new_ID_hspX_hspY in enumerate(fibrils['new_ID_hspX_hspY_num'])]


    fibril_IDS_coloc_df, fibrils=filter_colocal_fix_glitch(fibrils)
    fibril_IDS_coloc_df.rename(columns={'Contour ID':'Contour_ID'}, inplace=True)
    fibrils.rename(columns={'Contour ID':'Contour_ID'}, inplace=True)

    #make a list of the ID names that exist in this dataframe which is just colocalised spots
    IDS_list=list(fibril_IDS_coloc_df['ID_length_conc_image'])   

    #compare these dataframes and define the end region of each fibril
    new_fibrils_df, fibril_IDS_coloc_df=find_fibril_ends(fibril_IDS_coloc_df, fibrils, IDS_list)
    filepath=filepath.replace('\\', '/')
    treatment=filepath.split('/')[2]
    #nfow we have to classify how far the distance we want to classify as an 'end' is
    #define number of pixels that classifies an end
    radius = 2
    updated_distance_fibs=end_or_middle(new_fibrils_df=new_fibrils_df, radius=radius)
    updated_distance_fibs.to_csv(f'{output_folder}{Experiment_number}_{treatment}_end_vs_middle.csv')


Experiment_number='Experiment_100-SODcontrol-hspa8'
for filepath in filepaths:
    filepath
    workflow(filepath, output_folder, Experiment_number)


