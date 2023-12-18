import os, re
import threading
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger


input_top = 'Y:/Chaperone_subgroup/LaurenR/20230821_Experiment100-b/20230821_Experiment100b-FC2-fibrils-A8-B1-SOD1-5nM/'

output_folder='python_results/'

if not os.path.exists(output_folder):
        os.makedirs(output_folder)

throwout=['ROI.zip','TransformationMatrices.txt', 'Channel488 and 647 dual Ex_Seq0000.nd2','Channel488 and 647 dual Ex_Seq0001.nd2', 'coords_added', 'Thumbs.db']

tops=[folder for folder in os.listdir(f'{input_top}') if folder not in throwout]

Experiment_number='Experiment100b-FC2-1'
protein='HSPA8'
treatment = f'A8-B1-SOD1-5nM'

def read_trajectories(hsp_traj_paths,fibril_coloc_paths):
    new_trajs=[]
    new_fibs=[]
    for item in hsp_traj_paths:
        item
        # trajectory_hsp=hsp_traj_paths[[item]]
        item=item.replace('//','/')
        item=item.replace('\\','/')
        hspnumber=item.split('/')[-1].split('.tif')[0].split('_')[-1]
        experiment_number1=item.split('/')[-4]
        for item1 in fibril_coloc_paths:
            item1
            item1=item1.replace('//','/')
            item1=item1.replace('\\','/')
            # fibril_coloc=fibril_coloc_paths[[item1]]
            fibril_number=item1.split('/')[-1].split('.tif')[0].split('_')[-1]
            experiment_number2=item1.split('/')[-4]
            if (hspnumber == fibril_number) and (experiment_number1==experiment_number2):


                matching_trajectory_file=pd.read_csv(item)
                matching_trajectory_file.drop([col for col in matching_trajectory_file.columns.tolist() if ' ' in col], axis=1, inplace=True)
                #matching_trajectory_file.drop([col for col in matching_trajectory_file.columns.tolist() if 'Mean_0' in col], axis=1, inplace=True)
                fibrils_coloc_file=pd.read_csv(item1)
                # client_trajectories.columns = [str(col) + '_client' for col in client_trajectories.columns]

                fibrils_coloc_file['hsp_coords_X_Y']=list(zip(fibrils_coloc_file.X2, fibrils_coloc_file.Y2))
                fibrils_coloc_file.rename(columns={'Contour ID':'Contour_ID'}, inplace=True)
                hspcoords=fibrils_coloc_file['hsp_coords_X_Y'].tolist()

                new_coords_IDs=[]
                for whoops in hspcoords:
                    XY_new=str(whoops[0])+'_'+str(whoops[1])
                    new_coords_IDs.append(XY_new)
                fibrils_coloc_file['new_coords_IDs_X_Y']=new_coords_IDs
                fibrils_coloc_file["Contour_ID"] = fibrils_coloc_file["Contour_ID"].values.astype(str)
                fibrils_coloc_file['new_ID_hspX_hspY'] = fibrils_coloc_file["Contour_ID"] +"_"+ fibrils_coloc_file["new_coords_IDs_X_Y"]


                

                fibril_IDS_coloc_df=fibrils_coloc_file[fibrils_coloc_file['distance']>-1]
                


                trans_traj_file = matching_trajectory_file.T.reset_index()




                #trans_traj_file=trans_traj_file.tail(trans_traj_file.shape[0] -1)
                #this little conditional loop is necessary because the MACRO 'FILTER FOR COLOCAL'  in ImageJ has a glitch (to be fixed). It picks up EITEHR the first colocalised point on a fibril, OR if the first colocalsied point is not the first identified pixel on that Contour ID (fibril) then it just always takes the SECOND ROW (position two on the first fibril) of the fibril/colocalised hsp dataframe. THis messes us up here because a lot of the time the second row is not actually colocalised (i.e. distance is '-1') and so the filters i put in to find only colocalised spots get rid of this first row, however the macro uses that table to find the trajectories, so in my trajectories file i have an extra (mostly rubbish) trajectory at the start. We can perhaps filter this out if it's rubbish later, but for now don't want to mess up the order of the other ones and want keep them all consistent. Instead in this loop we just add that second row back onto the dataframe, so that we can make sure the coordinates match the correct trajectories (can't add a column onto different length)

                #checks whether the second row is equal to -1 (not colocalised)
                if fibrils_coloc_file['distance'][1]==-1:
                    #if it is, grabs that first row as an index and a series
                    fkn_hell=fibrils_coloc_file.loc[1,:]
                    #makes a new dataframe where this series is popped onto the end of our colocalised fibrils dataframe (which had been filtered by distance <1)
                    annoying=fibril_IDS_coloc_df.append(fkn_hell)
                    #makes sure this extra row goes at the start
                    annoying=annoying.sort_index(ascending=True)
                    #defines the hspcoordinates of the colocalised hsps (plus that annoying one row that will be (0,0)
                    coloc_hsp_coords=[x for x in annoying['hsp_coords_X_Y'].tolist()]

                #coloc_hsp_coords=[x for x in annoying['hsp_coords_X_Y'].tolist() if not x == (0,0)]
                #now we make another condition so that if one of the first two rows of the dataframe has a colocalised hsp, we just proceed as normal (using filtered dataframe)
                if fibrils_coloc_file['distance'][0]>-1 or fibrils_coloc_file['distance'][1]>-1:
                    coloc_hsp_coords=[x for x in fibril_IDS_coloc_df['hsp_coords_X_Y'].tolist()]
                
                #append the coordinates of colocalised hsps onto the trajectories so that we can identify their location on the fibril
                trans_traj_file['hsp_coords_X_Y']=coloc_hsp_coords


                #creates a dictionary that contains the coordinates of the colocalised hsp as the key (in a tuple) and the new ID which is the contour ID of the fibril it's colocalised with, plus the XY coordinates of the hsp, separated by an underscore. 
                fibril_IDs_dict=dict(zip(fibril_IDS_coloc_df.hsp_coords_X_Y, fibril_IDS_coloc_df.new_ID_hspX_hspY))

                trans_traj_file['fibril_ID']=trans_traj_file['hsp_coords_X_Y'].map(fibril_IDs_dict)


                #makes a new column that is the new ID, which is made of the traj name, and then the fibril ID (including x and y coords of hsps)
                trans_traj_file['new_ID_hspX_hspY'] = trans_traj_file["index"] +"_"+ trans_traj_file["fibril_ID"]

                #sets this new trajectory name as the index
                trans_traj_file['index']=trans_traj_file['new_ID_hspX_hspY']
                #drops the columns we don't need now the trajectory name has changed
                test=['fibril_ID','new_coords_IDs_X_Y','new_ID_hspX_hspY','hsp_coords_X_Y']
                testo=trans_traj_file.drop(columns=[col for col in trans_traj_file if col in test])
                #fixes the shape up for feeding into the py4bleaching
                testo=testo.T
                testo.columns = testo.iloc[0] 
                testo=testo.tail(testo.shape[0] -1)

                new_trajs.append(testo)
                new_fibs.append(fibrils_coloc_file)


    
    new_fibs=pd.concat(new_fibs)
    #IMOPRTANT STEP is enumerating these new names, so that we know for sure they're not matching any by coincidence (number from 1>len df)
    new_fibs['new_ID_hspX_hspY']=[f'{new_ID_hspX_hspY}_{x}' for x, new_ID_hspX_hspY in enumerate(new_fibs['new_ID_hspX_hspY'])]

    new_trajs=pd.concat(new_trajs, axis=1)
    
    return new_fibs, new_trajs

def grab_trajectories_paths(input_folder, timepoint_folders):
    fibril_coloc_paths=[]
    hsp_traj_paths=[]
    for treatment in timepoint_folders:
        treatment
        yikes_trajectories=[[f'{root}/{name}' for name in files if '_colocal_traj.csv' in name]for root, dirs, files in os.walk(f'{input_folder}{treatment}/')]
        yikes_trajectories=[item for sublist in yikes_trajectories for item in sublist]
        yikes_trajectories=[item for item in yikes_trajectories if 'AF488' not in item]
        yikes_trajectories=[item for item in yikes_trajectories if 'af488' not in item]
        


        yikes_fibrils=[[f'{root}/{name}' for name in files if 'Colocalisation_analysis' in root]for root, dirs, files in os.walk(f'{input_folder}/')]
        yikes_fibrils=[item for sublist in yikes_fibrils for item in sublist]
        # yikes_fibrils=[item for item in yikes_fibrils if 'AF488' not in item]
        # yikes_fibrils=[item for item in yikes_fibrils if 'af488' not in item]
        yikes_fibrils=[item for item in yikes_fibrils if 'colocalized' not in item]
        yikes_fibrils=[item for item in yikes_fibrils if 'Trajectories' not in item]
        fibril_coloc_paths.append(yikes_fibrils)
        hsp_traj_paths.append(yikes_trajectories)
    return fibril_coloc_paths, hsp_traj_paths

def workflow(input_top, output_folder, tops, treatment):

    for conc in tops: 
        input_folder=f'{input_top}{conc}/'
        conc=conc.split('_')[-1]
        
        timepoint_folders=[folder for folder in os.listdir(f'{input_folder}') if folder not in throwout]
        fibril_coloc_paths, hsp_traj_paths=grab_trajectories_paths(input_folder, timepoint_folders)

        hsp_traj_paths=[item for sublist in hsp_traj_paths for item in sublist]
        fibril_coloc_paths=[item for sublist in fibril_coloc_paths for item in sublist]
        output = []
        for x in fibril_coloc_paths:
            if x not in output:
                output.append(x)
        fibril_coloc_paths=output      
        #add a line to filter for when it's a 3 colour experiment and we need to only look at fibrils colocalised with chaprones, not each other.
        fibril_coloc_paths=[item for item in fibril_coloc_paths if 'AF488_AF647_coloc' not in item]
        #and filter so that the JB1 molecules aren't included either as I didn't photobleach them so we don't need to get the coloc data included in any trajectory names :) 
        fibril_coloc_paths=[item for item in fibril_coloc_paths if 'fibrils_AF488_coloc' not in item]
        fibril_coloc_paths=[item for item in fibril_coloc_paths if 'non-coloc' not in item]

    new_fibs, new_trajs=read_trajectories(hsp_traj_paths,fibril_coloc_paths)



    output_folder_fibs=f'{output_folder}/fibril_coloc_data/{Experiment_number}_{protein}_{treatment}/'
    output_folder_trajectories=f'{output_folder}/Trajectories/{Experiment_number}_{protein}_{treatment}/'

    if not os.path.exists(output_folder_fibs):
            os.makedirs(output_folder_fibs)

    if not os.path.exists(output_folder_trajectories):
            os.makedirs(output_folder_trajectories)

    new_trajs.to_csv(f'{output_folder_trajectories}HSPA8_colocal_traj.csv')

    new_fibs.to_csv(f'{output_folder_fibs}/HSPA8_fibril_colocal_data.csv')


workflow(input_top, output_folder, tops, treatment)