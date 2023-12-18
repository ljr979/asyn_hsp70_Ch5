import os, re
import pandas as pd
import numpy as np
from loguru import logger
import glob
import shutil

#change these for experiment
input_folder = 'Y:/Chaperone_subgroup/LaurenR/20230821_Experiment100-b/'
output_folder = 'imagejresults/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

#can change based on proteins in experiment
hsp_colocalised = 'JB1_fibril_colocal_data'
fibril_colocalised= 'HSPA8_fibril_colocal_data'
hsps_coloc='JB1_HSPA8_colocal_data'

# jb1_colocalised = 'DNAJB1_colocal_traj'
#hspa8_colocal ='HSPA8_colocal_traj'
hspa8_noncoloc ='HSPA8_non-colocal_traj'
#change these for the experiment you're using
folders = {
#flow cell 1
#treatment 1
#coloc
#'20230713_Experiment98_FC1-1_HSPA8-JB1-ATP/647/Colocalisation_analysis/Trajectories/AF647/coloc/':('FC1-1/Trajectories/Experiment98-A8-B1-flowout/coloc/', hspa8_colocal),




#non-coloc
'20230821_Experiment100b-FC1-fibrils-A8-B1-hsp110-5nM/647/non-coloc_trajectories/647/':('FC1-1/Trajectories/Experiment100-FC1-1_HSPA8_A8-B1-110-5nM/non-coloc/', hspa8_noncoloc),

#coloc data
'20230821_Experiment100b-FC1-fibrils-A8-B1-hsp110-5nM/647/Colocalisation_analysis/':('FC1-1/Colocalisation/Experiment100-FC1-1_HSPA8_A8-B1-110-5nM/fibrils-647/', fibril_colocalised),


# #non-coloc
'20230821_Experiment100b-FC1-fibrils-A8-B1-hsp110-0.5uM/647/non-coloc_trajectories/647/':('FC1-2/Trajectories/Experiment100-FC1-2_HSPA8_A8-B1-110-0.5uM/non-coloc/', hspa8_noncoloc),

#coloc data
'20230821_Experiment100b-FC1-fibrils-A8-B1-hsp110-0.5uM/647/Colocalisation_analysis/':('FC1-2/Colocalisation/Experiment100-FC1-2_HSPA8_A8-B1-110-0.5uM/fibrils-647/', fibril_colocalised),




#treatment2 
#non-coloc
    '20230821_Experiment100b-FC2-fibrils-A8-B1-SOD1-5nM/647/non-coloc_trajectories/647/': ('Trajectories/non-coloc/Experiment100-FC2-1_HSPA8_A8-B1-SOD1/non-coloc/', hspa8_noncoloc),

    #coloc data
    '20230821_Experiment100b-FC2-fibrils-A8-B1-SOD1-5nM/647/Colocalisation_analysis/': ('FC1-2/Colocalisation/Experiment100-FC2-1_HSPA8_A8-B1-SOD1/fibrils-647/', fibril_colocalised),




    }



for old_folder, (new_folder, filetype) in folders.items():
    old_files = [filename for filename in os.listdir(f'{input_folder}{old_folder}') if not '_colocalized.csv' in filename]
    if not os.path.exists(f'{output_folder}{new_folder}'):
        os.makedirs(f'{output_folder}{new_folder}')
    for x, filename in enumerate(old_files): 
        shutil.copyfile(f'{input_folder}{old_folder}{filename}', f'{output_folder}{new_folder}{filetype}{x}.csv')


