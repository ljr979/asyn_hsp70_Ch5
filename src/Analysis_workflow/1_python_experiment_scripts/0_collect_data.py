"""This script simply tranfers trajectories and colocalisation data from the experiment folder (where imageJ has output it) and sets up the folders and filenames so that they are identical for all experiments and the scripts that follow need less adjustment for each experiment.
"""
import os, re
import pandas as pd
import numpy as np
from loguru import logger
import glob
import shutil

def move_files(input_folder, output_folder, folders):
    """function to move and rename colocalisation and trajectory files from where they were output from image processing, into the repository with nomenclature that later scripts call on to carry treatment information etc

    Args:
        input_folder (str): top location of coloc and trajectory data
        output_folder (str): location to output the files
        folders (dict): dictionary with the folder the files exist in now, versus the folder nesting that you want them to be in in new location
    """
    for old_folder, (new_folder, filetype) in folders.items():
        old_files = [filename for filename in os.listdir(f'{input_folder}{old_folder}') if not '_colocalized.csv' in filename]
        if not os.path.exists(f'{output_folder}{new_folder}'):
            os.makedirs(f'{output_folder}{new_folder}')
        for x, filename in enumerate(old_files): 
            shutil.copyfile(f'{input_folder}{old_folder}{filename}', f'{output_folder}{new_folder}{filetype}{x}.csv')

if __name__ == "__main__":
        
    #change these for experiment
    #top level experiment folder path
    input_folder = '/'
    output_folder = 'imagejresults/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    #this is what they will be called.
    #can change based on proteins in experiment
    hsp_colocalised = 'JB1_fibril_colocal_data'
    fibril_colocalised = 'HSPA8_fibril_colocal_data'
    hsps_coloc = 'JB1_HSPA8_colocal_data'

    #jb1_colocalised = 'DNAJB1_colocal_traj'
    #hspa8_colocal ='HSPA8_colocal_traj'
    hspa8_noncoloc = 'HSPA8_non-colocal_traj'
    #change these for the experiment you're using
    folders = {

    #non-coloc
    'file_location_imagej_output/647/non-coloc_trajectories/647/':('treatment1/Trajectories/Experiment_number-treatment1_protein_treatment2-concentration/non-coloc/', hspa8_noncoloc),

    #coloc data
    'file_location_imagej_output/647/Colocalisation_analysis/':('treatment1/Colocalisation/Experiment_number-treatment1_protein_treatment2-concentration/fibrils-647/', fibril_colocalised),

        }


    move_files(input_folder, output_folder, folders)


