# A repository containing all code pertaining to chapter 5 of my PhD thesis (codes which generate the data in figures for chapter 5). These are detailed analyses of DNAJB1 and HSPA8 bound to fibrils of a-syn, imaged using TIRF microscopy.#

the scripts from this body of work are split into three categories within /analysis_scripts: 
0_ImageJscripts
1_analysis
2_thesis

# 0_ImageJscripts
These scripts do image processing on the RAW images coming from the microscope. 
They are annotated line by line, but a brief description of what each macro does is included below:

All scripts correct beam profile, find fibrils using a ridge detection plugin, find and save the XY coordinates of fibril (pixels) and Hsps. 

1_fibril_AF488_coloc_trajectories - finds fibrils and saves lengths, use this for experiments which have observed AF488-labelled proteins and fibrils labelled with ASCP dye (separate emission channels). Saves a file with the location and colocalisation state (coloc with a fibril pixel or not) of each AF488-labelled molecule. This also extracts and saves the trajectories (peak intensity in each frame of a stack) of the af488 molecules that ARE colocalised with a fibril. 

2_fibril_AF647_trajectories - This does the same as above ^ but for AF647 labelled proteins (i.e., in the same emission channel as the fibrils, but imaged at a different time using alternating excitation)

3_fibril_coloc_AF647_AF488_threecolour- This finds the colocalisation information as the above macros do, but this combines both of them. Finds fibrils, finds AF488 molecules, finds AF647 molecules. Finds colocalisation between AF647 molecules and fibrils, and af647 molecules and AF488 molecules, finds AF488 mols colocalised with fibrils. saves all of this. 

4_fibril_coloc_AF647_AF488_threecolour_trajectories - same as above, but also finds the peak intensity trajectories for the AF647 labelled molecules. 

each of these macros requires an ROI.zip, background images manually generated (gaussian blur of each channel generated in ImageJ), and a transformation matrices manually generated in ImageJ. These should be saved in the '647' folder of the raw data, where all of the images taken with the 647 laser on, reside. 

# 1_analysis
This repository contains the scripts that I used to analyse data that comes straight from image processing in ImageJ. So, for each experiment investigating fibril/Hsp interactions, these scripts were used to generate files which could then be collated for the 2_thesis scripts. 
To do for this:
- check annotation and clean up imports and functions
- simplify their names to actually what they do
- write a summary for each script detailing what it actually does.
- 
# 2_thesis
The data in the figures of my thesis is collated from many separate experiments. The scripts in this folder show the analysis that I did 
to do for these:
- delete unnecessary
- no hard coding: i.e., need to just show what I did with the collated, filtered file. use the collated filtered file as the 'input' in data, and that as the input directory in the actual scripts. 
- simplify to 'gathering data', 'analysing & plotting' and 'stats' for each figure (name them by figure.)
- simplify names of scripts
- write a simple 'what it does' for each script