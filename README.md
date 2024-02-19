# asyn_hsp70_Ch5 
A repository containing all code pertaining to chapter 5 of my PhD thesis (scripts which generate the data in all figures for chapter 5). These are detailed analyses of DNAJB1 and HSPA8 bound to fibrils of a-syn, imaged using TIRF microscopy

The scripts from this body of work are split into two categories within ```src```: 
```Analysis_workflow``` and ```Figures```. Analysis workflow includes all scripts used in the experimental workflow which generated the data which was eventually collated into the figures presented in Chapter 5 of the thesis. 

Unprocessed data, (and data from various stages of processing along the python workflow, as example outputs) as well as the data used to plot the figures in my thesis, are available in an online repository **(link here)** and can be downloaded into this workspace to explore the analysis workflow and plot the figures. To do so, use the download script in the ```utilities``` folder and follow the data directory that is downloaded with the repository. 

All paths within analysis scripts direct to the appropriate location within the ```data``` folder in order to run the script and produce the desired output. 

# *Analysis_workflow*
## 0_ImageJscripts

These scripts perform image processing on the RAW images coming from the microscope. 
They are annotated line by line, but a brief description of what each macro does is included below:

All scripts correct beam profile, find fibrils using a ridge detection plugin, find and save the XY coordinates of fibril (pixels) and Hsps. 

```1_fibril_AF488_coloc_trajectories``` - finds fibrils and saves lengths, use this for experiments which have observed AF488-labelled proteins and fibrils labelled with ASCP dye (separate emission channels). Saves a file with the location and colocalisation state (coloc with a fibril pixel or not) of each AF488-labelled molecule. This also extracts and saves the trajectories (peak intensity in each frame of a stack) of the af488 molecules that ARE colocalised with a fibril. 

```2_fibril_AF647_trajectories``` - This does the same as above ^ but for AF647 labelled proteins (i.e., in the same emission channel as the fibrils, but imaged at a different time using alternating excitation)

```3_fibril_coloc_AF647_AF488_threecolour```- This finds the colocalisation information as the above macros do, but this combines both of them. Finds fibrils, finds AF488 molecules, finds AF647 molecules. Finds colocalisation between AF647 molecules and fibrils, and af647 molecules and AF488 molecules, finds AF488 mols colocalised with fibrils. saves all of this. 

```4_fibril_coloc_AF647_AF488_threecolour_trajectories``` - same as above, but also finds the peak intensity trajectories for the AF647 labelled molecules. 

Each of these macros requires an ```ROI.zip```, background images manually generated (gaussian blur of each channel generated in ImageJ), and a ```transformation matrices``` manually generated in ImageJ. These should be saved in the '```647```' folder of the raw data, where all of the images taken with the 647 laser on, reside. 

## 1_python_workflow
This repository contains the scripts that I used to analyse data that comes straight from image processing in ImageJ. So, for each experiment investigating fibril/Hsp interactions, these scripts were used to generate files which could then be collated for the 2_thesis scripts. A brief description is included below to use as a directory. 

```0_add_co-ords``` 
- This is to be used when you plan to find the location of the fluorescent molecule binding the fibril. A script to add the coordinate of the molecule into the NAME of the trajectory before it is run through py4b to determine number of subunits. 

```0_collect_data``` 
- Use this to transfer all colocalisation and trajectory data out of where it is saved from imageJ, into the python repo for analysis. Renames everything so that it is compatible with latter scripts, and means you don't have to worry about what things are labelled prior to this point. but you can skip this one if you're organised and just have them named correctly and nested appropriately earlier on.

```1_hsp_stoichiometry```:
-  This script runs py4bleaching on trajectories and outputs 'molecule_counts.csv' which is the count of the number of subunits per trajectory/molecule. 

```2_coloc_length_```: 
- The 'colocalisation data' goes into this script, which has the length of each fibril, and whether it is colocalised at any point with another molecule. This outputs the length of fibrils that are colocalised and NOT colocalised at all with a hsp. It also outputs the percentage of fibrils that re colocalised at least in one location on the fibril.

```3_subunits_per_foci``` : 
- This just organises and plots the subunits per molecule data that has been output into molecule_counts.csv from py4bleaching.

```4_density_pixel``` : 
- Takes each fibril, counts the number of colocalised spots on the fibril, and divides this count by the number of pixels (length) of the fibril. as such, it is reporting the density of fluorescent spots on the fibril, per pixel length (this accounts for the differences in length of different fibrils). 

```5_end_vs_middle_init``` : 
- The first step in determining the location of binding of hsp on fibrils. this script finds the end region and middle region of each fibril, and then calculates the distance of each colocalised hsp from the end of that fibril it's bound to. then it assigns each molecule with 'END' or 'MIDDLE' and saves this data. 

```6_end_or_middle_plot``` : 
- Filters for fibrils which are long enough to (most) accurately determine the end region from the middle. accounts for the 'zoom' of the microscope (i.e, allows us to report binding per um for each region RATHER THAN PIXEL, which can't be combined between experiments if they're still in 'pixel' form). Plots as a bar plot the density/um in each region. 

```7_subunits_end_middle``` : 
- This uses the counts obtained via py4bleaching and combines it with the end_vs_middle output. assigns the name of the trajectory as the same name given to the colocalised spots on the fibrils, and matches them  with their location. i.e., we get out the subunits per molecule, and whether this is located at the end or middle. plots this as a boxplot per region, per treatment.

## utilities
```py4bleaching```:
- See ['add in link to py4bleahcing repository] to see the installable package. However, the package included here is for reference. It analyses trajectories extracted from images, and finds the number of fluorescently labelled subunits per Hsp molecule. This is saved in the py4bleaching output in the data repository.  

# *Figures*
- The data in the figures of my thesis is collated from many separate experiments. The scripts in this folder show the analysis that was performed on the collated output from the Analysis_workflow scripts. 

- Included in here are the scripts that perform all statistical analysis for each of the sections of this thesis chapter.

- By running each of the scripts within these folders, the corresponding figures can be regenerated and verified. 







