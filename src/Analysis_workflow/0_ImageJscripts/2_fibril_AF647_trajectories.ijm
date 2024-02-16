macro "Alignment macro [h]" {
//this script looks at FIBRIL AND AF647 colocalisation, and gets the 647 trajectories. so
//have to make sure to first select fibril, then select AF647 image.
//FIBRIL IMAGE SELECTED FIRST
var path = File.openDialog ("RAW image");
var dir = File.getParent(path);
var name = File.getName(path);

splitter = split(dir, File.separator) //splits image(s) of directory  (3) and separates images within folder
Array.print(splitter)//puts all things within array (the array defined by splitter) in a line
len = splitter.length // length of the array of splitter
trim = len -1 //needs to be 7 but we have 8 so uses only 7 of them (8-1)
ROIarr = Array.slice(splitter,0,trim) //the array of ROI path (our channel ROI'S) found from split images of directory defined in 'splitter' and only the 'trim' version (7 bits). not sure about the 0
Array.print(ROIarr)//puts all the ^^^^ into a line
//str = "My Passport/HONOURS/TIRF analysis/190711/"; //makes this a string to be added to other strings
str = ""
filesep = File.separator;
     for (i=0; i<ROIarr.length-1; i++) 
         str = str + ROIarr[i] + filesep; //literally no idea what's going on from 14-17
     str = str + ROIarr[ROIarr.length-1] + filesep; 
ROI = filesep + str + "ROI.zip"; //the path for the script to follow to find its' way to the ROIs
print("Raw Image path= " + path)//writes the raw image path and then the path that it follows ((3)
print("ROI path= " + ROI); //see above, but for ROI path (18)
output = dir + filesep 
print("Output path= " + output)
transformationfile = filesep + str + "TransformationMatrices.txt";
Hsp_results = output + "Hsp_results.csv";
fibril_results = output + "Fibril_results.csv";
Hsp_image = output + "Hsp.tif";

//SELECT AF647 IMAGE NEXT
//read in the image with the AF647 trajectories within them
var AF647_imagepath = File.openDialog ("AF647 image");
var dir = File.getParent(AF647_imagepath);
var name = File.getName(AF647_imagepath);

splitter = split(dir, File.separator) //splits image(s) of directory  (3) and separates images within folder
Array.print(splitter)//puts all things within array (the array defined by splitter) in a line
len = splitter.length // length of the array of splitter
trim = len -1 //needs to be 7 but we have 8 so uses only 7 of them (8-1)
ROIarr = Array.slice(splitter,0,trim) //the array of ROI path (our channel ROI'S) found from split images of directory defined in 'splitter' and only the 'trim' version (7 bits). not sure about the 0
Array.print(ROIarr)//puts all the ^^^^ into a line
//str = "My Passport/HONOURS/TIRF analysis/190711/"; //makes this a string to be added to other strings
str = ""
AF647_filesep = File.separator;
     for (i=0; i<ROIarr.length-1; i++) 
         str = str + ROIarr[i] + filesep; //literally no idea what's going on from 14-17
     str = str + ROIarr[ROIarr.length-1] + filesep; 
ROI = filesep + str + "ROI.zip"; //the path for the script to follow to find its' way to the ROIs
print("af647 Image path= " + AF647_imagepath)//writes the raw image path and then the path that it follows ((3)
print("ROI path= " + ROI); //see above, but for ROI path (18)
output_647 = dir + filesep 
print("Output path= " + output_647)
run("Set Measurements...", "centroid stack redirect=None decimal=3");
transformationfile = filesep + str + "TransformationMatrices.txt";

AF647_results = output_647 + "AF647_results.csv";

//define analysis images for both (i.e., corrected versions)
AF647_image = output_647 + "MAX_Background_corrected_AF647.tif";
bg_AF647= output_647+ "Background_corrected_AF647.tif"


run("Set Measurements...", "centroid stack redirect=None decimal=3");

//define background images for each
Hsp_background_image = filesep + str + "AF647_background.tif"
fib_background_image= filesep + str + "fibril_background.tif"

fibril = output + "fibril_MAX.tif";

//create output folder
Fibcolocalpath = filesep + str + "/Colocalisation_analysis/";
print(Fibcolocalpath);
File.makeDirectory(Fibcolocalpath);
//------------------------------------------------------
//now start actual analysis here
// Fibril channel
open(path) //select raw image
Imagetitle =getTitle()
rename("Raw"); 
//open ROI.zip where you defined channels and saved the file.
roiManager("Open", ROI);
roiManager("Select", 0);	//this will be the fibril channel
run("Duplicate...", "title=fibril duplicate");

//Background correction 
open(fib_background_image);
rename("Background_fib");
selectWindow("fibril");
run("32-bit");
run("Beam Profile Correction", "electronic_offset=500 background_image=Background_fib");
setMinAndMax(0, 65536);
run("16-bit");
//because these two are the same channel, the crop to fix the edges is called 1, not 2
roiManager("Select", 1);
run("Duplicate...", "title=fibril_background_corrected duplicate");
saveAs("tif", output + "Background_corrected_fib.tif");
selectWindow("Background_fib")
close();

//this will save the max fibril intensity to find fibril ridges
selectWindow("Background_corrected_fib.tif")
run("Z Project...", "stop=20 projection=[Max Intensity]");
run("Properties...", "channels=1 slices=1 frames=1 unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000");
saveAs("tif", output + "fibril_MAX.tif");
close()
roiManager("reset");
selectWindow("Background_corrected_fib.tif");
close()


//Hsp channel- don't need to align this because they are in the same channel but separate images
//(i.e. alternate excitation)
open(AF647_imagepath) //select raw image
Imagetitle =getTitle()
Imagetitle =getTitle()
rename("AF647_image"); 

roiManager("Open", ROI);
roiManager("Select", 0);	//this will be the AF647 channel
run("Duplicate...", "title=AF647 duplicate");

//Background correction 
open(Hsp_background_image);
rename("Background_AF647");
selectWindow("AF647");
run("32-bit");
run("Beam Profile Correction", "electronic_offset=350 background_image=Background_AF647");
setMinAndMax(0, 65536);
run("16-bit");
//crop edges again
roiManager("Select", 1);
run("Duplicate...", "title=AF647_background_corrected duplicate");
saveAs("tif", output_647 + "Background_corrected_AF647.tif");
selectWindow("Background_AF647")
close();

//creat max intensity file for peaks to be identified in AF647 channel
list = getFileList(dir); 
selectWindow("Background_corrected_AF647.tif");
run("Z Project...", "stop=20 projection=[Max Intensity] ");
run("Properties...", "channels=1 slices=1 frames=1 unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000");
saveAs("tiff", output_647+getTitle()); 
selectWindow("MAX_Background_corrected_AF647.tif");
close()
selectWindow("Background_corrected_AF647.tif");
close()

roiManager("reset");

//now do ridge detection of fibril image
open(fibril);
run("Enhance Contrast", "saturated=0.35");
run("8-bit");
run("Ridge Detection", "line_width=2.5 high_contrast=230 low_contrast=87 extend_line displayresults add_to_manager method_for_overlap_resolution=SLOPE sigma=1.2 lower_threshold=2 upper_threshold=11 minimum_line_length=4 maximum=0");
waitForUser("Delete the junctions from the ROI manager before pressing OK")
selectWindow("fibril_MAX.tif");
close()
selectWindow("Summary");
//save length data and x y coordinates of peaks
saveAs("Results", output + "Length_data.csv");
run("Close");
selectWindow("Results");
saveAs("Results", output + "Fibril_results.csv");
run("Close");
selectWindow("Junctions");
run("Close");
print("\\Clear")
roiManager("reset");


//Now do AF647 MOLECULE analysis 
open(AF647_image);
rename("Hsp");
//Find the peaks within the Hsp image
run("Enhance Contrast", "saturated=0.35");
run("Clear Results")
run("Peak Finder");// need to check the threshold here
waitForUser("Have you selected your peaks? \n \nPress OK to continue");
roiManager("Measure");	//this will allow for a table to be produced
saveAs("Results", output_647 + "Hsp_results.csv"); 
close();
selectWindow("Results"); 
     run("Close"); 
roiManager("reset");


//colocalization analysis
open(fibril_results)
IJ.renameResults("Fibril");
open("Hsp_results.csv")
IJ.renameResults("Hsp");
run("count colocalized peaks Fibril", "table_1=Fibril table_2=Hsp maximum_distance=3");
selectWindow("Colocalization");
saveAs("tiff", output + "Colocalization.tiff");
close()
selectWindow("Fibril");
saveAs("Results", output + "Fibril_colocalisation.csv"); 
saveAs("Results", Fibcolocalpath + Imagetitle + ".csv");
IJ.renameResults("Results");
selectWindow("Hsp");
		run("Close")
//now filter these data for only those colocalised not all the spots combined.
runMacro("C:/Users/ljr979/Desktop/Fiji.app/macros/Filter_for_colocal.ijm")
IJ.renameResults("Fibril_colocalisation_only")
selectWindow("Fibril_colocalisation_only");
saveAs("Results", output + "Fibril_colocalisation_only.csv"); 
saveAs("Results", Fibcolocalpath + Imagetitle + "_colocalized.csv");
run("Close");
roiManager("reset");

//create output folder for trajectories that are colocalised (only)
coloc_trajectory_output = filesep + str + "/coloc_trajectories/"
print(coloc_trajectory_output);
File.makeDirectory(coloc_trajectory_output);

//now find the coordinates of the hsps colocalised with fibrils, and get the trajectories at these coordinates
open(bg_AF647)
open(output + "Fibril_colocalisation_only.csv");
Table.rename("Fibril_colocalisation_only.csv", "Results");
runMacro("C:/Users/ljr979/Desktop/Fiji.app/plugins/Macros/v2_results_to_ROI_X2.txt");
selectWindow("Results");
        run("Close");
runMacro("C:/Users/ljr979/Desktop/Fiji.app/plugins/Macros/peak_intensity.txt"); // you need to put the location of your script
selectWindow("Results");
        saveAs("Results", coloc_trajectory_output + Imagetitle + "_HSP_colocal_traj.csv");
        run("Close");;
roiManager("reset");

selectWindow("Raw");
		run("Close");
selectWindow("fibril")
		run("Close")

selectWindow("AF647")
		run("Close")
selectWindow("AF647_image")
		run("Close")	

//NOW FIND NON COLCOALISED SPOTS TRAJECTORIES

//filter Non-colocal tables
noncoloc_output_folder = filesep + str  + "/non-coloc_trajectories/"
print(noncoloc_output_folder);
File.makeDirectory(noncoloc_output_folder);

hsp_noncoloc_trajectories = noncoloc_output_folder + "/647/"
print(hsp_noncoloc_trajectories);
File.makeDirectory(hsp_noncoloc_trajectories);

//Hsp
open(fibril_results);
IJ.renameResults("fibril");
open(output_647 + "Hsp_results.csv");
IJ.renameResults("Hsp");
run("count colocalized peaks Andrew", "table_1=Hsp table_2=fibril maximum_distance=3");

selectWindow("fibril");
		run("Close")
selectWindow("Hsp");
IJ.renameResults("Results");
runMacro("C:/Users/ljr979/Desktop/Fiji.app/macros/Filter_for_noncolocal.ijm")
IJ.renameResults("Hsp_non-colocalisation")
selectWindow("Hsp_non-colocalisation")
saveAs("Results", output + "647_non-colocalisation.csv"); 
selectWindow("647_non-colocalisation.csv");
     run("Close"); 

open(bg_AF647);	
open(output + "647_non-colocalisation.csv");
IJ.renameResults("Results")
runMacro("C:/Users/ljr979/Desktop/Fiji.app/macros/results_to_ROI.txt"); // you need to put the location of your script
selectWindow("Results");
        run("Close");
runMacro("C:/Users/ljr979/Desktop/Fiji.app/macros/peak_intensity.txt"); // you need to put the location of your script
//selectWindow("Results");
        saveAs("Results", output + "Hsp_non-colocal_traj.csv");
selectWindow("Results");
     	saveAs("Results", hsp_noncoloc_trajectories + Imagetitle + "_Hsp_non-colocal_traj.csv");
        run("Close");
selectWindow("Background_corrected_AF647.tif");
        run("Close");
selectWindow("Background_corrected_AF647-1.tif");
        run("Close");
roiManager("reset");
		
}


