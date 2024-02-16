macro "finding colocalisations with AF488-protein, AF647-protein, and fibrils" {
//setup for fibril and AF488 labelled protein in the same image

//FIBRIL/AF488 CONTAINING IMAGE FIRST
var path = File.openDialog ("RAW image");
var dir = File.getParent(path);
var name = File.getName(path);

splitter = split(dir, File.separator) //splits image(s) of directory  (3) and separates images within folder
Array.print(splitter)//puts all things within array (the array defined by splitter) in a line
len = splitter.length // length of the array of splitter
trim = len -1 //needs to be 7 but we have 8 so uses only 7 of them (8-1)
ROIarr = Array.slice(splitter,0,trim) //the array of ROI path (our channel ROI'S) found from split images of directory defined in 'splitter' and only the 'trim' version (7 bits). not sure about the 0
Array.print(ROIarr)//puts all the ^^^^ into a line

str = ""
filesep = File.separator;
     for (i=0; i<ROIarr.length-1; i++) 
         str = str + ROIarr[i] + filesep; //literally no idea what's going on from 14-17
     str = str + ROIarr[ROIarr.length-1] + filesep; 

ROI = filesep + str + "ROI.zip"; //the path for the script to follow to find its' way to the ROIs
print("Raw Image path= " + path)//writes the raw image path and then the path that it follows ((3)
print("ROI path= " + ROI); //see above, but for ROI path (18)
fib_488_output = dir + filesep 
print("Output path= " + fib_488_output)
transformationfile = filesep + str + "TransformationMatrices.txt";
Hsp_results = fib_488_output + "Hsp_results.csv";

Hsp_image = fib_488_output + "Hsp.tif";

run("Set Measurements...", "centroid stack redirect=None decimal=3");

//setup for AF647 labelled protein which has been imaged separately and shows up in the same channel as the FIBRILS but in a separate image
//HERE SELECT THE RAW IMAGE CONTAINING AF647 LABELLED MOLECULES
var AF647_imagepath = File.openDialog ("AF647 image");
var dir = File.getParent(AF647_imagepath);
var name = File.getName(AF647_imagepath);

splitter = split(dir, File.separator) //splits image(s) of directory  (3) and separates images within folder
Array.print(splitter)//puts all things within array (the array defined by splitter) in a line
len = splitter.length // length of the array of splitter
trim = len -1 //needs to be 7 but we have 8 so uses only 7 of them (8-1)
ROIarr = Array.slice(splitter,0,trim) //the array of ROI path (our channel ROI'S) found from split images of directory defined in 'splitter' and only the 'trim' version (7 bits). not sure about the 0
Array.print(ROIarr)//puts all the ^^^^ into a line
str = ""
AF647_filesep = File.separator;
     for (i=0; i<ROIarr.length-1; i++) 
         str = str + ROIarr[i] + filesep; //literally no idea what's going on from 14-17
     str = str + ROIarr[ROIarr.length-1] + filesep; 
     
print("af647 Image path= " + AF647_imagepath)//writes the raw image path and then the path that it follows ((3)
print("ROI path= " + ROI); //see above, but for ROI path (18)
output_647 = dir + filesep 
print("Output path= " + output_647)
transformationfile = filesep + str + "TransformationMatrices.txt";

AF647_results = output_647 + "AF647_results.csv";
fibril_results = fib_488_output + "Fibril_results.csv";

AF647_image = output_647 + "AF647.tif";

run("Set Measurements...", "centroid stack redirect=None decimal=3");

//define background images for each channel (made and saved before)
Hsp_background_image = filesep + str + "Hsp_background.tif"
fib_background_image= filesep + str + "647_background.tif"
AF647_background_image= AF647_filesep + str + "AF647_background.tif"


Hsp = fib_488_output + "MAX_Hsp.tif";
fibril = fib_488_output + "fibril_MAX.tif";
AF647= output_647+ "AF647_MAX.tif"



//define location for fibril colocalisation files to go to. first this folder is colocalisation analysis, then the 647 colocalisation with fibrils, and then the AF488 and AF647 colocalisation
Fibcolocalpath = filesep + str + "/Colocalisation_analysis/";
print(Fibcolocalpath);
File.makeDirectory(Fibcolocalpath);

coloc_fibril_488_path= Fibcolocalpath + "/fibrils_AF488_coloc/"
print(coloc_fibril_488_path);
File.makeDirectory(coloc_fibril_488_path);

coloc_fibril_647_path= Fibcolocalpath + "/fibrils_AF647_coloc/"
print(coloc_fibril_647_path);
File.makeDirectory(coloc_fibril_647_path);

coloc_488_647_path= Fibcolocalpath + "/AF488_AF647_coloc/"
print(coloc_488_647_path);
File.makeDirectory(coloc_488_647_path);


// Fibril channel
open(path) //select raw image
Imagetitle =getTitle()
rename("Raw"); 

roiManager("Open", ROI);
roiManager("Select", 0);	//this will be the fibril channel
run("Duplicate...", "title=fibril duplicate");

//Background correction 
open(fib_background_image);
rename("Background_fib");
selectWindow("fibril");
run("32-bit");
run("Beam Profile Correction", "electronic_offset=387 background_image=Background_fib");
setMinAndMax(0, 65536);
run("16-bit");
roiManager("Select", 2);
run("Duplicate...", "title=fibril_background_corrected duplicate");
saveAs("tif", fib_488_output + "Background_corrected_fib.tif");
selectWindow("Background_fib")
close();

//this will save the movie that can be used to get the intensity trajectories

selectWindow("Background_corrected_fib.tif")
run("Z Project...", "stop=20 projection=[Max Intensity]");
run("Properties...", "channels=1 slices=1 frames=1 unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000");
saveAs("tif", fib_488_output + "fibril_MAX.tif");
close()
roiManager("reset");
selectWindow("Background_corrected_fib.tif");
close()

//the AF647-labelled protein background correction 
open(AF647_imagepath) //select raw image
Imagetitle =getTitle()
rename("AF647_image"); 

roiManager("Open", ROI);
roiManager("Select", 0);	//this will be the AF647 channel
run("Duplicate...", "title=AF647 duplicate");

//Background correction 
open(AF647_background_image);
rename("Background_AF647");
selectWindow("AF647");
run("32-bit");
run("Beam Profile Correction", "electronic_offset=500 background_image=Background_AF647");
setMinAndMax(0, 65536);
run("16-bit");
roiManager("Select", 2);
run("Duplicate...", "title=AF647_background_corrected duplicate");
saveAs("tif", output_647 + "Background_corrected_AF647.tif");
selectWindow("Background_AF647")
close();

//this will save the movie that can be used to get the intensity trajectories

selectWindow("Background_corrected_AF647.tif")
run("Z Project...", "stop=20 projection=[Max Intensity]");
run("Properties...", "channels=1 slices=1 frames=1 unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000");
saveAs("tif", output_647 + "AF647_MAX.tif");
close()
roiManager("reset");
selectWindow("Background_corrected_AF647.tif");
close()


//Hsp channel- this is a little bit tricker as you need to align the channel using bead matrix

selectWindow("Raw")

roiManager("Open", ROI);
roiManager("Select", 1);	//this will be the Hsp channel
run("Duplicate...", "title=Hsp duplicate");

//Background correction 
open(Hsp_background_image);
rename("Background_Hsp");
selectWindow("Hsp");
run("32-bit");
run("Beam Profile Correction", "electronic_offset=500 background_image=Background_Hsp");
setMinAndMax(0, 65536);
run("16-bit");

selectWindow("Hsp");
saveAs("tif", fib_488_output + "Background_corrected_Hsp.tif");
selectWindow("Raw");
close();

//alignment of Hsp channel using bead image
splitDir= dir + "/Hsp/"
print(splitDir);
File.makeDirectory(splitDir);
list = getFileList(dir); 
selectWindow("Background_corrected_Hsp.tif");
run("Image Sequence... ", "format=TIFF use save=splitDir");
selectWindow("Background_corrected_Hsp.tif");
close();

outputFolder = splitDir

function action(outputFolder, outputFolder, filename) {
        open(outputFolder + filename);  //+ filesep
        run("MultiStackReg", "stack_1="+filename+" action_1=[Load Transformation File] file_1=["+ transformationfile +"] stack_2=None action_2=Ignore file_2=[] transformation=Affine");
        roiManager("Select", 2);
        run("Duplicate...", "duplicate");
        saveAs("tif", outputFolder + filesep + filename);
        close();
        close();
}


input = "outputFolder";
output = "alignedoutputFolder";

list = getFileList(outputFolder);
for (i = 0; i < list.length; i++)
        action(outputFolder, outputFolder, list[i]);
        output = dir + filesep
run("Image Sequence...", "open=[outputFolder] sort");
rename("Hsp");
run("Z Project...", "projection=[Average Intensity]");
run("Properties...", "channels=1 slices=1 frames=1 unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000");
selectWindow("AVG_Hsp");
saveAs("tiff", fib_488_output+getTitle()); 
close()
selectWindow("Hsp");
saveAs("tiff", fib_488_output+getTitle()); 
run("Z Project...", "stop=20 projection=[Max Intensity] ");
run("Properties...", "channels=1 slices=1 frames=1 unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000");
saveAs("tiff", fib_488_output+getTitle()); 
selectWindow("MAX_Hsp.tif");
run("Enhance Contrast", "saturated=0.35");
run("Properties...", "channels=1 slices=1 frames=1 unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000");
selectWindow("Hsp.tif")
close()
roiManager("reset");


//NOW look at the fibril image (saved earlier as fibril_MAX.tif) and do ridge detection to idnetify them and find their coords

open(fibril);
run("Enhance Contrast", "saturated=0.35");
run("8-bit");
run("Ridge Detection", "line_width=3.50 high_contrast=230 low_contrast=87 extend_line displayresults add_to_manager method_for_overlap_resolution=NONE sigma=1.2 lower_threshold=10 upper_threshold=11 minimum_line_length=4 maximum=0");
waitForUser("Delete the junctions from the ROI manager before pressing OK")
selectWindow("fibril_MAX.tif");
close()
selectWindow("Summary");
saveAs("Results", fib_488_output + "Length_data.csv");
run("Close");
selectWindow("Results");
saveAs("Results", fib_488_output + "Fibril_results.csv");
run("Close");
selectWindow("Junctions");
run("Close");
print("\\Clear")
roiManager("reset");


//open the Hsp image (which we saved earlier, called Max_Hsp.tif)
open(Hsp);
run("Enhance Contrast", "saturated=0.35");
run("Properties...", "channels=1 slices=1 frames=1 unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000");
run("Clear Results")

selectWindow("MAX_Hsp.tif");


//Find the peaks within the Hsp image
run("Peak Finder");// need to check the threshold here
waitForUser("Have you selected your peaks? \n \nPress OK to continue");
roiManager("Measure");	//this will allow for a table to be produced
saveAs("Results", fib_488_output + "Hsp_results.csv"); 
close();
selectWindow("Results"); 
     run("Close"); 
roiManager("reset");


//open the AF647 image (saved earlier, which we called AF647_MAX.tif)
open(AF647);
run("Enhance Contrast", "saturated=0.35");
run("Properties...", "channels=1 slices=1 frames=1 unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000");
run("Clear Results")

selectWindow("AF647_MAX.tif");


//Find the peaks within the Hsp image
run("Peak Finder");// need to check the threshold here
waitForUser("Have you selected your peaks? \n \nPress OK to continue");
roiManager("Measure");	//this will allow for a table to be produced
saveAs("Results", output_647 + "AF647_results.csv"); 
close();
selectWindow("Results"); 
     run("Close"); 
roiManager("reset");

//colocalization analysis FOR FIBRILS and AF488 hsp's
open(fibril_results)
IJ.renameResults("Fibril");
open(Hsp_results)
IJ.renameResults("Hsp");
run("count colocalized peaks Fibril", "table_1=Fibril table_2=Hsp maximum_distance=3");
selectWindow("Colocalization");
saveAs("tiff", fib_488_output + "Colocalization.tiff");
close()
selectWindow("Fibril");
saveAs("Results", fib_488_output + "Fibril_colocalisation.csv"); 
saveAs("Results", coloc_fibril_488_path + Imagetitle + ".csv");
IJ.renameResults("Results");
selectWindow("Hsp");
		run("Close")
runMacro("C:/Users/ljr979/Desktop/Fiji.app/macros/Filter_for_colocal.ijm")
IJ.renameResults("Fibril_colocalisation_only")
selectWindow("Fibril_colocalisation_only");
saveAs("Results", fib_488_output + "Fibril_colocalisation_only.csv"); 
saveAs("Results", coloc_fibril_488_path + Imagetitle + "_colocalized.csv");
run("Close");
roiManager("reset");


//now do colocalization analysis for fibrils and AF647 proteins
open(fibril_results)
IJ.renameResults("Fibril");
open(AF647_results)
IJ.renameResults("AF647");
run("count colocalized peaks Fibril", "table_1=Fibril table_2=AF647 maximum_distance=3");
selectWindow("Colocalization");
saveAs("tiff", output_647 + "Colocalization.tiff");
close()
selectWindow("Fibril");
saveAs("Results", output_647 + "Fibril_colocalisation.csv"); 
saveAs("Results", coloc_fibril_647_path + Imagetitle + ".csv");
IJ.renameResults("Results");
selectWindow("AF647");
		run("Close")
runMacro("C:/Users/ljr979/Desktop/Fiji.app/macros/Filter_for_colocal.ijm")
IJ.renameResults("Fibril_colocalisation_only")
selectWindow("Fibril_colocalisation_only");
saveAs("Results", output_647 + "Fibril_colocalisation_only.csv"); 
saveAs("Results", coloc_fibril_647_path + Imagetitle + "_colocalized.csv");
run("Close");
roiManager("reset");




//now do colocalization analysis for AF647 and AF488 proteins

open(AF647_results);	//this needs to be changed, need to put the path into the begginning of the macro
//selectWindow("client_results.csv");
IJ.renameResults("AF647");
open(Hsp_results);	//this needs to be changed, need to put the path into the begginning of the macro
//selectWindow("Hsp_results.csv");
IJ.renameResults("Hsp");
run("count colocalized peaks Andrew", "table_1=AF647 table_2=Hsp maximum_distance=3");

selectWindow("AF647");
saveAs("Results", output_647 + "AF647_colocalisation.csv"); 
saveAs("Results", coloc_488_647_path + Imagetitle + ".csv");
IJ.renameResults("Results");
selectWindow("AF647");
		run("Close")
runMacro("C:/Users/ljr979/Desktop/Luke/fiji-win32/Fiji.app/macros/Filter_for_colocal.ijm")
selectWindow("Colocalization");
saveAs("tiff", coloc_488_647_path + "Colocalization.tiff");
close()

IJ.renameResults("AF647_colocalisation")
selectWindow("AF647_colocalisation")
saveAs("Results", coloc_488_647_path + Imagetitle + "_AF647_colocalisation.csv"); 
selectWindow("Hsp");
     run("Close"); 
selectWindow("fibril");
     run("Close"); 

}


