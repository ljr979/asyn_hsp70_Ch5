macro "AF488 channel + ASCP labelled fibrils" {

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
output = dir + filesep 
print("Output path= " + output)
transformationfile = filesep + str + "TransformationMatrices.txt";
Hsp_results = output + "Hsp_results.csv";
fibril_results = output + "Fibril_results.csv";
Hsp_image = output + "Hsp.tif";

run("Set Measurements...", "centroid stack redirect=None decimal=3");

//defining the images that are used at the background & image to find peaks for each channel
Hsp_background_image = filesep + str + "Hsp_background.tif"
fib_background_image= filesep + str + "647_background.tif"
Hsp = output + "MAX_Hsp.tif";
fibril = output + "fibril_MAX.tif";
//create folder to save colocalisation in
Fibcolocalpath = filesep + str + "/Colocalisation_analysis/";
print(Fibcolocalpath);
File.makeDirectory(Fibcolocalpath);


// Fibril channel
open(path) //select raw image with fibril in one channel and 488 in other
Imagetitle =getTitle()
rename("Raw"); 
//open ROI.zip file you defined earlier
roiManager("Open", ROI);
//crop to just fibril channel
roiManager("Select", 0);	//this will be the fibril channel
run("Duplicate...", "title=fibril duplicate");

//Background correction 
open(fib_background_image);
rename("Background_fib");
selectWindow("fibril");
run("32-bit");
//correct beam profile using 'background' image
run("Beam Profile Correction", "electronic_offset=750 background_image=Background_fib");
setMinAndMax(0, 65536);
run("16-bit");
//crop edges 
roiManager("Select", 2);
run("Duplicate...", "title=fibril_background_corrected duplicate");
saveAs("tif", output + "Background_corrected_fib.tif");
selectWindow("Background_fib")
close();
//re-open this corrected image
selectWindow("Background_corrected_fib.tif")
//create z project / maximum so that fibrils are bright and easily defined for ridge detection
run("Z Project...", "stop=20 projection=[Max Intensity]");
run("Properties...", "channels=1 slices=1 frames=1 unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000");
saveAs("tif", output + "fibril_MAX.tif");
close()
roiManager("reset");
selectWindow("Background_corrected_fib.tif");
close()


//Hsp channel- this is a little bit tricker as you need to align the channel using bead matrix
selectWindow("Raw")
roiManager("Open", ROI);
//now do the same with other ROI which is AF488 channel 
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
saveAs("tif", output + "Background_corrected_Hsp.tif");
selectWindow("Raw");
close();

//alignment of Hsp channel using bead image
//first split the STACK of AF488 channel into individual images, and save separately in 'Hsp' folder
splitDir= dir + "/Hsp/"
print(splitDir);
File.makeDirectory(splitDir);
list = getFileList(dir); 
selectWindow("Background_corrected_Hsp.tif");
run("Image Sequence... ", "format=TIFF use save=splitDir");
selectWindow("Background_corrected_Hsp.tif");
close();

outputFolder = splitDir
//define function which aligns each individual frame from the Hsp stack, to the 647 channel.
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
//runs the function on the stack
list = getFileList(outputFolder);
for (i = 0; i < list.length; i++)
        action(outputFolder, outputFolder, list[i]);
        output = dir + filesep
run("Image Sequence...", "open=[outputFolder] sort");
rename("Hsp");
//make a z stack of the ALIGNED frames
run("Z Project...", "projection=[Average Intensity]");
run("Properties...", "channels=1 slices=1 frames=1 unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000");
selectWindow("AVG_Hsp");
saveAs("tiff", output+getTitle()); 
close()

selectWindow("Hsp");
saveAs("tiff", output+getTitle()); 
//maximum intensity projection to identify peaks
run("Z Project...", "stop=20 projection=[Max Intensity] ");
run("Properties...", "channels=1 slices=1 frames=1 unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000");
saveAs("tiff", output+getTitle()); 
selectWindow("MAX_Hsp.tif");


run("Enhance Contrast", "saturated=0.35");
run("Properties...", "channels=1 slices=1 frames=1 unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000");
close()
selectWindow("Hsp.tif")
close()
roiManager("reset");

//now find the fibril ridges
open(fibril);
run("Enhance Contrast", "saturated=0.35");
run("8-bit");
//run ridge detection (check this manually first for appropriate parameters)
run("Ridge Detection", "line_width=3.50 high_contrast=230 low_contrast=87 extend_line displayresults add_to_manager method_for_overlap_resolution=NONE sigma=1.2 lower_threshold=10 upper_threshold=11 minimum_line_length=4 maximum=0");
waitForUser("Delete the junctions from the ROI manager before pressing OK")
selectWindow("fibril_MAX.tif");
close()
selectWindow("Summary");
//saves lengths of fibrils in pixels
saveAs("Results", output + "Length_data.csv");
run("Close");
selectWindow("Results");
saveAs("Results", output + "Fibril_results.csv");
run("Close");
selectWindow("Junctions");
run("Close");
print("\\Clear")
roiManager("reset");

//now find chaperones in 488 channel
open(Hsp);
rename("Hsp");
run("Enhance Contrast", "saturated=0.35");
run("Properties...", "channels=1 slices=1 frames=1 unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000");
run("Clear Results")

//Find the peaks within the Hsp image
run("Peak Finder");// need to check the threshold here
waitForUser("Have you selected your peaks? \n \nPress OK to continue");
roiManager("Measure");	//this will allow for a table to be produced
saveAs("Results", output + "Hsp_results.csv"); 
close();
selectWindow("Results"); 
     run("Close"); 
roiManager("reset");


//colocalization analysis
open(fibril_results)
IJ.renameResults("Fibril");
open("Hsp_results.csv")
IJ.renameResults("Hsp");
//this next line runs the plugin that finds 488-channel peaks that are colocalised with the fibril pixels
run("count colocalized peaks Fibril", "table_1=Fibril table_2=Hsp maximum_distance=3");
selectWindow("Colocalization");
saveAs("tiff", output + "Colocalization.tiff");
close()
//save all your results
selectWindow("Fibril");
saveAs("Results", output + "Fibril_colocalisation.csv"); 
saveAs("Results", Fibcolocalpath + Imagetitle + ".csv");
IJ.renameResults("Results");
selectWindow("Hsp");
		run("Close")

		//filter this table for only those colocalised w fibrils
runMacro("C:/Users/ljr979/Desktop/Fiji.app/macros/Filter_for_colocal.ijm")
IJ.renameResults("Fibril_colocalisation_only")
selectWindow("Fibril_colocalisation_only");
saveAs("Results", output + "Fibril_colocalisation_only.csv"); 
saveAs("Results", Fibcolocalpath + Imagetitle + "_colocalized.csv");
run("Close");
roiManager("reset");

//now find trajectories of all of the af488 molecules and save them to this folder
coloc_trajectory_output = filesep + str + "/coloc_trajectories/"
print(coloc_trajectory_output);
File.makeDirectory(coloc_trajectory_output);

//now find the coordinates of the hsps colocalised with fibrils, and get the trajectories at these coordinates
open(Hsp_image)
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
//
selectWindow("Hsp.tif");
		run("Close");

}


