
Readme file for the CRB Catchment Cluster Analysis in R.
R version 3.5.3 (2019-03-11) -- "Great Truth"


Brief Description:
There are two R Scripts that will allow the user to run the cluster analysis:
	>> The "1_DataPreparation.R" script will read the data into R and the user can choose which subset to use for the analysis
	>> The "2_ClusterAnalysis.R" script will run the kmeans cluster analysis


Required Packages:
tidyverse
sf


Required Data Files:
SEE "R_Inputs" FOLDER


User Instructions:
1) Open the "1_DataPreparation.R" script and update the input_path (Line 8) to correspond with the directory where all "R_Inputs" files are located on your machine

2) Run Lines 3 through 14 - This will read data into R and load the necessary packages
	>> NOTE: The 'centers' data that is read in correspond to the K-means centers that were used for the previous runs of this analysis. Therefore, if you would like to reproduce these results, you can use the data from 'centers'. If you would like R to re-generate new centers, or run the analysis with a different number of clusters, then the 'centers' data will not be necessary. This will be explained further in Step 7. 

3) Determine which subset you want to run the analysis on
	Set 0 = All statistical moments for all variables
	Set 1 = All variables, MEAN statistical moments only
	Set 2 = Set 1, but only months 3, 7, and 11 for FPAR, LAI, ET, and PPT
	Set 3 = Set 1, but only annual median for FPAR, LAI, ET, and PPT
	Set 4 = Set 3 with only FPAR and PPT
	Set 5 = Set 3 with only LAI and PPT
	Set 6 = Set 3 with only ET and PPT
	Set 7 = All variables, STD statistical moments only
	Set 8 = All variables, MEAN and STD statistical moments only
	Set 9 = Set 2, but only FPAR and PPT
	Set 10 = Set 2, but only LAI and PPT
	Set 11 = Set 2, but only ET and PPT

4) Run the script for the chosen subset

5) Open the "2_ClusterAnalysis.R" script and update the input_path (Line 8) to correspond with the directory where all "R_Inputs" files are located on your machine. Update the output_path (Line 9) to correspond with the directory where analysis results will be stored on your machine. 

6) Define the model parameters and data (Lines 13 through 16): 
	>> "n_clusters" is the number of clusters you would like to use (NOTE: If using the pre-defined centers file, this should be 6)
	>> "Data" should be the the data set that you chose (e.g., Data_norm_8)
	>> "Centers" can be one of two things. If you want to use the pre-defined centers file, then this should be set to that file (e.g., centers_8). However, if you want R to regenerate centers, then you simply need to set this as whatever number of clusters you want to run (e.g., 6)

7) Run the cluster analysis (Lines 20 through 23)

8) Extract the results (Lines 27 through 45)

9) Join cluster data to the shapefile (Lines 48 through 51)

10) Save results (Lines 57 through 69)
	>> NOTE: Running these lines will create a new folder in your "R_Outputs" directory and all results files will be written in that new folder
	>> NOTE: If you save multiple results from multiple runs in a single day, you could risk overwriting the previous results. Therefore, once results are written, be sure to save a backup copy elsewhere. If you run analyses on different days, then a new folder with that date will be generated, and those new results will go in that new folder. 








