# YRB_Water_Column_Respiration
(DRAFT)

Summary of files in this repository

This repository contains code (R scripts) and data associated with the manuscript “Yakima River Basin Water Column Respiration is a Minor Component of River Ecosystem Respiration” submitted to EGU Biogeochemistry (Fulton et al. 202X; submitted to Geochemistry). 

The repository is comprised of three main folders, one containing five analysis-specific R script subfolders, one containing six data-specific subfolders, and one containing map figures for the manuscript. Related code and data subfolders use the same naming convention to easily identify the data input files associated with each R script.

Code subfolders listed below contain the following R scripts:

•	Appling_ERtot_analysis
1)	Appling_ERtot_analysis.R: R script to process and analyze daily reach-averaged ecosystem respiration (ERtot) data for 356 rivers and streams as described in Appling et al. (2018b).  This script creates a subset of ERtot data from 208 StreamPULSE sites and converts these data from the original areal units (g O2 m^-2 d^-1) to volumetric units (mg O2 L^-1 d^-1) for comparison with water column respiration (ERwc; mg O2 L^-1 d^-1) data observed in the Yakima River basin (this study).
2)	Creates an output data file of ERtot values used as an R input data file to generate descriptive statistics and the kernel density plot in Figure 5b. in the as discussed in the manuscript and use in plotting Figure 6.

•	Map_Layers
1)	ERwc_Manuscript_Map.R: R script to create map figures Figure 1 and Figure S1 using site location and water column respiration (ERwc) data from this study and publicly available geospatial data.

•	Multiple_linear_regression
1)	RC2_spatial_study_MLR.R: R script to process and analyze ERwc data observed in the Yakima River basin (this study), ERwc values from the published literature, and reach-averaged ecosystem respiration data from Appling et al. (2018b,c). The script performs 1) a multiple linear regression using step-wise forward selection to evaluate the relationship between ERwc and physical and chemical surface water parameters and watershed characteristics; 2) generates kernel density plots and/or descriptive statistics and for ERwc rates observed in the Yakima River basin (this study), ERwc rates from the published literature, and ERtot rates from 208 rivers and streams across the CONUS from Appling et al. (2018b,c); and 3) plots and creates Figures 3-6, and Figure S1 in the manuscript and Supplementary Information.

•	OM_transformation_analysis
1)	Transformation_Analysis_by_Sample_sapply.R: R code to generate organic matter (OM) biochemical transformation data using Fourier-transform ion cyclotron resonance mass spectrometry (FTICR-MS) data as described in the manuscript.
2)	Creates an output data file of OM transformation data used as an R script input data file in the multiple regression analysis (i.e., RC2_spatial_study_MLR.R).

•	Site_selection_cluster_analysis
1)	1_DataPreparation.R: The first of two R codes to run the Columbia River basin catchment cluster analysis as described in the manuscript. This R code will read the data into R and the user can choose which variable set to use for the cluster analysis performed in the R script 2_ClusterAnalysis.R.
2)	2_ClusterAnalysis.R: The second of two R codes to run the Columbia River basin catchment cluster analysis to group Columbia River NHDPlusV2.1 catchments into six classes sharing similar characteristics as described in the manuscript. This R code runs the k-means cluster analysis using the normalized data for the variable-statistical moment combinations of 16 key biophysical and hydrological attributes.

Data subfolders listed below contain the following datasets required by each R script:

•	Appling_ERtot_analysis
1)	daily_predictions_ER_depth.csv: Stream depth and predicted daily reach-scale river ecosystem respiration data (ERtot) for 490,907 observations on 356 rivers and streams as described in Appling et al. (2018b,c). 
2)	StreamPULSE_bestSiteIDs.csv: Site information data for 222 StreamPULSE sites as described in Bernhardt et al. (2022). StreamPULSE sites are a subset of sites from the Appling et al. (2018b,c) dataset for 222 rivers across the CONUS created by Bernhardt et al. (2022).
3)	mean_ERtot_bestSiteIDs.csv: Data output file containing mean stream depth (m) and ERtot in areal (g O2 m^-2 d^-1) and volumetric (mg O2 L^-1 d^-1) units for 208 StreamPULSE sites on streams and rivers across the CONUS using data from Appling et al. (2018a,b) and Bernhardt et al. (2022) as described in the manuscript.

•	OM_transformation_analysis
1)	SPS_Total_and_Normalized_Transformations_01-03-23.csv: (VGC needs to describe)
2)	Processed_RC2_SPS_11-2021_Data_Clean_VGC.csv: (VGC needs to describe)
3)	Transformation_Database_07-2020.csv: Data file containing the name or molecular formula and mass for 1,255 commonly observed biochemical transformations.

•	Multiple_linear_regression
1)	Minidot_Summary_Statistics_v2.csv: (Xinming needs to move this file from 2021_spatial_study_data folder and edit code, then delete the 2021_spatial_study_data if no longer needed) 
2)	spatial_data.csv: Biological, physical, and geospatial data used to evaluate drivers of spatial variation in water column respiration at 47 sites throughout in the Yakima River basin, Washington, USA. Dataset includes mean values for water column respiration (ERwc); water temperature; upstream drainage area; and stream order for each site. ERwc and temperature data are from the full published sensor dataset (Fulton et al. 2022; https://data.ess-dive.lbl.gov/view/doi:10.15485/1892052).
3)	SPS_Total_and_Normalized_Transformations_01-03-23.csv: Output file from R script for the OM transformation analysis (see Transformation_Analysis_by_Sample_sapply.R). 
4)	v2_SPS_NPOC_TN_DIC_TSS_Ions_Summary.csv: Surface water chemistry data used to evaluate drivers of spatial variation in water column respiration at 47 sites throughout in the Yakima River basin, Washington, USA (this study). Dataset includes mean values for nitrate (NO3-N; mg L^-1), dissolved organic carbon (DOC; mg L^-1), dissolved inorganic carbon (DIC; mg L^-1), and total suspended solids (TSS; mg L^-1). Mean data values were calculated from the full published surface water chemistry dataset (Grieger et al. 2022; https://data.ess-dive.lbl.gov/view/doi:10.15485/1898914).

•	Water_column_respiration_published
1)	Water_column_respiration_published_values.csv: Mean water column respiration rate (ERwc) data from the published literature for 25 sites along rivers across the conterminous United States and the Amazon River Basin (see Table 2 in Fulton et al. 202X).

•	Site_selection_cluster_analysis
1)	AllDataNormalized_zonal.csv: R code data input file containing the normalized values for the variable-statistical moment combinations used to group NHDPlusV2.1 catchments in the Columbia River basin into six classes sharing similar characteristics using cluster analysis as described in the manuscript. The variable-statistical moment combinations consist of four statistical moments (mean; minimum; maximum; and standard deviation) calculated for each of 16 key biophysical and hydrological attributes. The 16 key biophysical and hydrological attributes included climate, vegetation structure and function, topography, and wildfire potential that were generated using existing, readily available geospatial data. NOTE: The file is too large to store on GitHub. Download the file directly from ESS-DIVE at (URL or DOI or cite data package directly: Fulton et al. 202Xx). 
2)	6clusters_centers.csv: R code data input file containing the six cluster centroid values for each variable-statistical moment combination. The k-means clustering algorithm uses these six cluster centroids to group NHDPlusV2.1 catchments in the Columbia River basin into six classes sharing similar characteristics as described in the manuscript.
3)	CB_catchments_clustered: Subfolder containing CB_catchments_clustered.shp and related GIS files (*.shx, *.dbf, *.prj, *.sbn, *.sbx, and *.xml) for the GIS shapefile of all the NHDPlusV2.1 catchments for the Columbia River basin (Need Jerry to clarify if/how it’s used in cluster analysis).

•	Map_Layers: Subfolder containing the QGIS project and subfolders for the individual GIS data layers used to create Figure 1 and Figure S1 and archived intermediate figure files.
1)	Spatial21_relief.qgis: QGIS file to….(Brie needs to define)
2)	dem30m_yrb: (Brie needs from Kyongho)
3)	landcover16_yk: https://www.usgs.gov/centers/eros/science/national-land-cover-database
4)	SPS_Sampling_Sites: GIS shapefile of the 47 sampling locations for this study created using the published metadata (cite data package).
5)	tl_2021_us_state: https://www2.census.gov/geo/tiger/TIGER2021/STATE/
6)	YakimaRiverBasin_Boundary: https://ecology.wa.gov/research-data/data-resources/geographic-information-systems-gis/data
7)	YRB_Cluster: Catchment cluster map for the Yakima River basin clipped from the Columbia River basin catchment cluster map parent file. 
8)	YRB_River_4326: (Brie needs from Sophia)
9)	Archive_Intermediate_Files: Subfolder containing the following files:
    (a) SPS_ER_Water_Column_Map.pdf
  	(b) SPS_ER_Water_Column_Map_Cluster.pdf
  	(c) SPS_ER_Water_Column_Map_Formatted.pdf

The main Figure folder listed below comprises the following three map figures from the manuscript:
1)	Figure1_LandUse_LandCover_Relief_Maps.pdf: PDF file of Figure 1 in the manuscript.
2)	Figure1_LandUse_LandCover_Relief_Maps.png: PNG file of Figure 1 in the manuscript.
3)	FigureS1_Cluster_Respiration_Maps.pdf: PDF file of Figure S1 from the Supplementary Information section of the manuscript.

Citations

Appling, A. P., Hall, R. O., Yackulic, C. B., and Arroita, M.: Overcoming Equifinality: Leveraging Long Time Series for Stream Metabolism Estimation, Journal of Geophysical Research: Biogeosciences, 123, 624-645, 10.1002/2017jg004140, 2018a.

Appling, A. P., Read, J. S., Winslow, L. A., Arroita, M., Bernhardt, E. S., Griffiths, N. A., Hall, R. O., Jr., Harvey, J. W., Heffernan, J. B., Stanley, E. H., Stets, E. G., and Yackulic, C. B.: The metabolic regimes of 356 rivers in the United States, Scientific Data, 5, 180292, 10.1038/sdata.2018.292, 2018b.

Appling, A. P., Read, J. S., Winslow, L. A., Arroita, M., Bernhardt, E. S., Griffiths, N. A., Hall, R. O., Jr., Harvey, J. W., Heffernan, J. B., Stanley, E. H., Stets, E. G., and Yackulic, C. B.: Metabolism estimates for 356 U.S. rivers (2007-2017): U.S. Geological Survey data release [dataset], https://doi.org/10.5066/F70864KX, 2018c.

Bernhardt, E. S., Savoy, P., Vlah, M. J., Appling, A. P., Koenig, L. E., Hall, R. O., Arroita, M., Blaszczak, J. R., Carter, A. M., Cohen, M., Harvey, J. W., Heffernan, J. B., Helton, A. M., Hosen, J. D., Kirk, L., McDowell, W. H., Stanley, E. H., Yackulic, C. B., and Grimm, N. B.: Light and flow regimes regulate the metabolism of rivers, Proceedings of the National Academy of Sciences, 119, e2121976119, 10.1073/pnas.2121976119, 2022.

Fulton, S. G., Barnes, M., Borton, M. A., Chen, X., Farris, Y., Forbes, B., Garayburu-Caruso, V. A., Goldman, A. E., Grieger, S., Hall, R. O., Jr., Kaufman, M. H., Lin, X., McCann, E., McKever, S. A., Myers-Pigg, A., Otenburg, O., Pelly, A. C., Ren, H., Renteria, L., Scheibe, T. D., Son, K., Tagestad, J., Torgeson, J. M., and Stegen, J. C.: Yakima River Basin Water Column Respiration is a Minor Component of River Ecosystem Respiration, Submitted to Biogeochemistry, XXXXX.

Fulton, S. G., Barnes, M., Borton, M. A., Chen, X., Farris, Y., Forbes, B., Garayburu-Caruso, V. A., Goldman, A. E., Grieger, S., Kaufman, M. H., Lin, X., McKever, S. A., Myers-Pigg, A., Otenburg, O., Pelly, A., Ren, H., Renteria, L., Scheibe, T. D., Son, K., Torgeson, J. M., and Stegen, J. C.: Spatial Study 2021: Sensor-Based Time Series of Surface Water Temperature, Specific Conductance, Total Dissolved Solids, Turbidity, pH, and Dissolved Oxygen from across Multiple Watersheds in the Yakima River Basin, Washington, USA (v2) [dataset], 10.15485/1892052, 2022.

Grieger, S., Barnes, M., Borton, M. A., Chen, X., Chu, R., Farris, Y., Forbes, B., Fulton, S. G., Garayburu-Caruso, V. A., Goldman, A. E., Gonzalez, B. I., Kaufman, M. H., McKever, S. A., Myers-Pigg, A., Otenburg, O., Pelly, A., Renteria, L., Scheibe, T. D., Son, K., Torgeson, J. M., Toyoda, J. G., and Stegen, J. C.: Spatial Study 2021: Sample-Based Surface Water Chemistry and Organic Matter Characterization across Watersheds in the Yakima River Basin, Washington, USA (v2) [dataset], 10.15485/1898914, 2022.
