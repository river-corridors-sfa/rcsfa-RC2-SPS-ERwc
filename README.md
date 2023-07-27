# YRB_Water_Column_Respiration
(Incomplete DRAFT)

A public GitHub repository containing data and code (R scripts) to generate the statistics and plots (i.e., figures) to be published in the manuscript "Yakima River Basin Water Column Respiration is a Minor Component of River Ecosystem Respiration" for Fulton et al. (20XX). See Fulton et al. (20XX) for references cited below. The four R scripts and input files (i.e., data) listed below will run all the statistical analyses and plot the figures described in the following Methods sections of the manuscript:

Section 2.2 Watershed characterization and site selection

Script: 1_DataPreparation.R: First step in a 2-step to prep data and run cluster analysis grouping NHDPlusV2 catchments into classes that share similar biophysical characteristics. This first R script processes geospatial data layers for input into R script "2_ClusterAnalysis.R" (see below).
Input files and data sources: See Methods Section 3.2 for publicly available data sources and web links.

Script: 2_ClusterAnalysis.R: Second step in a 2-step process to prep data and run cluster analysis. This R script runs the cluster analysis and XXXXXX (Steph).
Input files: Describe the dataset (Steph). Source: Where can this data be found? (Steph)

Section 2.5 OM chemistry via ultra-high resolution mass spectrometry and biochemical transformations

Script: Transformation_Analysis_by_Sample_sapply.R: Describe what the script does.
Input files: 
-Processed_RC2_SPS_11-2021_Data_Clean_VGC.csv: Describe the dataset (VGC). Source: Where can this data be found? (VGC)
-Transformation_Database_07-2020.csv: Describe dataset (VGC). Source: Where can this data be found? (VGC)

Section 2.6 Relationship of water column respiration rates to watershed characteristics and surface water chemistry

Script: RC2_spatial_study_MLR_v3.R, which runs the multiple linear regression models and plots Figures 3, 4, and 5 (Sect. 3.1); and generates kernel density statistics and plots Figure 6 (Sect. 3.2) in the Results section of the manuscript.
Input files: 
-spatial_data.csv: Describe the dataset (Xinming). Source: Where can this data be found? (Xinming)
-SPS_Total_and_Normalized_Transformations_01-03-23.csv: Output file (VGC???) from R script "Transformation_Analysis_by_Sample_sapply.R" (see #2, above). Source: Where can this data be found? (VGC). 
-v2_SPS_NPOC_TN_DIC_TSS_Ions_Summary.csv: Water chemistry data from 2021 spatial study (NOTE: include exact variables names from data package in the R script). Source: Grieger, S., Barnes, M., Borton, M. A., Chen, X., Chu, R., Farris, Y., et al.: Spatial Study 2021: Sample-Based Surface Water Chemistry and Organic Matter Characterization across Watersheds in the Yakima River Basin, Washington, USA (v2) [dataset], 10.15485/1898914, 2022.
-mean_ERvolumetric_best_streamPULSEsites.csv: Output file of whole river ecosystem respiration (ERtot) from R script "Appling_ER_data.R" (see #4, below). These data are a subset of 208 of the 222 StreamPULSE sites from Bernhardt et al. (2022) that have been subject to rigorous quality control to cull sites potentially affected by process or observation error. StreamPULSE sites where 0 <= ERtot All ERtot values >= 0 mg O2 L^-1 d^-1 have been converted to zero (i.e., 0.00). See Methods Sect. 2.7 of the manuscript for a more detailed description of the derivation of the dataset. 
-ERwc_combined_lit_valuesV2.csv: Published values of water column respiration rates (ERwc). Source: Table 2, Methods Sect. 3.1 of the manuscript.


Section 2.7 Comparison to published water column respiration rates

Script: Appling_ER_data.R: Describe what the script does (Steph).
Input files: 
-daily_predictions.csv: Predicted daily values of GPP and ER for 356 sites across the continental United States (CONUS) from Appling et al. (2018b, c). Source: Appling AP, Read JS, Winslow LA, Arroita M, Bernhardt ES, Griffiths NA, et al. Metabolism estimates for 356 U.S. rivers (2007-2017): U.S. Geological Survey data release. 2018c. https://doi.org/10.5066/F70864KX
-streampulse_synthesis_statset.csv: Predicted values of GPP and ER for 222 sites from across the CONUS from Bernhardt et al. (2022) (i.e., "StreamPULSE" sites). StreamPULSE sites are a subset of sites from Appling et al. (2018b) that have been subject to rigorous quality control to cull sites potentially affected by process or observation error. Source: Dr. Robert Hall, Jr., Flathead Lake Biological Station, University of Montana, Polson, MT, USA.
