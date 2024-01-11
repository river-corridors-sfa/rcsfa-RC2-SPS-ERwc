# ------------------------------------------------------------------
#
# Script name: 1_DataPreparation.R
#
# Author: Erin L McCann
#
# Contact: erin.mccann@pnnl.gov
#
# Date Created: 2020-02-04
#
# Description: Read in AllDataNormalized_zonal.csv and cluster_centers.csv files 
#              Create subsets to be used for cluster analysis.
#              Set descriptions are commented below and can be found in the README
#
# -------------------------------------------------------------------------


# Load Packages:
library(tidyverse)

# Set working directory
current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd("../..")
getwd()

# Define Input Data Path:
input_path = "./Data/Site_selection_cluster_analysis/" # Directory where input files are located


# Read in Normalized Data and Centers:
Data_norm = read.csv(paste0(input_path, "AllDataNormalized_zonal.csv")) # due to the file size, this file is not included in the GitHub repository. it can be located in the ESS-DIVE Data Package
centers = read.csv(paste0(input_path, "6clusters_centers.csv"))



#---------------------------------------------------------------------
#----------------------- Create Data Subset: -------------------------
#---------------------------------------------------------------------



# SET 0 --------------------------------------------------------------

## Set 0 Description: All statistical moments for all variables
Data_norm_0 = Data_norm %>% 
  select(-AreaSqKM)
centers_0 = centers



# SET 1 --------------------------------------------------------------

## Set 1 Description: Using only MEAN moments for all variables
Data_norm_1 = Data_norm %>% 
  select(catchment_ID, ends_with("MEAN"))
centers_1 = centers %>% 
  select(ends_with("MEAN"))



# SET 2 --------------------------------------------------------------

## Set 2 Description: Using only MEAN moments for all variables, 
# but only months 3, 7, and 11 for FPAR, LAI, ET, and PPT
Data_norm_2 = Data_norm %>% 
  select(catchment_ID, ends_with("MEAN"), -c(contains("_1_MEAN"), contains("_01_MEAN"), contains("2_MEAN"), 
                                             contains("4_MEAN"), contains("5_MEAN"), contains("6_MEAN"), 
                                             contains("8_MEAN"), contains("9_MEAN"), contains("10_MEAN"), 
                                             contains("12_MEAN"), contains("norm_annual_MEAN"),
                                             contains("2009_2019_median_MEAN")))
centers_2 = centers %>% 
  select(ends_with("MEAN"), -c(contains("_1_MEAN"), contains("_01_MEAN"), contains("2_MEAN"), 
                               contains("4_MEAN"), contains("5_MEAN"), contains("6_MEAN"), 
                               contains("8_MEAN"), contains("9_MEAN"), contains("10_MEAN"), 
                               contains("12_MEAN"), contains("norm_annual_MEAN"),
                               contains("2009_2019_median_MEAN")))



# SET 3 --------------------------------------------------------------

## Set 3 Description: Using only MEAN moments for all variables,
# but only the annual median for FPAR, LAI, ET, and PPT
Data_norm_3 = Data_norm %>% 
  select(catchment_ID, ends_with("MEAN"), -c(contains("1_MEAN"), contains("2_MEAN"), contains("3_MEAN"), 
                                             contains("4_MEAN"), contains("5_MEAN"), contains("6_MEAN"),
                                             contains("_7_MEAN"), contains("07_MEAN"), contains("8_MEAN"), 
                                             contains("9_MEAN"), contains("10_MEAN"), contains("12_MEAN")))
centers_3 = centers %>% 
  select(ends_with("MEAN"), -c(contains("1_MEAN"), contains("2_MEAN"), contains("3_MEAN"), 
                               contains("4_MEAN"), contains("5_MEAN"), contains("6_MEAN"),
                               contains("_7_MEAN"), contains("07_MEAN"), contains("8_MEAN"), 
                               contains("9_MEAN"), contains("10_MEAN"), contains("12_MEAN")))



# SET 4 --------------------------------------------------------------

## Set 4 Description: Using only MEAN moments for all variables,
# but only the annual median for FPAR and PPT
Data_norm_4 = Data_norm %>% 
  select(catchment_ID, ends_with("MEAN"), -c(contains("1_MEAN"), contains("2_MEAN"), contains("3_MEAN"), 
                                             contains("4_MEAN"), contains("5_MEAN"), contains("6_MEAN"),
                                             contains("_7_MEAN"), contains("07_MEAN"), contains("8_MEAN"), 
                                             contains("9_MEAN"), contains("10_MEAN"), contains("12_MEAN"),
                                             contains("MCD15A3H_LAI_2009_2019_median_MEAN"),
                                             contains("MOD16A2_ET_2009_2019_median_MEAN")))
centers_4 = centers %>% 
  select(ends_with("MEAN"), -c(contains("1_MEAN"), contains("2_MEAN"), contains("3_MEAN"), 
                               contains("4_MEAN"), contains("5_MEAN"), contains("6_MEAN"),
                               contains("_7_MEAN"), contains("07_MEAN"), contains("8_MEAN"), 
                               contains("9_MEAN"), contains("10_MEAN"), contains("12_MEAN"),
                               contains("MCD15A3H_LAI_2009_2019_median_MEAN"),
                               contains("MOD16A2_ET_2009_2019_median_MEAN")))



# SET 5 --------------------------------------------------------------

## Set 5 Description: Using only MEAN moments for all variables,
# but only the annual median for LAI and PPT
Data_norm_5 = Data_norm %>% 
  select(catchment_ID, ends_with("MEAN"), -c(contains("1_MEAN"), contains("2_MEAN"), contains("3_MEAN"), 
                                             contains("4_MEAN"), contains("5_MEAN"), contains("6_MEAN"),
                                             contains("_7_MEAN"), contains("07_MEAN"), contains("8_MEAN"), 
                                             contains("9_MEAN"), contains("10_MEAN"), contains("12_MEAN"),
                                             contains("MCD15A3H_FPAR_2009_2019_median_MEAN"),
                                             contains("MOD16A2_ET_2009_2019_median_MEAN")))
centers_5 = centers %>% 
  select(ends_with("MEAN"), -c(contains("1_MEAN"), contains("2_MEAN"), contains("3_MEAN"), 
                               contains("4_MEAN"), contains("5_MEAN"), contains("6_MEAN"),
                               contains("_7_MEAN"), contains("07_MEAN"), contains("8_MEAN"), 
                               contains("9_MEAN"), contains("10_MEAN"), contains("12_MEAN"),
                               contains("MCD15A3H_FPAR_2009_2019_median_MEAN"),
                               contains("MOD16A2_ET_2009_2019_median_MEAN")))



# SET 6 --------------------------------------------------------------

## Set 6 Description: Using only MEAN moments for all variables,
# but only the annual median for ET and PPT
Data_norm_6 = Data_norm %>% 
  select(catchment_ID, ends_with("MEAN"), -c(contains("1_MEAN"), contains("2_MEAN"), contains("3_MEAN"), 
                                             contains("4_MEAN"), contains("5_MEAN"), contains("6_MEAN"),
                                             contains("_7_MEAN"), contains("07_MEAN"), contains("8_MEAN"), 
                                             contains("9_MEAN"), contains("10_MEAN"), contains("12_MEAN"),
                                             contains("MCD15A3H_FPAR_2009_2019_median_MEAN"),
                                             contains("MCD15A3H_LAI_2009_2019_median_MEAN")))
centers_6 = centers %>% 
  select(ends_with("MEAN"), -c(contains("1_MEAN"), contains("2_MEAN"), contains("3_MEAN"), 
                               contains("4_MEAN"), contains("5_MEAN"), contains("6_MEAN"),
                               contains("_7_MEAN"), contains("07_MEAN"), contains("8_MEAN"), 
                               contains("9_MEAN"), contains("10_MEAN"), contains("12_MEAN"),
                               contains("MCD15A3H_FPAR_2009_2019_median_MEAN"),
                               contains("MCD15A3H_LAI_2009_2019_median_MEAN")))



# SET 7 --------------------------------------------------------------

## Set 7 Description: Using only STD moments for all variables
Data_norm_7 = Data_norm %>% 
  select(catchment_ID, ends_with("STD"))
centers_7 = centers %>% 
  select(ends_with("STD"))



# SET 8 --------------------------------------------------------------

## Set 8 Description: Using only MEAN and STD moments for all variables
Data_norm_8 = Data_norm %>% 
  select(catchment_ID, ends_with("MEAN"), ends_with("STD"))
centers_8 = centers %>% 
  select(ends_with("MEAN"), ends_with("STD"))



# SET 9 --------------------------------------------------------------

## Set 9 Description: Using only MEAN moments for all variables, 
# but only months 3, 7, and 11 for FPAR and PPT
Data_norm_9 = Data_norm %>% 
  select(catchment_ID, ends_with("MEAN"), -c(contains("_1_MEAN"), contains("_01_MEAN"), contains("2_MEAN"),
                                             contains("ET_3"), contains("LAI_3"), contains("4_MEAN"), 
                                             contains("5_MEAN"), contains("6_MEAN"), contains("ET_7"),
                                             contains("LAI_7"), contains("8_MEAN"), contains("9_MEAN"), 
                                             contains("10_MEAN"), contains("ET_11"), contains("LAI_11"),
                                             contains("12_MEAN"), contains("norm_annual_MEAN"),
                                             contains("2009_2019_median_MEAN")))
centers_9 = centers %>% 
  select(ends_with("MEAN"), -c(contains("_1_MEAN"), contains("_01_MEAN"), contains("2_MEAN"),
                               contains("ET_3"), contains("LAI_3"), contains("4_MEAN"), 
                               contains("5_MEAN"), contains("6_MEAN"), contains("ET_7"),
                               contains("LAI_7"), contains("8_MEAN"), contains("9_MEAN"), 
                               contains("10_MEAN"), contains("ET_11"), contains("LAI_11"),
                               contains("12_MEAN"), contains("norm_annual_MEAN"),
                               contains("2009_2019_median_MEAN")))



# SET 10 --------------------------------------------------------------

## Set 10 Description: Using only MEAN moments for all variables, 
# but only months 3, 7, and 11 for LAI and PPT
Data_norm_10 = Data_norm %>% 
  select(catchment_ID, ends_with("MEAN"), -c(contains("_1_MEAN"), contains("_01_MEAN"), contains("2_MEAN"),
                                             contains("ET_3"), contains("FPAR_3"), contains("4_MEAN"), 
                                             contains("5_MEAN"), contains("6_MEAN"), contains("ET_7"),
                                             contains("FPAR_7"), contains("8_MEAN"), contains("9_MEAN"), 
                                             contains("10_MEAN"), contains("ET_11"), contains("FPAR_11"),
                                             contains("12_MEAN"), contains("norm_annual_MEAN"),
                                             contains("2009_2019_median_MEAN")))
centers_10 = centers %>% 
  select(ends_with("MEAN"), -c(contains("_1_MEAN"), contains("_01_MEAN"), contains("2_MEAN"),
                               contains("ET_3"), contains("FPAR_3"), contains("4_MEAN"), 
                               contains("5_MEAN"), contains("6_MEAN"), contains("ET_7"),
                               contains("FPAR_7"), contains("8_MEAN"), contains("9_MEAN"), 
                               contains("10_MEAN"), contains("ET_11"), contains("FPAR_11"),
                               contains("12_MEAN"), contains("norm_annual_MEAN"),
                               contains("2009_2019_median_MEAN")))



# SET 11 --------------------------------------------------------------

## Set 11 Description: Using only MEAN moments for all variables, 
# but only months 3, 7, and 11 for ET and PPT
Data_norm_11 = Data_norm %>% 
  select(catchment_ID, ends_with("MEAN"), -c(contains("_1_MEAN"), contains("_01_MEAN"), contains("2_MEAN"),
                                             contains("FPAR_3"), contains("LAI_3"), contains("4_MEAN"), 
                                             contains("5_MEAN"), contains("6_MEAN"), contains("FPAR_7"),
                                             contains("LAI_7"), contains("8_MEAN"), contains("9_MEAN"), 
                                             contains("10_MEAN"), contains("FPAR_11"), contains("LAI_11"),
                                             contains("12_MEAN"), contains("norm_annual_MEAN"),
                                             contains("2009_2019_median_MEAN")))
centers_11 = centers %>% 
  select(ends_with("MEAN"), -c(contains("_1_MEAN"), contains("_01_MEAN"), contains("2_MEAN"),
                               contains("FPAR_3"), contains("LAI_3"), contains("4_MEAN"), 
                               contains("5_MEAN"), contains("6_MEAN"), contains("FPAR_7"),
                               contains("LAI_7"), contains("8_MEAN"), contains("9_MEAN"), 
                               contains("10_MEAN"), contains("FPAR_11"), contains("LAI_11"),
                               contains("12_MEAN"), contains("norm_annual_MEAN"),
                               contains("2009_2019_median_MEAN")))

