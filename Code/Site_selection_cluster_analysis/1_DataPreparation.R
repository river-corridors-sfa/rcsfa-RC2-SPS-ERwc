

# Load Packages:
library(tidyverse)


# Define Input Data Path:
input_path = "//pnl/projects/ColumbiaGIS/Tables/zonal/_R_/ClusterAnalysis_Files_Scripts/R_Inputs/"


# Read in Normalized Data and Centers:
Data_norm = read.csv(paste0(input_path, "AllDataNormalized_zonal.csv"))
#centers = read.csv(paste0(input_path, "5clusters_centers.csv"))
centers = read.csv(paste0(input_path, "6clusters_centers.csv"))



#---------------------------------------------------------------------
#----------------------- Create Data Subset: -------------------------
#---------------------------------------------------------------------



# SET 0 --------------------------------------------------------------

## Set 0: All statistical moments for all variables
Data_norm_0 = Data_norm %>% 
  select(-AreaSqKM)
centers_0 = centers



# SET 1 --------------------------------------------------------------

## Set 1: Using only MEAN moments for all variables
Data_norm_1 = Data_norm %>% 
  select(catchment_ID, ends_with("MEAN"))
centers_1 = centers %>% 
  select(ends_with("MEAN"))



# SET 2 --------------------------------------------------------------

## Set 2: Using only MEAN moments for all variables, 
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

## Set 3: Using only MEAN moments for all variables,
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

## Set 4: Using only MEAN moments for all variables,
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

## Set 5: Using only MEAN moments for all variables,
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

## Set 6: Using only MEAN moments for all variables,
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

## Set 7: Using only STD moments for all variables
Data_norm_7 = Data_norm %>% 
  select(catchment_ID, ends_with("STD"))
centers_7 = centers %>% 
  select(ends_with("STD"))



# SET 8 --------------------------------------------------------------

## Set 8: Using only MEAN and STD moments for all variables
Data_norm_8 = Data_norm %>% 
  select(catchment_ID, ends_with("MEAN"), ends_with("STD"))
centers_8 = centers %>% 
  select(ends_with("MEAN"), ends_with("STD"))



# SET 9 --------------------------------------------------------------

## Set 9: Using only MEAN moments for all variables, 
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

## Set 10: Using only MEAN moments for all variables, 
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

## Set 11: Using only MEAN moments for all variables, 
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

