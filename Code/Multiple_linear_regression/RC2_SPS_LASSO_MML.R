# Code to run LASSO on SPS ERwc data 


# Libraries ---------------------------------------------------------------
library(tidyverse)
library(dplyr)
#library(corrplot)
#library(ggpubr)
#library(ggpmisc)
# library(factoextra)
# library(stringr)
# library(glmnet)
# library(magick)

rm(list=ls());graphics.off()

cube_root <- function(x) sign(x) * (abs(x))^(1/3)

current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd("../..")
getwd()


# Read in Data ------------------------------------------------------------

mean_erwc = read.csv("C:/GitHub/rcsfa-RC2-SPS-ERwc/Data/Multiple_linear_regression/ERwc_Mean.csv") %>% 
  select(-X)

## are these the only variables I want?
geo = read.csv("C:/GitHub/rcsfa-RC2-SPS-ERwc/Data/Multiple_linear_regression/v2_RCSFA_Extracted_Geospatial_Data_2023-06-21 (1).csv") %>% 
  select(c(site, streamorde, totdasqkm)) %>% 
  dplyr::rename(Site_ID = site)

# need to get DIC from VGC, check Ions from Sophia, and pull together in summary file code, add ultrameter water chemistry to this? add manta data to this? 
sample = read.csv("C:/GitHub/rcsfa-RC2-SPS-ERwc/Data/Multiple_linear_regression/Summary_Not_Cleaned.csv") %>% 
  mutate(Sample_Name = str_remove(Parent, "_MEAN")) %>%   select(-c(X, Parent, X00000_NH4_mg_per_L_as_NH4_mean, X01130_Li_mg_per_L_mean, X00653_PO4_mg_per_L_as_PO4_mean, X71856_NO2_mg_per_L_as_NO2_mean)) 
  
mapping = read.csv("C:/GitHub/rcsfa-RC2-SPS-ERwc/Data/Multiple_linear_regression/v2_SPS_Sensor_Field_Metadata.csv") %>% 
  select(c(Site_ID, Sample_Name))

om <- read.csv(file.path('./Data/OM_transformation_analysis/SPS_Total_and_Normalized_Transformations_01-03-23.csv')) 

merged_data_OM <- merge(merged_data, sdata,by=c("Sample_Name"))

# Merge Data --------------------------------------------------------------

all_data = left_join(mean_erwc, mapping, by = "Site_ID") %>% 
  left_join(geo, by = "Site_ID") %>% 
  full_join(sample, by = "Sample_Name") %>% 
  select(-c(Sample_Name))


# Transform data ----------------------------------------------------------

#decide how you want to do this, eg, cube or log transform


