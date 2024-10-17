library(tidyverse)
library(dplyr)
library(corrplot)
#library(ggpubr)
#library(ggpmisc)
# library(factoextra)
# library(stringr)
library(glmnet)
# library(magick)

# Working Directory -------------------------------------------------------
current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
getwd()

rm(list=ls());graphics.off()

ultra = read.csv("./Multiple_linear_regression/v3_SFA_SpatialStudy_2021_SensorData/v2_SPS_Ultrameter_WaterChem.csv", skip = 18) %>% 
  select(c(Site_ID, Ultrameter_Specific_Conductance_1, Ultrameter_Specific_Conductance_2, Ultrameter_Specific_Conductance_3)) %>% 
  rowwise() %>% 
  mutate(mean_spc = mean(c_across(c("Ultrameter_Specific_Conductance_1", "Ultrameter_Specific_Conductance_2", "Ultrameter_Specific_Conductance_3")), na.rm = T)) %>% 
  select(c(Site_ID, mean_spc))

manta = read.csv("./Multiple_linear_regression/v3_SFA_SpatialStudy_2021_SensorData/MantaRiver/Plots_and_Summary_Statistics/v2_Manta_Summary_Statistics.csv", skip = 16) %>% 
  select(c(Site_ID, Specific_Conductance_Mean))

mapping = read.csv("./Multiple_linear_regression/v2_SPS_Sensor_Field_Metadata.csv") %>% 
  select(c(Site_ID, Sample_Name))

all_spc = merge(mapping, manta) %>% 
  left_join(ultra) %>% 
  mutate(spc = ifelse(Specific_Conductance_Mean == "-9999", mean_spc, Specific_Conductance_Mean)) %>% 
  select(c(Sample_Name, spc)) %>% 
  arrange(spc)

write.csv(all_spc, "Ions_SpC_order.csv", row.names = F)

