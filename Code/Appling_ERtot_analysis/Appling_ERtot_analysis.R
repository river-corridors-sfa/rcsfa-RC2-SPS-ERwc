# Appling_ERtot_analysis.R
# Stephanie G Fulton
# stephgfulton@gmail.com

# Input data files:
# 1) daily_predictions_ERtot_depth.csv. This reduced dataset contains the site ID, reach-averaged daily predicted ecosystem respiration (ERtot; g O2 m^-2 d^-1), and depth (m) data for the 356 streams and rivers across the conterminous United States from Appling et al. (2018b, c). The complete data file can also be downloaded directly here: https://www.sciencebase.gov/catalog/item/59eb9c0ae4b0026a55ffe389.
# 2) StreamPULSE_bestSiteIDs.csv (StreamPULSE sites). This reduced dataset contains the site ID, site name, and data source for the 222 sites in the StreamPULSE dataset from Bernhardt et al. (2018). The complete dataset can be downloaded directly here: https://data.streampulse.org/download_bulk.

# Set working directory
current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd("../..")
getwd()

# Read data
library(readr)
library(tidyverse)

daily_predictions_ERtot_depth <- read.csv(file.path('./Data/Appling_ERtot_analysis/daily_predictions_ERtot_depth.csv'))
colnames(daily_predictions_ERtot_depth)

# Convert daily areal respiration rates (ER; g O2 m^-2 d^-1) to volumetric rates (mg O2 L^-1 d^-1) by multiplying ER by 1/depth (m)
daily_predictions_ERtot_depth$Total_Ecosystem_Respiration_Volumetric <- daily_predictions_ERtot_depth$Total_Ecosystem_Respiration_Areal * 1/daily_predictions_ERtot_depth$Depth
daily_predictions_ERtot_depth$Total_Ecosystem_Respiration_Volumetric

# Calculate mean ERtot and depth by site
mean_ERtot_depth_by_site <- aggregate.data.frame(daily_predictions_ERtot_depth, list(Site_ID = daily_predictions_ERtot_depth$Site_ID), mean)
mean_ERtot_depth_by_site

# Remove Appling et al. (2018) sites potentially affected by process or observation error by matching site ID codes in the Appling dataset ("Site_ID") with StreamPULSE site ID codes ("Site_ID")

# 1) read streamPULSE data 
StreamPULSE_bestSiteIDs <- read.csv(file.path('./Data/Appling_ERtot_analysis/StreamPULSE_bestSiteIDs.csv'))
colnames(StreamPULSE_bestSiteIDs)
# 2) Subset sites from Appling dataset by matching site IDs (i.e., mean_ERtot_depth_by_site$Site_ID = StreamPULSE_bestSiteIDs$Site_ID) and remove unnecessary data (columns 1 and 3)
mean_ERtot_bestSiteIDs <- subset(mean_ERtot_depth_by_site, subset =  mean_ERtot_depth_by_site$Site_ID %in% StreamPULSE_bestSiteIDs$Site_ID)
mean_ERtot_bestSiteIDs[,-c(1,3)]
colnames(mean_ERtot_bestSiteIDs)
# 3) Check how many sites still have positive respiration rates (i.e., Total_Ecosystem_Respiration_Volumetric > 0)
numBestSitesGtrThan <- sum(mean_ERtot_bestSiteIDs$Total_Ecosystem_Respiration_Volumetric > 0) 
numBestSitesGtrThan  # Number of sites where Total_Ecosystem_Respiration_Volumetric > 0 = 1 (Site_ID = nwis_12100490)
# 4) Export data file to csv to create kernel density plots (Fig. 6) using "RC2_spatial_study_MLR_v3.R" 
write.csv(mean_ERtot_bestSiteIDs)


# Read data
library(readr)
daily_predictions_ERtot_depth <- read.csv("./Data/Appling_ERtot_analysis/daily_predictions_ERtot_depth.csv")
colnames(daily_predictions_ERtot_depth)

# Convert daily areal respiration rates (ER; g O2 m^-2 d^-1) to volumetric rates (mg O2 L^-1 d^-1) by multiplying ER by 1/depth (m)
daily_predictions_ERtot_depth$Total_Ecosystem_Respiration_Volumetric <- daily_predictions_ERtot_depth$ER * 1/daily_predictions_ERtot_depth$Depth
daily_predictions_ERtot_depth$Total_Ecosystem_Respiration_Volumetric

# Calculate mean ERtot and depth by site
mean_ERtot_depth_by_site <- aggregate.data.frame(daily_predictions_ERtot_depth, list(Site_ID = daily_predictions_ERtot_depth$Site_ID), mean)
mean_ERtot_depth_by_site

# Remove Appling et al. (2018) sites potentially affected by process or observation error by matching site ID codes in the Appling dataset ("Site_ID") with StreamPULSE site ID codes ("Site_ID")
# 1) read streamPULSE data (Source: B. Hall)
StreamPULSE_bestSiteIDs <- read.csv("./Data/Appling_ERtot_analysis/StreamPULSE_bestSiteIDs.csv")
colnames(StreamPULSE_bestSiteIDs)
# 2) Subset sites from Appling dataset by matching site IDs (i.e., mean_ERtot_depth_by_site$Site_ID = StreamPULSE_bestSiteIDs$Site_ID) and remove unnecessary data (columns 1 and 3)
mean_ERtot_bestSiteIDs <- subset(mean_ERtot_depth_by_site, subset =  mean_ERtot_depth_by_site$Site_ID %in% StreamPULSE_bestSiteIDs$Site_ID)
mean_ERtot_bestSiteIDs[,-c(1,3)]
colnames(mean_ERtot_bestSiteIDs)
# 3) Check how many sites still have positive respiration rates (i.e., Total_Ecosystem_Respiration_Volumetric > 0)
numBestSitesGtrThan <- sum(mean_ERtot_bestSiteIDs$Total_Ecosystem_Respiration_Volumetric > 0) 
numBestSitesGtrThan  # Number of sites where Total_Ecosystem_Respiration_Volumetric > 0 = 1 (Site_ID = nwis_12100490)
# 4) Export data file to csv for kernel density plots, etc.
write.csv(mean_ERtot_bestSiteIDs, "./Data/Appling_ERtot_analysis/mean_ERtot_bestSiteIDs.csv")

# Maggi - Bernhardt compare ####

bernhardt = readRDS(url("https://raw.githubusercontent.com/streampulse/metabolism_synthesis/master/output_data/lotic_gap_filled.rds")) %>% #url to Bernhardt (2022) github file with gap filled data
  bind_rows() %>%
  mutate(Date = format(Date, format = "%Y-%m-%d")) %>% 
  filter(grepl("nwis", Site_ID)) %>% 
  select(c(Site_ID, Date, ER_filled)) 

mean_bernhardt_er = bernhardt %>% 
  group_by(Site_ID) %>% # take average ERtot by site
  summarise(ERtot_Areal = mean(ER_filled))  

# Read in Appling to get depth data

appling = read.csv("./Data/Appling_ERtot_analysis/daily_predictions.csv") %>% 
  mutate(Date = as.POSIXct(date, format = "%m/%d/%Y")) %>% 
  mutate(Date = format(Date, format = "%Y-%m-%d")) %>% 
  rename(Site_ID = site_name) %>% 
  select(c(Site_ID, Date, depth)) 

mean_appling_depth = appling %>% 
  group_by(Site_ID) %>% 
  summarise(mean_depth = mean(depth))

# Merge Bernhardt and Appling data 

ERtot = full_join(mean_bernhardt_er, mean_appling_depth, by = c("Site_ID")) %>% 
  drop_na(ERtot_Areal) %>%  # remove sites with no ER_filled, leaving 208 site
  mutate(ERtot_Volumetric = ERtot_Areal * (1/mean_depth)) # divide by mean depth to get into mg O2/L/day

write.csv(ERtot, "./Data/Appling_ERtot_analysis/mean_ERtot_cleaned.csv")
