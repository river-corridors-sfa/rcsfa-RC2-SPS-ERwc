# Input data files:
# 1) daily_predictions_ERtot_depth.csv (Source: Appling et al. (2018b, c). GPP and ER units in g O2 m^-2 d^-1. Depth in m.
# NOTE: The data file (156 MB) is available from the online data repository associated with the manuscript. See the "Code and data availability" statement in the manuscript for information about accessing the data repository. The data file can also be downloaded directly here: https://www.sciencebase.gov/catalog/item/59eb9c0ae4b0026a55ffe389.
# 2) StreamPULSE_bestSiteIDs.csv (StreamPULSE sites) (Source: B. Hall, Jr)

# Read data
library(readr)
daily_predictions_ERtot_depth <- read.csv(file.path('./Data/Appling_ERtot_analysis/daily_predictions_ERtot_depth.csv'))
colnames(daily_predictions_ERtot_depth)

# Convert daily areal respiration rates (ER; g O2 m^-2 d^-1) to volumetric rates (mg O2 L^-1 d^-1) by multiplying ER by 1/depth (m)
daily_predictions_ERtot_depth$ERvolumetric <- daily_predictions_ERtot_depth$ER * 1/daily_predictions_ERtot_depth$depth
daily_predictions_ERtot_depth$ERvolumetric

# Calculate mean ERtot and depth by site
mean_ERtot_depth_by_site <- aggregate.data.frame(daily_predictions_ERtot_depth, list(Site_ID = daily_predictions_ERtot_depth$Site_ID), mean)
mean_ERtot_depth_by_site

# Remove Appling et al. (2018) sites potentially affected by process or observation error by matching site ID codes in the Appling dataset ("Site_ID") with StreamPULSE site ID codes ("Site_ID")

# 1) read streamPULSE data (Source: B. Hall)
StreamPULSE_bestSiteIDs <- read.csv(file.path('./Data/Appling_ERtot_analysis/StreamPULSE_bestSiteIDs.csv'))
colnames(StreamPULSE_bestSiteIDs)
# 2) Subset sites from Appling dataset by matching site IDs (i.e., mean_ERtot_depth_by_site$Site_ID = StreamPULSE_bestSiteIDs$Site_ID) and remove unnecessary data (columns 1 and 3)
mean_ERtot_bestSiteIDs <- subset(mean_ERtot_depth_by_site, subset =  mean_ERtot_depth_by_site$Site_ID %in% StreamPULSE_bestSiteIDs$Site_ID)
mean_ERtot_bestSiteIDs[,-c(1,3)]
colnames(mean_ERtot_bestSiteIDs)
# 3) Check how many sites still have positive respiration rates (i.e., ERvolumetric > 0)
numBestSitesGtrThan <- sum(mean_ERtot_bestSiteIDs$ERvolumetric > 0) 
numBestSitesGtrThan  # Number of sites where ERvolumetric > 0 = 1 (Site_ID = nwis_12100490)
# 4) Export data file to csv to create kernel density plots (Fig. 6) using "RC2_spatial_study_MLR_v3.R" 
write.csv(mean_ERtot_bestSiteIDs)


# Read data
library(readr)
daily_predictions_ERtot_depth <- read.csv("C:/Users/fult771/OneDrive - PNNL/Documents/GitHub/YRB_Water_Column_Respiration/Data/Appling_ERtot_analysis/daily_predictions_ERtot_depth.csv")
colnames(daily_predictions_ERtot_depth)

# Convert daily areal respiration rates (ER; g O2 m^-2 d^-1) to volumetric rates (mg O2 L^-1 d^-1) by multiplying ER by 1/depth (m)
daily_predictions_ERtot_depth$ERvolumetric <- daily_predictions_ERtot_depth$ER * 1/daily_predictions_ERtot_depth$depth
daily_predictions_ERtot_depth$ERvolumetric

# Calculate mean ERtot and depth by site
mean_ERtot_depth_by_site <- aggregate.data.frame(daily_predictions_ERtot_depth, list(Site_ID = daily_predictions_ERtot_depth$Site_ID), mean)
mean_ERtot_depth_by_site

# Remove Appling et al. (2018) sites potentially affected by process or observation error by matching site ID codes in the Appling dataset ("Site_ID") with StreamPULSE site ID codes ("Site_ID")
# 1) read streamPULSE data (Source: B. Hall)
StreamPULSE_bestSiteIDs <- read.csv("C:/Users/fult771/OneDrive - PNNL/Documents/GitHub/YRB_Water_Column_Respiration/Data/Appling_ERtot_analysis/StreamPULSE_bestSiteIDs.csv")
colnames(StreamPULSE_bestSiteIDs)
# 2) Subset sites from Appling dataset by matching site IDs (i.e., mean_ERtot_depth_by_site$Site_ID = StreamPULSE_bestSiteIDs$Site_ID) and remove unnecessary data (columns 1 and 3)
mean_ERtot_bestSiteIDs <- subset(mean_ERtot_depth_by_site, subset =  mean_ERtot_depth_by_site$Site_ID %in% StreamPULSE_bestSiteIDs$Site_ID)
mean_ERtot_bestSiteIDs[,-c(1,3)]
colnames(mean_ERtot_bestSiteIDs)
# 3) Check how many sites still have positive respiration rates (i.e., ERvolumetric > 0)
numBestSitesGtrThan <- sum(mean_ERtot_bestSiteIDs$ERvolumetric > 0) 
numBestSitesGtrThan  # Number of sites where ERvolumetric > 0 = 1 (Site_ID = nwis_12100490)
# 4) Export data file to csv for kernel density plots, etc.
write.csv(mean_ERtot_bestSiteIDs, "C:/Users/fult771/OneDrive - PNNL/Documents/GitHub/YRB_Water_Column_Respiration/Data/Appling_ERtot_analysis/mean_ERtot_bestSiteIDs.csv")

