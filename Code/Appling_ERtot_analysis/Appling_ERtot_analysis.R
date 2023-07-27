# Input data files:
# 1) daily_predictions.csv (Source: Appling et al. (2018b, c). GPP and ER units in g O2 m^-2 d^-1. Depth in m.
# NOTE: The data file (156 MB) is available in the online data repository associated with the manuscript. Read the "Code and data availability" statement in the manuscript for information about accessing the data repository. The data file can also be downloaded directly here: https://www.sciencebase.gov/catalog/item/59eb9c0ae4b0026a55ffe389.
# 2) streampulse_synthesis_statset.csv (StreamPULSE sites) (Source: B. Hall, Jr)

# Read data
library(readr)
daily_predictions <- read.csv("C:/Users/fult771/OneDrive - PNNL/Documents/GitHub/YRB_Water_Column_Respiration/Data/Appling_ERtot_analysis/daily_predictions.csv")
colnames(daily_predictions)

# Convert daily aereal respiration rates (ER; g O2 m^-2 d^-1) to volumetric rates (mg O2 L^-1 d^-1) by multiplying ER by 1/depth (m)
daily_predictions$ERvolumetric <- daily_predictions$ER * 1/daily_predictions$depth
daily_predictions$ERvolumetric

# Calculate mean values by site for all variables in daily_predictions
mean_daily_predictions_by_site <- aggregate.data.frame(daily_predictions, list(SiteID = daily_predictions$site_name), mean)
mean_daily_predictions_by_site

# Calculate the number of sites and with ERvolumetric 1) > 0 and 2) > 0.5
numSitesPositiveERtot <- sum(mean_daily_predictions_by_site$ERvolumetric > 0)
numSitesPositiveERtot  ##ERtot > 0 mg O2 L^-1 d^-1 = 15
numSitesERtotGrtThanHalfGram <- sum(mean_daily_predictions_by_site$ERvolumetric > 0.5)
numSitesERtotGrtThanHalfGram  #ERtot > 0.5 mg O2 L^-1 d^-1 = 11
#Write data file to csv for export
write.csv(mean_daily_predictions_by_site, "C:/Users/fult771/OneDrive - PNNL/Documents/GitHub/YRB_Water_Column_Respiration/Data/Appling_ERtot_analysis/mean_daily_predictions_by_site.csv")

# Remove Appling et al. (2018) sites potentially affected by process or observation error by matching site ID codes in the Appling dataset ("SiteID") with StreamPULSE site ID codes ("sitecode")
# 1) read streamPULSE data (Source: B. Hall)
streampulse_synthesis_statset <- read.csv("C:/Users/fult771/OneDrive - PNNL/Documents/GitHub/YRB_Water_Column_Respiration/Data/Appling_ERtot_analysis/streampulse_synthesis_statset.csv")
colnames(streampulse_synthesis_statset)
# 2) Subset sites from mean_daily_predictions_by_site by matching sites where mean_daily_predictions_by_site$SiteID = 20210902_streampulse_synthesis_statset$sitecode
mean_ERvolumetric_best_streamPULSEsites <- subset(mean_daily_predictions_by_site, subset =  mean_daily_predictions_by_site$SiteID %in% streampulse_synthesis_statset$sitecode)
colnames(mean_ERvolumetric_best_streamPULSEsites)
# 3) Count the number of sites where ERvolumetric > 0
numBestSitesGtrThan <- sum(mean_ERvolumetric_best_streamPULSEsites$ERvolumetric > 0) 
numBestSitesGtrThan  # Number of sites where ERvolumetric > 0 = 1 (SiteID = nwis_12100490)
# 4) Export data file to csv for kernel density plots, etc. by Xinming
write.csv(mean_ERvolumetric_best_streamPULSEsites, "C:/Users/fult771/OneDrive - PNNL/Documents/GitHub/YRB_Water_Column_Respiration/Data/Appling_ERtot_analysis/mean_ERvolumetric_best_streamPULSEsites.csv")

