# Reads in gap-filled, cleaned ERtot data from Bernhardt et al (2022) and depth data from Appling 2018

# Maggi Laan 
# maggi.laan@gmail.com

# Libraries ---------------------------------------------------------------
library(readr)
library(tidyverse)

rm(list=ls());graphics.off()

# Set working directory
current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd("./..")
getwd()

# Bernhardt compare ####

# get data here: https://github.com/streampulse/metabolism_synthesis
# url to Bernhardt (2022) github file with gap filled data
bernhardt = readRDS(url("https://raw.githubusercontent.com/streampulse/metabolism_synthesis/master/output_data/lotic_gap_filled.rds")) %>% 
  bind_rows() %>%
  mutate(Date = format(Date, format = "%Y-%m-%d")) %>% 
  filter(grepl("nwis", Site_ID)) %>% 
  select(c(Site_ID, Date, ER_filled)) 

mean_bernhardt_er = bernhardt %>% 
  group_by(Site_ID) %>% # take average ERtot by site
  summarise(ERtot_Areal = mean(ER_filled))  

# Read in Appling to get depth data
# get data here and add to Published_Data folder: https://www.sciencebase.gov/catalog/item/59eb9c0ae4b0026a55ffe389
appling = read.csv("./Data/Published_Data/daily_predictions.csv") %>% 
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

write.csv(ERtot, "./Data/mean_ERtot_cleaned.csv", row.names = F)
