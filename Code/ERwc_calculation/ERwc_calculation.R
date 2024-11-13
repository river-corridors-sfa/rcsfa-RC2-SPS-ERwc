# RC2 spatial study 2021
# calculate ERwc
# M Laan July 19 2024
rm(list=ls(all=TRUE))

library(readxl)
library(tidyverse)
library(reshape)
library(ggpmisc)

DO_Path = 'C:/GitHub/rcsfa-RC2-SPS-ERwc/Data/Multiple_linear_regression/v3_Minidot_Summary_Statistics.csv'

Out_Path = 'C:/GitHub/rcsfa-RC2-SPS-ERwc/Data/Multiple_linear_regression/'

#sdir<- 'C:/Brieanne/GitHub/SSS_metabolism/'


##########################################################################
## DO slope

# Read in miniDOT manual chamber summary file
DO_summary<-read_csv(DO_Path, comment = '#', na = c('N/A', -9999))


# Pull out slope/NRMSE from each site 
DO_slope<-DO_summary %>%
  select(Site_ID,Dissolved_Oxygen_1_Slope,
                              Dissolved_Oxygen_2_Slope,Dissolved_Oxygen_3_Slope,
                              Dissolved_Oxygen_1_NRMSE,Dissolved_Oxygen_2_NRMSE,Dissolved_Oxygen_3_NRMSE, Dissolved_Oxygen_1_RMSE, Dissolved_Oxygen_2_RMSE, Dissolved_Oxygen_3_RMSE, Dissolved_Oxygen_1_Rsquared, Dissolved_Oxygen_2_Rsquared, Dissolved_Oxygen_3_Rsquared, Temperature_1_Mean, Temperature_2_Mean, Temperature_3_Mean)

DO_pos = DO_slope %>% 
  select(c(Site_ID, Dissolved_Oxygen_1_Slope, Dissolved_Oxygen_2_Slope, Dissolved_Oxygen_3_Slope)) %>% 
  reshape2::melt(id.vars = "Site_ID") %>% 
  filter(value > 0 ) %>% 
  distinct(Site_ID)
# If NRMSE > 0.03, remove slope (largest is 0.003 in this dataset, which is a 0.3% error)

for (i in c(1:3)){
  k = which(DO_slope[,paste0('Dissolved_Oxygen_',i,'_NRMSE')]>0.03)
  DO_slope[k,paste0('Dissolved_Oxygen_',i,'_Slope')]<-NA
}

# No depth measurements, so will stay as mg O2/L/min

# Look at other possible NRMSE issues 
# df_melt <- reshape2::melt(DO_slope, id.vars = "Site_ID")
# 
# df_melt <- df_melt %>%
#     mutate(Rep = sub("Dissolved_Oxygen_(\\d)_(.*)", "\\1", variable),
#            Measure = sub("Dissolved_Oxygen_(\\d)_(.*)", "\\2", variable)) %>%
#     select(-variable) %>%
#     pivot_wider(names_from = Measure, values_from = value)
# 
# slope_hist <- ggplot(data = df_melt, aes(x = Slope ))+
#     geom_histogram()
# 
# nrmse_hist <- ggplot(data = df_melt, aes(x = NRMSE ))+
#     geom_histogram()
# 
# r2_hist <- ggplot(data = df_melt, aes(x = Rsquared ))+
#     geom_histogram()
# 
# 
# scatter_slope_r2 <- ggplot(data = df_melt, aes(x = log10(abs(Slope)), y = Rsquared)) +
#     geom_point() +
#     geom_smooth(method = "lm", formula = y ~ x)+
#     stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
#                  formula = y ~ x, parse = TRUE, size = 5)
# 
# scatter_slope_nrmse_log <- ggplot(data = df_melt, aes(x = log10(abs(Slope)), y = NRMSE)) +
#     geom_point() +
#     geom_smooth(method = "lm", formula = y ~ x)+
#     stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
#                  formula = y ~ x, parse = TRUE, size = 5)
# 
#   scatter_slope_nrmse <- ggplot(data = df_melt, aes(x = Slope, y = NRMSE)) +
#     geom_point() +
#     geom_smooth(method = "lm", formula = y ~ x)+
#     stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
#                  formula = y ~ x, parse = TRUE, size = 5)
# 
# scatter_nrmse_r2 <- ggplot(data = df_melt, aes(x = Rsquared, y = NRMSE)) +
#   geom_point() +
#   geom_smooth(method = "lm", formula = y ~ x)+
#   stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
#                formula = y ~ x, parse = TRUE, size = 5)

## Convert to mg O2/L/day

ERwc <- DO_slope %>%
  mutate(Water_Column_Respiration_1 = signif(Dissolved_Oxygen_1_Slope*(60*24), 3),
         Water_Column_Respiration_2 = signif(Dissolved_Oxygen_2_Slope*(60*24), 3),
         Water_Column_Respiration_3 = signif(Dissolved_Oxygen_3_Slope*(60*24), 3)) %>%
  select(-contains('Slope')) %>% 
  rowwise() %>% 
  mutate(Mean_ERwc = mean(c(Water_Column_Respiration_1, Water_Column_Respiration_2, Water_Column_Respiration_3), na.rm = TRUE)) %>% 
  mutate(SD_ERwc = sd(c(Water_Column_Respiration_1, Water_Column_Respiration_2, Water_Column_Respiration_3), na.rm = TRUE)) %>% 
  mutate(CV_ERwc = abs((SD_ERwc/Mean_ERwc)*100))%>% 
  mutate(range_ERWc = max(c_across(Water_Column_Respiration_1:Water_Column_Respiration_3)) - min(c_across(Water_Column_Respiration_1:Water_Column_Respiration_3))) %>%
  mutate(Mean_Temp = mean(c(Temperature_1_Mean, Temperature_2_Mean, Temperature_3_Mean), na.rm = TRUE)) %>% 
  mutate(SD_Temp = sd(c(Temperature_1_Mean, Temperature_2_Mean, Temperature_3_Mean), na.rm = TRUE)) %>% 
  mutate(CV_Temp = (SD_Temp/Mean_Temp)*100) %>% 
  relocate(Mean_ERwc, .after = Site_ID) %>% 
  relocate(CV_ERwc, .after = Mean_ERwc) %>% 
    relocate(range_ERWc, .after = CV_ERwc)

ERwc_mean_melt =  reshape2::melt(ERwc, id.vars = c("Site_ID", "CV_ERwc")) %>% 
  filter(grepl("Respiration", variable)) %>% 
  arrange(CV_ERwc) 

ordered_site_ids <- ERwc_mean_melt %>%
  dplyr::select(c(Site_ID, CV_ERwc)) %>%
  distinct() %>%
  arrange(CV_ERwc) %>%
  pull(Site_ID)

ggplot(ERwc_mean_melt, aes(x = factor(Site_ID, levels = ordered_site_ids), y = value)) + 
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(ERwc, aes(x = arrange(CV_ERwc))) + 
  geom_boxplot(y = )

## pull out mean temperature

ERwc_mean = ERwc %>% 
  select(c(Site_ID, Mean_ERwc, Mean_Temp))

write.csv(ERwc_mean, paste0(Out_Path, '/ERwc_Mean.csv'))

## Plot S20R and T02

s20r = read.csv("C:/GitHub/rcsfa-RC2-SPS-ERwc/Data/Multiple_linear_regression/v3_SFA_SpatialStudy_2021_SensorData/MinidotManualChamber/Data/v2_Minidot_S20R_Cle_Elum_2021-09-02.csv", skip = 15)
  