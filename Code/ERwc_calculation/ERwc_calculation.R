# Calculates Mean ERwc per site

# Makes Figure S3

# Maggi Laan
#maggi.laan@gmail.com


# Libraries ---------------------------------------------------------------

library(readxl)
library(tidyverse)
library(reshape)
library(ggpmisc)
library(ggpubr)

rm(list=ls(all=TRUE))

current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd("../..")
getwd()

##########################################################################
## DO slope

# get data here: https://data.ess-dive.lbl.gov/view/doi%3A10.15485%2F1892052, v3_SFA_SpatialStudy_2021_SensorData.zip

# Read in miniDOT manual chamber summary file
DO_summary<-read_csv('./Data/Multiple_linear_regression/v3_Minidot_Summary_Statistics.csv', comment = '#', na = c('N/A', -9999))

# Pull out slope/NRMSE from each site 
DO_slope<-DO_summary %>%
  select(Site_ID,Dissolved_Oxygen_1_Slope,
                              Dissolved_Oxygen_2_Slope,Dissolved_Oxygen_3_Slope,
                              Dissolved_Oxygen_1_NRMSE,Dissolved_Oxygen_2_NRMSE,Dissolved_Oxygen_3_NRMSE, Dissolved_Oxygen_1_RMSE, Dissolved_Oxygen_2_RMSE, Dissolved_Oxygen_3_RMSE, Dissolved_Oxygen_1_Rsquared, Dissolved_Oxygen_2_Rsquared, Dissolved_Oxygen_3_Rsquared, Temperature_1_Mean, Temperature_2_Mean, Temperature_3_Mean) %>% 
  rowwise() %>% 
  mutate(SD_ERwc = sd(c(Dissolved_Oxygen_1_Slope, Dissolved_Oxygen_2_Slope, Dissolved_Oxygen_3_Slope), na.rm = TRUE))

# check positive slopes
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

write.csv(ERwc_mean, './Data/Multiple_linear_regression/', '/ERwc_Mean.csv')

## Plot S20R and T02 ####

s20r = read.csv("./Data/Multiple_linear_regression/v3_SFA_SpatialStudy_2021_SensorData/MinidotManualChamber/Data/v2_Minidot_S20R_Cle_Elum_2021-09-02.csv", skip = 15) %>% 
  mutate(DateTime = as.POSIXct(DateTime, format = "%Y-%m-%d %H:%M:%S")) %>% 
  mutate(Hour = as.numeric(format(DateTime, "%H"))) %>% 
  mutate(Minute = as.numeric(format(DateTime, "%M"))) %>% 
  mutate(Minute = if_else(Hour == 12, 60  + Minute, Minute))

lm_s20r_1 = lm(Dissolved_Oxygen_1 ~ Minute, data = s20r)
slope_s20r_1 = coef(lm_s20r_1)[2]
label_s20r_1 = paste0("slope = ", round(slope_s20r_1, 5))

lm_s20r_2 = lm(Dissolved_Oxygen_2 ~ Minute, data = s20r)
slope_s20r_2 = coef(lm_s20r_2)[2]
label_s20r_2 = paste0("slope = ", round(slope_s20r_2, 5))

lm_s20r_3 = lm(Dissolved_Oxygen_3 ~ Minute, data = s20r)
slope_s20r_3 = coef(lm_s20r_3)[2]
label_s20r_3 = paste0("slope = ",round(slope_s20r_3,5),"0")


t02 = read.csv("./Data/Multiple_linear_regression/v3_SFA_SpatialStudy_2021_SensorData/MinidotManualChamber/Data/v2_Minidot_T02_Yakima_2021-09-07.csv", skip = 15)%>% 
  mutate(DateTime = as.POSIXct(DateTime, format = "%Y-%m-%d %H:%M:%S")) %>% 
  mutate(Hour = as.numeric(format(DateTime, "%H"))) %>% 
  mutate(Minute = as.numeric(format(DateTime, "%M"))) %>% 
  mutate(Minute = if_else(Hour == 12, 60  + Minute, Minute))

lm_t02_1 = lm(Dissolved_Oxygen_1 ~ Minute, data = t02)
slope_t02_1 = coef(lm_t02_1)[2]
label_t02_1 = paste0("slope = ", round(slope_t02_1, 5))

lm_t02_2 = lm(Dissolved_Oxygen_2 ~ Minute, data = t02)
slope_t02_2 = coef(lm_t02_2)[2]
label_t02_2 = paste0("slope = ", round(slope_t02_2, 5))

lm_t02_3 = lm(Dissolved_Oxygen_3 ~ Minute, data = t02)
slope_t02_3 = coef(lm_t02_3)[2]
label_t02_3 = paste0("slope = ", round(slope_t02_3, 5))

s20r_do_plot = ggplot(s20r) + 
  geom_point(aes(x = DateTime, y = Dissolved_Oxygen_1), color = "red", shape = 21, size = 5) +
  geom_smooth(method = "lm", se = FALSE, aes(x = DateTime, y = Dissolved_Oxygen_1), color = "red") +
  annotate("text", x = s20r$DateTime[81], y = 9.17, label = label_s20r_1, color = "red", size = 3)+
  geom_point(aes(x = DateTime, y = Dissolved_Oxygen_2), color = "blue", shape = 21, size = 5) +
  geom_smooth(method = "lm", se = FALSE, aes(x = DateTime, y = Dissolved_Oxygen_2), color = "blue") +
  annotate("text", x = s20r$DateTime[81], y = 9.14, label = label_s20r_2, color = "blue", size = 3)+
  geom_point(aes(x = DateTime, y = Dissolved_Oxygen_3), color = "black", shape = 21, size = 5) +
  geom_smooth(method = "lm", se = FALSE, aes(x = DateTime, y = Dissolved_Oxygen_3), color = "black") +
  annotate("text", x = s20r$DateTime[81], y = 9.11, label = label_s20r_3, color = "black", size = 3)+
  ylim(8.7, 9.2)+
  scale_x_datetime(date_labels = "%m/%d %H:%M") +
  theme_bw() + 
  annotate("text", x = s20r$DateTime[10], y = 9.17, label = "Site ID S20R", size = 4) +
  labs(x = "Date and Time in 2021", y = bquote("Dissolved Oxygen (mg L"^-1*")")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10), size = 12), 
        axis.title.y = element_text(margin = margin(r = 10), size = 12), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10))

s20r_do_plot

t02_do_plot = ggplot(t02) + 
  geom_point(aes(x = DateTime, y = Dissolved_Oxygen_1), color = "red", shape = 21, size = 5) +
  geom_smooth(method = "lm", se = FALSE, aes(x = DateTime, y = Dissolved_Oxygen_1), color = "red") +
  annotate("text", x = t02$DateTime[81], y = 8.77, label = label_t02_1, color = "red", size = 3)+
  geom_point(aes(x = DateTime, y = Dissolved_Oxygen_2), color = "blue", shape = 21, size = 5) +
  geom_smooth(method = "lm", se = FALSE, aes(x = DateTime, y = Dissolved_Oxygen_2), color = "blue") +
  annotate("text", x = t02$DateTime[81], y = 8.74, label = label_t02_2, color = "blue", size = 3)+
  geom_point(aes(x = DateTime, y = Dissolved_Oxygen_3), color = "black", shape = 21, size = 5) +
  geom_smooth(method = "lm", se = FALSE, aes(x = DateTime, y = Dissolved_Oxygen_3), color = "black") +
  annotate("text", x = t02$DateTime[81], y = 8.71, label = label_t02_3, color = "black", size = 3)+
  ylim(8.3, 8.8)+
  scale_x_datetime(date_labels = "%m/%d %H:%M") +
  theme_bw() + 
  annotate("text", x = t02$DateTime[10], y = 8.77, label = "Site ID T02", size = 4) +
  labs(x = "Date and Time in 2021", y = bquote("Dissolved Oxygen (mg L"^-1*")")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10), size = 12), 
        axis.title.y = element_text(margin = margin(r = 10), size = 12), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10))

t02_do_plot

s20r_temp_plot = ggplot(s20r) + 
  geom_point(aes(x = DateTime, y = Temperature_1), color = "red", shape = 21, size = 5) +
  geom_point(aes(x = DateTime, y = Temperature_2), color = "blue", shape = 21, size = 5) +
  geom_point(aes(x = DateTime, y = Temperature_3), color = "black", shape = 21, size = 5) +
  ylim(17.4, 18.0)+
  scale_x_datetime(date_labels = "%m/%d %H:%M") +
  theme_bw() + 
  annotate("text", x = s20r$DateTime[10], y = 17.97, label = "Site ID S20R", size = 4) +
  labs(x = "Date and Time in 2021", y = bquote("Temperature ("~degree*"C)")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10), size = 12), 
        axis.title.y = element_text(margin = margin(r = 10), size = 12), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10))

s20r_temp_plot

t02_temp_plot = ggplot(t02) + 
  geom_point(aes(x = DateTime, y = Temperature_1), color = "red", shape = 21, size = 5) +
  geom_point(aes(x = DateTime, y = Temperature_2), color = "blue", shape = 21, size = 5) +
  geom_point(aes(x = DateTime, y = Temperature_3), color = "black", shape = 21, size = 5) +
  ylim(20.0, 20.6)+
  scale_x_datetime(date_labels = "%m/%d %H:%M") +
  theme_bw() + 
  annotate("text", x = t02$DateTime[10], y = 20.57, label = "Site ID T02", size = 4) +
  labs(x = "Date and Time in 2021", y = bquote("Temperature ("~degree*"C)")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10), size = 12), 
        axis.title.y = element_text(margin = margin(r = 10), size = 12), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10))

t02_temp_plot

ex_combine = ggarrange(s20r_temp_plot, t02_temp_plot, s20r_do_plot, t02_do_plot, nrow = 2, ncol = 2)

ex_combine

ggsave(file.path('./Figures',"Figure_S3_Temp_DO_Plots.png"),
       plot = ex_combine, width = 10, height = 6, dpi = 300,device = "png")
