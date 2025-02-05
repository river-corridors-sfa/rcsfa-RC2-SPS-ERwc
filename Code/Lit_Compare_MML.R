
# Code to run ERlit comparison

# Makes Figures 2, S4

# Author: Maggi Laan (maggi.laan@gmail.com)


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(moments)
library(gtable)
library(gridExtra)
library(ggpubr)
library(readxl)

rm(list=ls());graphics.off()

current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd("./..")
getwd()


# Read in ERwc Data -------------------------------------------------------

erwc_sps = read.csv(file.path('./Data/ERwc_Mean.csv')) %>% 
  mutate(Mean_ERwc = round(Mean_ERwc, 3)) %>% 
  mutate(Mean_Temp = round(Mean_Temp, 3)) %>% 
  filter(Mean_ERwc < 0.5) # remove positive respiration rates > 0.5

# read in  ERtotal data
ERriv = read.csv(file.path('./Data/','mean_ERtot_cleaned.csv'))

# Keeps Devol Min/Max 
#download .xlsx file from SI and add sheet 2 as .csv to Published_Data folder
erwc_lit = read.csv(file.path("./Data/Published_Data/","Table_S5_Calculations.csv"), skip = 1) %>% 
  select(c(Paper, River, `Basin.Station`, Corrected_Value)) %>% 
  rename(Water_Column_Respiration_Literature = Corrected_Value)

median(erwc_sps$Mean_ERwc) # our median: -0.579
mean(erwc_sps$Mean_ERwc) # our mean: -0.817
sd(erwc_sps$Mean_ERwc)
sd(erwc_sps$Mean_ERwc)/sqrt(length((erwc_sps$Mean_ERwc))) # our SE: 0.19

# Density Plots -----------------------------------------------------------

## FIGURE 2a - Density Plot all data ####

## Density Plot Smoothed
all_density = ggplot() + 
  geom_density(erwc_sps, mapping = aes(x = Mean_ERwc, color = "wc", fill = "wc"), alpha = 0.75, adjust = 4) +
  geom_density(ERriv, mapping = aes(x=ERtot_Volumetric, y = ..density..*3,color='tot',fill="tot"), alpha = 0.5) +
  geom_density(erwc_lit, mapping = aes(x = Water_Column_Respiration_Literature, color = "lit",  fill = "lit"), alpha = 0.75, adjust = 4) +
  scale_y_continuous(
    name = expression("ER"[wc]*" (Lit) and ER"[wc]*" (this study) Density"),
    sec.axis = sec_axis(~./3, name = expression("ER"[tot]*" Density"))
  ) +
  scale_colour_manual("",breaks = c("tot", "wc", "lit"),labels = c(expression("ER"[tot]*""), expression("ER"[wc]*" (this study)"), expression("ER"[wc]*" (Lit) ")),
                      values = c("black", "blue", "#F9847B"),aesthetics = c("colour"))+
  scale_fill_manual("",breaks = c("tot", "wc", "lit"),labels = c(expression("ER"[tot]*""), expression("ER"[wc]*" (this study)"), expression("ER"[wc]*" (Lit) ")),
                    values = c("grey", "lightblue", "#fbb1ac"))+
  geom_vline(ERriv, mapping = aes(xintercept=median(ERtot_Volumetric)), color='black', size=0.8)+ 
  geom_vline(erwc_lit, mapping = aes(xintercept = median(Water_Column_Respiration_Literature)), color = "#F9847B", size = 0.8) +
  geom_vline(erwc_sps, mapping = aes(xintercept = median(Mean_ERwc)), color = "blue", size = 0.8) +
  theme_classic()+
  labs(x = expression("ER"[tot]*" and ER"[wc]*" (mg O"[2]*" L"^-1*" d"^-1*")"), y = 'Density', color = "Legend")+
  theme(
    legend.position = c(.175, .95),
    legend.justification = c( "top"),
    legend.margin = margin(-12, 0, 3, 3),
    legend.text = element_text(size=8, hjust = 0, margin = margin(l = 2, r = 5, unit = "pt")),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    legend.key = element_rect(fill = "white", color = "black", linewidth = 0.4),
    legend.box.just = "right", 
    axis.title.y = element_text(size = 10))

all_density

#ggsave(file.path('./Figures',"Figure2a_ERlit_Density.png"), plot= all_density, width = 8, height = 3, dpi = 300,device = "png") 

## FIGURE 2b - Literature Comparison by Paper ####

median_lit = erwc_lit %>% 
  group_by(Paper) %>% 
  summarise(median_lit = median(Water_Column_Respiration_Literature), 
            min_lit = max(Water_Column_Respiration_Literature), 
            max_lit = min(Water_Column_Respiration_Literature))

# Boxplots of Medians -----------------------------------------------------

box_df = erwc_sps %>% 
  rename(Water_Column_Respiration_Literature = Mean_ERwc) %>% 
  mutate(Paper = "Yakima River basin") %>% 
  mutate(River = "N/A") %>% 
  mutate(`Basin.Station` = "N/A") %>% 
  select(c(Paper, River, `Basin.Station`, Water_Column_Respiration_Literature)) %>% 
  rbind(erwc_lit) %>% 
  mutate(group = ifelse(Paper == "Yakima River basin", "YRB", "Lit"))

sps_median = median(erwc_sps$Mean_ERwc)
lit_median = median(erwc_lit$Water_Column_Respiration_Literature)


box_plot = ggplot() +
  geom_boxplot(box_df, mapping = aes(x = Paper, y = Water_Column_Respiration_Literature, fill = group, color = group)) +
  geom_hline(aes(yintercept = sps_median, color = "YRB"), size = 0.5, linetype = "dashed") + 
  geom_hline(linetype = "dashed", size = 0.5, aes(yintercept = lit_median, color = "Lit")) +
  theme_bw () +
  xlab("") +
  scale_color_manual("", breaks = c("YRB", "Lit"), labels = c(expression("ER"[wc]*" (this study)"), expression("ER"[wc]*" (Lit) ")),
           values = c("blue", "#F9847B"))+
  scale_fill_manual("", breaks = c("YRB", "Lit"), labels = c(expression("ER"[wc]*" (this study)"), expression("ER"[wc]*" (Lit) ")),
                    values = c("lightblue", "#fbb1ac"))+
  ylab(expression("ER"[wc]*" (mg O"[2]*" L"^-1*" d"^-1*")"))+
  #scale_color_manual(values = c("black", "#FF61CC"), 
                    # labels = c("All Literature Median", "Yakima River basin Median"))+
  theme (axis.text.x = element_text(angle = 90, vjust = 0.5), 
         legend.position = c(.225, .25),
         legend.justification = c( "top"),
         legend.margin = margin(-12, 0, 3, 3),
         legend.text = element_text(size=8, hjust = 0, margin = margin(l = 2, r = 5, unit = "pt")),
         legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
         #legend.key = element_rect(fill = "white", color = "black", linewidth = 0.4),
         legend.box.just = "right", 
         plot.margin = unit (c(1, 2, 1, 1), "lines"))


box_plot

### FIGURE 2a + 2b ####

comb_box_dens = ggarrange(all_density, box_plot, nrow = 1, widths = c(2,  1), labels = c("(a)", "(b)"), label.x = c(0.925, 0.9), label.y = c(0.95, 0.95))

comb_box_dens

ggsave(file.path('./Plots',"density_box_plot.png"), plot=comb_box_dens, width = 12, height = 4.5, dpi = 300,device = "png") 


## FIGURE S4 ####

# read in raw data
raw_data <- read.csv(file.path('./Data/Published_Data/v3_SFA_SpatialStudy_2021_SensorData/MinidotManualChamber/Plots_and_Summary_Statistics/v3_Minidot_Summary_Statistics.csv'), skip=52) %>% 
  select(c("Site_ID", "Dissolved_Oxygen_1_Slope", "Dissolved_Oxygen_2_Slope", "Dissolved_Oxygen_3_Slope")) %>%
  pivot_longer(
    cols = !Site_ID, 
    names_to = "replicates",
    values_to = "ERwc") %>% 
  filter(!grepl("-9999", ERwc)) %>% 
  mutate(ERwc = ERwc * 60 *24)

# clean data

clean_data = raw_data %>%
  group_by(Site_ID) %>%
  mutate(mean = mean(ERwc),median = median(ERwc)) %>% 
  filter(!grepl("S38|S83", Site_ID)) %>% # remove positive values > 0.5
  arrange(., mean) 

DotPlot <- ggplot(clean_data, aes(x= reorder(Site_ID, ERwc), y=ERwc)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "red",alpha=0.8)+
  stat_summary(fun=mean, geom="point", shape=18,size=3, color="red") +
  geom_dotplot(binaxis='y', stackdir='center',binwidth = 0.2,dotsize = 0.4) +
  ylab(expression("ER"[wc]*" (mg O"[2]*" L"^-1*"day"^-1*")")) + xlab("Site ID") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = c(0, 0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black"))

DotPlot

ggsave(file.path("./Figures",paste0('FigureS4_Dot_Plot_Mean_Rank',".png")), 
       plot=DotPlotFin, width = 6, height = 3, dpi = 300,device = "png") #grid.arrange(p1,p2, nrow=1)
