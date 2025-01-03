library(tidyverse)
library(moments)
library(gtable)
library(gridExtra)
library(ggpubr)
library(readxl)

rm(list=ls());graphics.off()

current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd("../..")
getwd()

################################################
# read in ERwc data 

data = read.csv(file.path('./Data/Multiple_linear_regression/spatial_data.csv')) %>% 
  select(-c(Water_Column_Respiration, Temperature_Mean)) ### this is old .csv, has S68 while published version does not

erwc_sps = read.csv(file.path('./Data/Multiple_linear_regression/ERwc_Mean.csv')) %>% 
  select(-c(X)) %>% 
  mutate(Mean_ERwc = round(Mean_ERwc, 3)) %>% 
  mutate(Mean_Temp = round(Mean_Temp, 3)) %>% 
  left_join(data, by = "Site_ID") %>% 
  filter(Mean_ERwc < 0.5) # remove positive respiration rates > 0.5

# read in  ERtotal data

ERriv = read.csv(file.path('./Data/Appling_ERtot_analysis','mean_ERtot_cleaned.csv')) %>% select(-c(X))

# Keeps Devol min/max, 
erwc_lit = read_xlsx(file.path("C:/Users/laan208/OneDrive - PNNL/Water_Column_Respiration_Spatial/Drafts/Maggi/Table_2_Calculations.xlsx"), skip = 1,sheet = 2) %>% 
  #filter(!grepl("Reisinger", Paper)) %>% 
  filter(!grepl("Gagne", Paper)) %>%  #same as Ward 2018?
  select(c(Paper, River, `Basin/Station`, Corrected_Value)) %>% 
  rename(Water_Column_Respiration_Literature = Corrected_Value)

median(erwc_sps$Mean_ERwc) # our median: -0.579
mean(erwc_sps$Mean_ERwc) # our mean: -0.817
sd(erwc_sps$Mean_ERwc)/sqrt(length((erwc_sps$Mean_ERwc))) # our SE: 0.19

## calculate the skewness for ERtotal and ER water in this study
sk1= round(skewness(ERriv$ERtot_Volumetric),2) # - 1.66
sk2= round(skewness(erwc_sps$Mean_ERwc),2) # -3.42

##################### Density Plots ###########################

## Density Plot for ERwc + ERlit (v lines) with legend
# might need to change geom_density bounds here?

p1 <- ggplot() + 
  geom_density(data=erwc_sps, aes(x=Mean_ERwc,fill='wc'),color='blue',adjust = 4, bounds = c(min(erwc_sps$Mean_ERwc), max(erwc_sps$Mean_ERwc)))+
  geom_vline(aes(xintercept=median(erwc_sps$Mean_ERwc)),color="blue",  size=1)+
  geom_vline(data=erwc_lit, aes(xintercept=Water_Column_Respiration_Literature,color='lit'),linetype="dashed")+
  labs(x = expression("ER"[wc]*" (mg O"[2]*" L"^-1*" d"^-1*")"), y = 'Density')+
  scale_fill_manual("",breaks = c("wc"),labels = c(expression("ER"[wc]*" (this study)")),
                    values = c("skyblue"))+
  scale_colour_hue("",breaks = c("lit"),labels = c( expression("ER"[wc]*" (Lit)"))
  )+
  scale_linetype_manual("",breaks = c("lit"),labels = c( expression("ER"[wc]*"(Lit)")),
                        values = c("dashed"))+
  theme_classic()+
  #xlim(min(erwc_sps$Mean_ERwc), 0)+
  theme(legend.position = c(0.25,0.95),
        legend.justification = c("top"),
        legend.margin = margin(-12, 0, 3, 3),
        legend.text = element_text(size=8, hjust = 0, margin = margin(l = 2, r = 5, unit = "pt")),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.spacing.y = unit(0, 'cm'),
        legend.box.background = element_rect(fill = "white", linewidth = 0.4,colour = "black"),
        legend.key = element_rect(fill = "white", color = "black", linewidth = 0.5))

p1
# Extract the colour legend - leg1
leg1 <- gtable_filter(ggplot_gtable(ggplot_build(p1)), "guide-box") 

ggsave(file.path('./Plots',"hist_density_plot_ERwater_legend_ML.png"),
       plot=p1, width = 4, height = 3, dpi = 300,device = "png") #grid.arrange(p1,p2, nrow=1)


## Updated histograms/kernel density plots?

lit_sps_density = ggplot() + 
  geom_density(data=erwc_sps, aes(x=Mean_ERwc,fill='wc'),color='blue',adjust = 4, #bounds = c(min(erwc_sps$Mean_ERwc), max(erwc_sps$Mean_ERwc))
               )+
  geom_vline(aes(xintercept=median(erwc_sps$Mean_ERwc)),color="blue",  size=1)+
  geom_density(data = erwc_lit, aes(x = Water_Column_Respiration_Literature, fill = "lit"), color = "#F9847B", adjust = 4, alpha= 0.5) +
  geom_vline(aes(xintercept=median(erwc_lit$Water_Column_Respiration_Literature)),color="#F9847B",  size=1)+
  labs(x = expression("ER"[wc]*" (mg O"[2]*" L"^-1*" d"^-1*")"), y = 'Density')+
  scale_fill_manual("",breaks = c("wc", "lit"),labels = c(expression("ER"[wc]*" (this study)"), expression("ER"[wc]*" (Lit)")),
                    values = c("skyblue", "#F9847B"))+
  theme_classic()+
  theme(legend.position = c(0.25,0.95),
        legend.justification = c("top"),
        legend.margin = margin(-12, 0, 3, 3),
        legend.text = element_text(size=8, hjust = 0, margin = margin(l = 2, r = 5, unit = "pt")),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.spacing.y = unit(0, 'cm'),
        legend.box.background = element_rect(fill = "white", linewidth = 0.4,colour = "black"),
        legend.key = element_rect(fill = "white", color = "black", linewidth = 0.5))

lit_sps_density

ggsave(file.path('./Plots',"sps_lit_density_plot_ML.png"),
       plot=lit_sps_density, width = 4, height = 3, dpi = 300,device = "png")

lit_sps_hist = ggplot() + 
  geom_histogram(data = erwc_lit, aes(x = Water_Column_Respiration_Literature, fill = "lit"), color = "#F9847B", alpha = 0.7) +
  geom_vline(aes(xintercept=median(erwc_lit$Water_Column_Respiration_Literature)),color="#F9847B",  size=1)+
  geom_histogram(data=erwc_sps, aes(x=Mean_ERwc,fill='wc'),color='blue', alpha = 0.7)+
  geom_vline(aes(xintercept=median(erwc_sps$Mean_ERwc)),color="blue",  size=1)+
  #geom_vline(data=erwc_lit, aes(xintercept=Water_Column_Respiration_Literature,color='lit'),linetype="dashed")+
  labs(x = expression("ER"[wc]*" (mg O"[2]*" L"^-1*" d"^-1*")"))+
  scale_fill_manual("",breaks = c("wc", "lit"),labels = c(expression("ER"[wc]*" (this study)"), expression("ER"[wc]*" (Lit)")),
                    values = c("skyblue", "#F9847B"))+
  #scale_colour_hue("",breaks = c("lit"),labels = c( expression("ER"[wc]*" (Lit)"))
  #)+
  #("",breaks = c("lit"),labels = c( expression("ER"[wc]*"(Lit)")),
  #                      values = c("dashed"))+
  theme_classic()+
  #xlim(min(erwc_sps$Mean_ERwc), 0)+
  theme(legend.position = c(0.25,0.95),
        legend.justification = c("top"),
        legend.margin = margin(-12, 0, 3, 3),
        legend.text = element_text(size=8, hjust = 0, margin = margin(l = 2, r = 5, unit = "pt")),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.spacing.y = unit(0, 'cm'),
        legend.box.background = element_rect(fill = "white", linewidth = 0.4,colour = "black"),
        legend.key = element_rect(fill = "white", color = "black", linewidth = 0.5))

lit_sps_hist

ggsave(file.path('./Plots',"sps_lit_histograms_plot_ML.png"),
       plot=lit_sps_hist, width = 4, height = 3, dpi = 300,device = "png")





# density plot for ERtotal 
# figure out ymins/maxes for ERwc_sps and ERwc_lit
p2 <- ggplot(ERriv, aes(x=ERtot_Volumetric,color='tot',fill="tot")) + 
  geom_density(bounds = c(min(ERriv$ERtot_Volumetric), max(ERriv$ERtot_Volumetric)))+ 
  geom_vline(aes(xintercept=median(ERtot_Volumetric)), color='black', size=1)+ 
  # xlab(expression("ER"[tot]*"")) +
  # ylab('Density') + scale_fill_grey() + 
  theme_classic()+
  labs(x = expression("ER"[tot]*" (mg O"[2]*" L"^-1*" d"^-1*")"), y = 'Density', color = "Legend")+
  geom_rect(aes(xmin = min(erwc_sps$Mean_ERwc), xmax = 0,ymin=0.001,ymax=0.03,colour="wc",fill='wc'),alpha=0.1)+ #ER WC
  geom_rect(aes(xmin = min(erwc_lit$Water_Column_Respiration_Literature),xmax = max(erwc_lit$Water_Column_Respiration_Literature), ymin=0.001, ymax=0.03, colour="lit", fill='lit'), alpha = 0.1) + #ER lit
  scale_colour_manual("",breaks = c("tot", "wc", "lit"),labels = c(expression("ER"[tot]*""), expression("ER"[wc]*" range (this study)"), expression("ER"[wc]*" range (Lit) ")),
                      values = c("black", "blue", "#F9847B"),aesthetics = c("colour"))+
  scale_fill_manual("",breaks = c("tot", "wc", "lit"),labels = c(expression("ER"[tot]*""), expression("ER"[wc]*" range (this study)"), expression("ER"[wc]*" range (Lit) ")),
                    values = c("grey", "skyblue", "#F9847B"), guide = guide_legend(override.aes = list(alpha = .5)))+
  #xlim(min(ERriv$Total_Ecosystem_Respiration_Volumetric), max(ERriv$Total_Ecosystem_Respiration_Volumetric))+ # this was originally -20 to max(ERriv$ERtot), not sure what ERtot was
  theme(
    legend.position = c(.35, .95),
    legend.justification = c( "top"),
    legend.margin = margin(5, 5, 3, 5),
    legend.text = element_text(size=10,hjust = 0, margin = margin(l = 0, r = 5, unit = "pt")),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),
    legend.key = element_rect(fill = "white", color = "black", linewidth = 0.2),
    legend.box.just = "right"
  )

p2

ggsave(file.path('./Plots',"hist_density_plot_ERtot_legend_ML.png"),
       plot=p2, width = 4, height = 3, dpi = 300,device = "png") #grid.arrange(p1,p2, nrow=1)


## Density plot for ERtotal with legend
# Bounds might need to be changed here as well?

p3 <- ggplot(ERriv, aes(x=ERtot_Volumetric,color='tot',fill="tot")) + 
  geom_density(bounds = c(min(ERriv$ERtot_Volumetric), max(ERriv$ERtot_Volumetric)))+ 
  geom_vline(aes(xintercept=median(ERtot_Volumetric)), color='black', size=0.8)+ 
  # xlab(expression("ER"[tot]*"")) +
  # ylab('Density') + scale_fill_grey() + 
  theme_classic()+
  labs(x = expression("ER"[tot]*" (mg O"[2]*" L"^-1*" d"^-1*")"), y = 'Density', color = "Legend")+
  # geom_rect(aes(xmin=-4.63,xmax=-0.02,ymin=0,ymax=0.08,colour="lit",fill='lit'))+ #ER lit
  geom_rect(aes(xmin= min(erwc_sps$Mean_ERwc),xmax=max(erwc_sps$Mean_ERwc),ymin=0,ymax=0.026,colour="wc",fill='wc'), alpha=0.5)+ # changed this from stagnant values to pull from df
  geom_rect(aes(xmin = min(erwc_lit$Water_Column_Respiration_Literature),xmax= max(erwc_lit$Water_Column_Respiration_Literature),ymin=0,ymax=0.026, color = "lit", fill="lit"), alpha = 0.5)+
   scale_colour_manual("",breaks = c("tot", "wc", "lit"),labels = c(expression("ER"[tot]*""), expression("ER"[wc]*" range (this study)"), expression("ER"[wc]*" range (Lit) ")),
                      values = c("black", "blue", "#F9847B"),aesthetics = c("colour"))+
  scale_fill_manual("",breaks = c("tot", "wc", "lit"),labels = c(expression("ER"[tot]*""), expression("ER"[wc]*" range (this study)"), expression("ER"[wc]*" range (Lit) ")),
                    values = c("grey", "lightblue", "#fbb1ac"))+
  xlim(-20, 1)+#max(ERriv$ERtot_Volumetric))+
theme(
  legend.position = c(.275, .95),
  legend.justification = c( "top"),
  legend.margin = margin(-12, 0, 3, 3),
  legend.text = element_text(size=8, hjust = 0, margin = margin(l = 2, r = 5, unit = "pt")),
  legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
  legend.key = element_rect(fill = "white", color = "black", linewidth = 0.4),
  legend.box.just = "right"
)

p3

combined = grid.arrange(p1,p3, nrow=1)
combined = ggarrange(p1, p3, nrow = 1, labels = c("(a)", "(b)"), label.x = c(0.9, 0.85), label.y = c(0.95, 0.95))

combined

ggsave(file.path('./Plots',"bounded_combined_hist_density_plot_ML_pos_removed.png"), plot=combined, width = 8, height = 3, dpi = 300,device = "png") 

## Density Plot all data ####

## Density Plot Smoothed? This is weird?
all_density = ggplot() + 
  geom_density(erwc_sps, mapping = aes(x = Mean_ERwc, color = "wc", fill = "wc"), alpha = 0.75, adjust = 4) +
  geom_density(ERriv, mapping = aes(x=ERtot_Volumetric, y = ..density..*3,color='tot',fill="tot"), alpha = 0.5) +
  geom_density(erwc_lit, mapping = aes(x = Water_Column_Respiration_Literature, color = "lit",  fill = "lit"), alpha = 0.75, adjust = 4) +
  scale_y_continuous(
    name = expression("ER"[lit]*" and ER"[wc]*" Density"),
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
  labs(x = expression("ER"[tot]*" (mg O"[2]*" L"^-1*" d"^-1*")"), y = 'Density', color = "Legend")+
  theme(
    legend.position = c(.175, .95),
    legend.justification = c( "top"),
    legend.margin = margin(-12, 0, 3, 3),
    legend.text = element_text(size=8, hjust = 0, margin = margin(l = 2, r = 5, unit = "pt")),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    legend.key = element_rect(fill = "white", color = "black", linewidth = 0.4),
    legend.box.just = "right")

all_density

ggsave(file.path('./Plots',"all_density_plot_ML.png"), plot= all_density, width = 8, height = 3, dpi = 300,device = "png") 

## Histogram comparison of values

all_hist = ggplot() + 
  geom_histogram(data = ERriv, aes(x = ERtot_Volumetric, color = "tot", fill = "tot"), binwidth = 1) +
  geom_histogram(data = erwc_lit, aes(x = Water_Column_Respiration_Literature, color = "lit", fill = "lit"), binwidth = 0.5) +
  geom_histogram(data = erwc_sps, aes(x = Mean_ERwc, color = "wc", fill = "wc"), binwidth = 0.25) +
scale_colour_manual("",breaks = c("tot", "wc", "lit"),labels = c(expression("ER"[tot]*""), expression("ER"[wc]*" (this study)"), expression("ER"[wc]*" (Lit) ")),
                    values = c("black", "blue", "#F9847B"),aesthetics = c("colour"))+
  scale_fill_manual("",breaks = c("tot", "wc", "lit"),labels = c(expression("ER"[tot]*""), expression("ER"[wc]*" (this study)"), expression("ER"[wc]*" (Lit) ")),
                    values = c("grey", "lightblue", "#fbb1ac"))+
  geom_vline(ERriv, mapping = aes(xintercept=median(ERtot_Volumetric)), color='black', size=0.5)+ 
  geom_vline(erwc_lit, mapping = aes(xintercept=median(Water_Column_Respiration_Literature)), color="#F9847B", size=0.5)+ 
  geom_vline(erwc_sps, mapping = aes(xintercept=median(Mean_ERwc)), color="blue", size=0.5)+ 
  theme_classic()+
  labs(x = expression("ER"[tot]*" (mg O"[2]*" L"^-1*" d"^-1*")"), color = "Legend")+
  theme(
    legend.position = c(.175, .95),
    legend.justification = c( "top"),
    legend.margin = margin(-12, 0, 3, 3),
    legend.text = element_text(size=8, hjust = 0, margin = margin(l = 2, r = 5, unit = "pt")),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    legend.key = element_rect(fill = "white", color = "black", linewidth = 0.4),
    legend.box.just = "right")

ggsave(file.path('./Plots',"all_hist_plot_ML.png"), plot= all_hist, width = 8, height = 3, dpi = 300,device = "png") 


### Literature Comparison by Paper ####

ward = erwc_lit %>% 
  filter(grepl("Ward", Paper))

ellis = erwc_lit %>%
  filter(grepl("Ellis", Paper)) 

genzoli = erwc_lit %>% 
  filter(grepl("Genzoli", Paper))

quay = erwc_lit %>% 
  filter(grepl("Quay", Paper))

reisinger = erwc_lit %>% 
  filter(grepl("Reisinger", Paper))

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
  mutate(`Basin/Station` = "N/A") %>% 
  select(c(Paper, River, `Basin/Station`, Water_Column_Respiration_Literature)) %>% 
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
  ylab(expression("ER"[wc]*(" mg O"[2]*" L"^-1*" d"^-1)))+
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

comb_box_dens = ggarrange(all_density, box_plot, nrow = 1, widths = c(2,  1), labels = c("(a)", "(b)"), label.x = c(0.925, 0.9), label.y = c(0.95, 0.95))

comb_box_dens

ggsave(file.path('./Plots',"density_box_plot.png"), plot=comb_box_dens, width = 12, height = 4.5, dpi = 300,device = "png") 



##############################################################
# use lme4 to estimate the means among the sites
rdata <- read.csv(file.path('./Data/Multiple_linear_regression/v3_Minidot_Summary_Statistics.csv'), skip=52)

rdata <- rdata[c('Site_ID','Dissolved_Oxygen_1_Slope','Dissolved_Oxygen_2_Slope','Dissolved_Oxygen_3_Slope')]

rdata2 <- rdata %>%
  pivot_longer(
    cols = names(rdata)[2:4], 
    names_to = "replicates",
    values_to = "ERwc") %>% 
  filter(!grepl("-9999", ERwc)) %>% 
  mutate(ERwc = ERwc * 60 *24)

## within-sample variation
wsummary<-rdata2 %>%
  group_by(Site_ID) %>%
  summarise(mean = mean(ERwc),
            sd = sd(ERwc),
            var = var(ERwc))
mean(wsummary$sd) #average standard deviation by site?
mean(wsummary$var) #average variance by site?

## lme4 fitting

lfit <- lmer( ERwc ~  (1 | Site_ID), data=rdata2)  
summary(lfit)


###############################################################
# plot to visualize the value of ERwater in each Site
#odata=data[c('Site_ID','ERwc')]
odata = rdata2 %>%
  group_by(Site_ID) %>%
  mutate(mean = mean(ERwc),median = median(ERwc)) %>% 
  filter(!grepl("S38|S83", Site_ID))

##
#order by mean values
odata <- odata[order(odata$mean),]
sites<-unique(odata$Site_ID)
odata['mean_rank'] <-'rank1'
for (i in 1:length(sites)){
  odata$mean_rank[odata$Site_ID==sites[i]] = paste0('rank',sprintf("%02d", i))
}

# odata <- odata%>%arrange(median(ERwc))%>%
#   mutate(Site_ID = factor(Site_ID, levels = unique(Site_ID)))

DotPlot <- ggplot(odata, aes(x=mean_rank, y=ERwc)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "red",alpha=0.8)+
  stat_summary(fun=mean, geom="point", shape=18,size=3, color="red") +
  geom_dotplot(binaxis='y', stackdir='center',binwidth = 0.2,dotsize = 0.4) +
  ylab(expression("ER"[wc]*" (mg O"[2]*" L"^-1*"day"^-1*")")) + xlab("Site ID") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = c(0, 0))+
  scale_x_discrete(labels=sites)

DotPlotFin <- DotPlot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(DotPlotFin)

ggsave(file.path("./Plots",paste0('ERwc_dotplot_mean_rank_ML',".png")), 
       plot=DotPlotFin, width = 6, height = 3, dpi = 300,device = "png") #grid.arrange(p1,p2, nrow=1)

##
#order by median values
odata <- odata[order(odata$median),]
sites<-unique(odata$Site_ID)
odata['median_rank'] <-'rank1'
for (i in 1:length(sites)){
  odata$median_rank[odata$Site_ID==sites[i]] = paste0('rank',sprintf("%02d", i))
}

# odata <- odata%>%arrange(median(ERwc))%>%
#   mutate(Site_ID = factor(Site_ID, levels = unique(Site_ID)))

DotPlot <- ggplot(odata, aes(x=median_rank, y=ERwc)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "red",alpha=0.8)+
  stat_summary(fun=median, geom="point", shape=18,size=3, color="red") +
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.2, dotsize = 0.4) +
  ylab(expression("ER"[wc]*" (mg O"[2]*" L"^-1*"day"^-1*")")) + xlab("Site ID") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = c(0, 0))+
  scale_x_discrete(labels=sites)

DotPlotFin <- DotPlot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(DotPlotFin)

ggsave(file.path("./Plots",paste0('ERwc_dotplot_median_rank_ML',".png")), 
       plot=DotPlotFin, width = 6, height = 3, dpi = 300,device = "png") #grid.arrange(p1,p2, nrow=1)
