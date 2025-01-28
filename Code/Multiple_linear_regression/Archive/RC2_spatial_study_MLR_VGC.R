# RC2 spatial study - Multiple linear regression for manual_chamber_data
# Xinming Lin Nov 22 2022
# xinming.lin@pnnl.gov

################################################
# Read in data 
################################################
rm(list=ls(all=TRUE))
library(MASS)
#library(relaimpo)
# library(visreg)
# library(ggstatsplot)
# library(caret)
# library(leaps)
# library(car)
# library(patchwork)
# library(scales)
# library(moments)
# library(plotrix)
library(ggplot2)
# library(gridExtra)
# library(gtable)
# library(grid)
# library(ggbreak)
 library(tidyverse)
# library(latex2exp)
# library(lme4)
# library("PerformanceAnalytics")

# Set working directory
current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd("../..")
getwd()

# read in ERwc(ERwater),"T_mean","StreamOrde","Total_Drainage_Area"data
data<- read.csv('//pnl/projects/SBR_SFA/RC2/04_Spatial_Study_2021/Water_Column_Respiration_Re-Workup/SPS_ERwc_corrected.csv')
names(data) <-c('Site_ID','ERwc')
# set positive ERwater to 0
#data$ERwc[data$ERwc>0]<-0 # Note that this data is volumetric

################################################
# read in  ERtotal data
ERriv <- read.csv(file.path('./Data/Appling_ERtot_analysis','mean_ERtot_bestSiteIDs.csv'))
# ERriv$Total_Ecosystem_Respiration_Volumetric[ERriv$Total_Ecosystem_Respiration_Volumetric>0]<-0
names(ERriv)[4] = 'ERvolumetric'
ERriv$VGC = ERriv$Total_Ecosystem_Respiration_Areal*(1/ERriv$Depth)
#ERwater from SS2021
ERwc<-data

# ERwater data from literture
#ERwc2 <- read.csv(file.path('./Data/Multiple_linear_regression','ERwc_combined_lit_valuesV3.csv'))
ERwc2 <- read.csv(file.path('./Data/Water_column_respiration_published','Water_column_respiration_published_values.csv'))
names(ERwc2)[4] = 'ERwc1'

library(moments); library(ggbreak)
## calculate the skewness for ERtotal and ER water in this study
sk1= round(skewness(ERriv$ERvolumetric),2)
sk2= round(skewness(ERwc$ERwc),2)

################################################
# make the density plots
## density plot for ERwater in this study
## density plot for ERwater in this study
p0 <- ggplot() + 
  geom_density(data=ERwc, aes(x=ERwc,fill='wc'),color='blue',adjust = 6)+
  geom_vline(aes(xintercept=median(ERwc$ERwc)),color="blue",  size=1)+
  geom_vline(data=ERwc2, aes(xintercept=ERwc1,color='lit'),linetype="dashed")+
  #scale_x_cut(breaks=c(-0.12), which=c(1), scales=c(0.25, 1),space = 0.2)+ theme_bw()+ 
  # xlab(expression("ER"[wc]*"")) +
  # ylab('Density')  + theme_classic()+ #+ scale_fill_grey()
  labs(x = expression("ER"[wc]*" (mg O"[2]*" L"^-1*" d"^-1*")"), y = 'Density')+
  scale_fill_manual("",breaks = c("wc"),labels = c(expression("ER"[wc]*"")),
                    values = c("skyblue"))+
  # scale_colour_manual("",breaks = c("wc"),labels = c(expression("ER"[wc]*"")),
  #                   values = c("blue"))+
  scale_colour_hue("",breaks = c("lit"),labels = c( expression("ER"[wc]*" Lit"))
  )+
  scale_linetype_manual("",breaks = c("lit"),labels = c( expression("ER"[wc]*" Lit")),
                        values = c("dashed"))+theme_classic()+
  theme(legend.position ="none")

ggsave(file.path('./Plots',"hist_density_plot_ERwater.pdf"),
       plot=p0, width = 4, height = 3, dpi = 300,device = "pdf") #grid.arrange(p1,p2, nrow=1)


## density plot for ERwater with legend
p1 <- ggplot() + 
  geom_density(data=ERwc, aes(x=ERwc,fill='wc'),color='blue',adjust = 4)+
  geom_vline(aes(xintercept=median(ERwc$ERwc)),color="blue",  size=1)+
  geom_vline(data=ERwc2, aes(xintercept=ERwc1,color='lit'),linetype="dashed")+
 # scale_x_cut(breaks=c(-0.12), which=c(1), scales=c(0.25, 1),space = 0.2)+ theme_bw()+ 
  # xlab(expression("ER"[wc]*"")) +
  # ylab('Density')  + theme_classic()+ #+ scale_fill_grey()
  labs(x = expression("ER"[wc]*" (mg O"[2]*" L"^-1*" d"^-1*")"), y = 'Density')+
  scale_fill_manual("",breaks = c("wc"),labels = c(expression("ER"[wc]*" (this study)")),
                    values = c("skyblue"))+
  scale_colour_hue("",breaks = c("lit"),labels = c( expression("ER"[wc]*" (Lit)"))
                      )+
  scale_linetype_manual("",breaks = c("lit"),labels = c( expression("ER"[wc]*"(Lit)")),
                        values = c("dashed"))+theme_classic()+
 theme(legend.position = 'right',#legend.box = "horizontal",
 legend.spacing.y = unit(0, "mm"),
 legend.text = element_text(size=12),
        #legend.background = element_rect(size=0.5, linetype="solid",colour = "black"),
        #legend.justification = c("right", "top"),
        #legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
legend.key = element_rect(fill = "white", color = "black", linewidth = 0.5),
        legend.box.just = "left")
library(gtable)
# Extract the colour legend - leg1
leg1 <- gtable_filter(ggplot_gtable(ggplot_build(p1)), "guide-box")

ggsave(file.path('./Plots',"hist_density_plot_ERwater_legend.pdf"),
       plot=p1, width = 4, height = 3, dpi = 300,device = "pdf") #grid.arrange(p1,p2, nrow=1)

# plotNew <- p0 + 
#   annotation_custom(grob = leg1, xmin = -0.075, xmax = -0.045, ymin = 20, ymax = 30)

# density plot for ERtotal 
p2 <- ggplot(ERriv, aes(x=ERvolumetric,color='tot',fill="tot")) + 
  geom_density()+ 
  geom_vline(aes(xintercept=median(ERvolumetric)), color='black', size=1)+ 
  # xlab(expression("ER"[tot]*"")) +
  # ylab('Density') + scale_fill_grey() + 
  theme_classic()+
  labs(x = expression("ER"[tot]*" (mg O"[2]*" L"^-1*" d"^-1*")"), y = 'Density', color = "Legend")+
  geom_rect(aes(xmin=min(data$ERwc),xmax=max(data$ERwc),ymin=0.001,ymax=0.03,colour="wc",fill='wc'),alpha=0.1)+ #ER WC
  geom_rect(aes(xmin=min(ERwc2$ERwc1),xmax=max(ERwc2$ERwc1),ymin=0.001,ymax=0.03,colour="lit",fill='lit'))+ #ER lit
  scale_colour_manual("",breaks = c("tot", "wc", "lit"),labels = c(expression("ER"[tot]*""), expression("ER"[wc]*" range (this study)"), expression("ER"[wc]*" range (Lit) ")),
                      values = c("black", "blue", "#F9847B"),aesthetics = c("colour"))+
  scale_fill_manual("",breaks = c("tot", "wc", "lit"),labels = c(expression("ER"[tot]*""), expression("ER"[wc]*" range (this study)"), expression("ER"[wc]*" range (Lit) ")),
                      values = c("grey", "skyblue", "#F9847B"),guide = guide_legend(override.aes = list(alpha = .5)))+
  #xlim(-20, max(ERriv$ERtot))+
  theme(
    legend.position = c(.35, .95),
    legend.justification = c( "top"),
    legend.margin = margin(5, 5, 3, 5),
    legend.text = element_text(size=10,hjust = 0, margin = margin(l = 0, r = 5, unit = "pt")),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),
    legend.key = element_rect(fill = "white", color = "black", linewidth = 0.2),
    legend.box.just = "right"
  )
ggsave(file.path('./Plots',"hist_density_plot_ERtot_legend.pdf"),
       plot=p2,  width = 4, height = 3,dpi = 300,device = "pdf") #grid.arrange(p1,p2, nrow=1)


