# RC2 spatial study - Multiple linear regression for manual_chamber_data
# X Lin Nov 22 2022
################################################
# Read in data 
################################################
rm(list=ls(all=TRUE))
library(MASS)
library(relaimpo)
library(visreg)
library(ggstatsplot)
library(caret)
library(leaps)
library(car)
library(patchwork)
library(scales)
library(moments)
library(plotrix)
library(ggplot2)
library(gridExtra)
library(gtable)
library(grid)
library(ggbreak)
library("PerformanceAnalytics")
# read in DO_slope(ERwater),"T_mean","StreamOrde","TOT_BASIN_AREA"data
data<- read.csv(file.path('./Data/spatial_data.csv'))
names(data)[c(1,2,7)]<-c('Site_ID','Parent_ID','DO_slope')
data<-data[c(1,2,7:10)]
# set positive ERwater to 0
data$DO_slope[data$DO_slope>0]<-0

## Transformatio data
sdata<- read.csv(file.path('./Data','SPS_Total_and_Normalized_Transformations_01-03-23.csv'))
names(sdata)<-c('Parent_ID','Transformations','Peaks','Normalized_Transformations')
data <- merge(data,sdata,by=c("Parent_ID"))

## chemical data from 'v2_SFA_SpatialStudy_2021_Sample_Based_Surface_Water_DataPackage'
chemdata <- read.csv(file.path('./Data/2021_spatial_study_data/v2_SFA_SpatialStudy_2021_Sample_Based_Surface_Water_DataPackage/data','v2_SPS_NPOC_TN_DIC_TSS_Ions_Summary.csv'),skip=2)
chemdata <-chemdata[grep('SPS',chemdata$Sample_Name),]
names(chemdata)
chemdata <-chemdata[,c(2,4,18:20)]; 
names(chemdata)<-c('Parent_ID','DIC','NPOC','TN','TSS')
chemdata[chemdata=='-9999'] =NA
chemdata[c('DIC','NPOC','TN','TSS')] <- sapply(chemdata[c('DIC','NPOC','TN','TSS')],as.numeric)

# merge all data
data <-merge(data,chemdata,by=c("Parent_ID"))
# remove data point with NA
cdata <- na.omit(data)
###############################################################
## plot correlation matrix
vars <- c('DIC','NPOC', 'TN','TSS','T_mean','TOT_BASIN_AREA','StreamOrde','Normalized_Transformations')
png(file.path('.',paste0('exploratory_variables_correlation_matrix',".png")),
    width = 10, height = 10, units = 'in', res = 600)
#par(mfrow=c(2,2)) 
chart.Correlation(data[vars], histogram=TRUE, pch=19)
dev.off()

###############################################################
# plot to visualize the value of ERwater in each Site
DotPlot <- ggplot(data, aes(x=Site_ID, y=DO_slope)) + 
  stat_summary(fun=mean, geom="point", shape=18,size=3, color="red") +
  geom_dotplot(binaxis='y', stackdir='center',binwidth = 0.002,) +
  ylab(expression("ER"[wc]*" (mg O"[2]*" L"^-1*"day"^-1*")")) + xlab("Site ID") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = c(0, 0))
DotPlotFin <- DotPlot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(DotPlotFin)

################################################
# Stepwise Regression

#define intercept-only model
intercept_only <- lm(DO_slope ~ 1, data=cdata)

#define model with all predictors
all <- lm(DO_slope ~ DIC + NPOC + TN+TSS+T_mean+TOT_BASIN_AREA+StreamOrde+Normalized_Transformations, data = cdata)

###################################
#perform forward stepwise regression
forward <- step(intercept_only, direction='forward', scope=formula(all), trace=1)
forward$anova
forward$coefficients
###################################
#perform forward stepwise regression
backward <- step(all, direction='backward', scope=formula(all), trace=0)
backward$anova
backward$coefficients
#############
#  lm fitting using selected variables from forward stepwise selection
bfit<- lm(DO_slope ~TN +TOT_BASIN_AREA+ T_mean, data = data)
#bfit<- lm(DO_slope ~TN +TOT_BASIN_AREA+ T_mean+StreamOrde+Transformations, data = data)
summary(bfit)

## partial-regression plot
png(file.path('.',paste0('partial_regression_avPlots_3v_new',".png")),
    width = 3, height = 7, units = 'in', res = 600)
par(mfrow=c(3,1),mgp=c(2,1,0),mar=c(3.4,4.1,1,1.5))
#
avPlots(bfit, ~ TN,id=FALSE,main='',xlab='', #xlab=expression("Residuals - TDN (mg L"^-1*" )"),
        ylab=expression(paste("Residuals - ","ER"[wc]*" (mg O"[2]*" L"^-1*" day"^-1*")")))
mtext(expression("Residuals - TDN (mg L"^-1*" )"),side=1, line=2.5,cex =0.7)
avPlots(bfit, ~ TOT_BASIN_AREA,id=FALSE,main='',xlab='', #xlab=expression("Residuals - Drainage Area (km"^2*")"),
        ylab=expression(paste("Residuals - ","ER"[wc]*" (mg O"[2]*" L"^-1*" day"^-1*")")))
mtext(expression("Residuals - Drainage Area (km"^2*")"),side=1, line=2.5,cex =0.7)
avPlots(bfit, ~ T_mean,id=FALSE,main='', xlab='',#xlab=expression("Residuals - Temperature (Â°C)"),
        ylab=expression(paste("Residuals - ","ER"[wc]*" (mg O"[2]*" L"^-1*" day"^-1*")")))
mtext(expression("Residuals - Temperature (°C)"),side=1, line=2.5,cex =0.7)
dev.off()
# 

# partial-residual plots
png(file.path('.',paste0('partial_residual_crPlots_3v_new',".png")),
    width = 3, height = 8, units = 'in', res = 600)
par(mfrow=c(3,1),mgp=c(2,1,0),mar=c(3.4,4.1,1,1.5))
#par(mfrow=c(1,3))
crPlots(bfit, ~ TN,id=FALSE,main='',smooth=FALSE,xlab='',#xlab=expression("TDN (mg L"^-1*" )"),
        ylab=expression(paste("Partial Residuals - ","ER"[wc]*" (mg O"[2]*" L"^-1*" day"^-1*")")))
mtext(expression("TDN (mg L"^-1*" )"),side=1, line=2.5,cex =0.7)
crPlots(bfit, ~ TOT_BASIN_AREA,id=FALSE,main='',smooth=FALSE,xlab='', #xlab=expression("Drainage Area (km"^2*")"),
        ylab=expression(paste("Partial Residuals - ","ER"[wc]*" (mg O"[2]*" L"^-1*" day"^-1*")")))
mtext(expression("Drainage Area (km"^2*")"),side=1, line=2.5,cex =0.7)
crPlots(bfit, ~ T_mean,id=FALSE,main='',smooth=FALSE,xlab='',#xlab=expression("Temperature (?C)"),
        ylab=expression(paste("Partial Residuals - ","ER"[wc]*" (mg O"[2]*" L"^-1*" day"^-1*")")))
mtext(expression("Temperature (°C)"),side=1, line=2.5,cex =0.7)
#par(mfrow=c(1,3))
dev.off()

################################################
# remove TN due to the missing value
################################################
#define intercept-only model
intercept_only <- lm(DO_slope ~ 1, data=na.omit(data))

#define model with all predictors
all <- lm(DO_slope ~ DIC + NPOC + TSS+T_mean+TOT_BASIN_AREA+StreamOrde+Transformations, data = na.omit(data))

###################################
#perform forward stepwise regression
forward <- step(intercept_only, direction='forward', scope=formula(all), trace=0)
forward$anova
forward$coefficients


#perform forward stepwise regression
backward <- step(all, direction='backward', scope=formula(all), trace=0)
backward$anova
backward$coefficients

#
bfit<- lm(DO_slope ~T_mean+DIC+TOT_BASIN_AREA, data = cdata)
summary(bfit)

## partial-regression plot
png(file.path('.',paste0('partial_regression_avPlots_3v_new',".png")),
    width = 3, height = 7, units = 'in', res = 600)
par(mfrow=c(3,1),mgp=c(2,1,0),mar=c(3.4,4.1,1,1.5))
#
avPlots(bfit, ~ DIC,id=FALSE,main='',xlab='', #xlab=expression("Residuals - TDN (mg L"^-1*" )"),
        ylab=expression(paste("Residuals - ","ER"[wc]*" (mg O"[2]*" L"^-1*" day"^-1*")")))
mtext(expression("Residuals - DIC (mg L"^-1*" )"),side=1, line=2.5,cex =0.7)
avPlots(bfit, ~ TOT_BASIN_AREA,id=FALSE,main='',xlab='', #xlab=expression("Residuals - Drainage Area (km"^2*")"),
        ylab=expression(paste("Residuals - ","ER"[wc]*" (mg O"[2]*" L"^-1*" day"^-1*")")))
mtext(expression("Residuals - Drainage Area (km"^2*")"),side=1, line=2.5,cex =0.7)
avPlots(bfit, ~ T_mean,id=FALSE,main='', xlab='',#xlab=expression("Residuals - Temperature (Â°C)"),
        ylab=expression(paste("Residuals - ","ER"[wc]*" (mg O"[2]*" L"^-1*" day"^-1*")")))
mtext(expression("Residuals - Temperature (°C)"),side=1, line=2.5,cex =0.7)
dev.off()
# 

################################################
# read in  ERtotal data
ERriv <- read.csv(file.path('./Data','mean_ERvolumetric_best_streamPULSEsites.csv'))
ERriv$ERvolumetric[ERriv$ERvolumetric>0]<-0

#ERwater from SS2021
ERwc<-data

# ERwater data from literture
ERwc2 <- read.csv(file.path('data','ERwc_combined_lit_valuesV2.csv'))

## calculate the skewness for ERtotal and ER water in this study
sk1= round(skewness(ERriv$ERvolumetric),2)
sk2= round(skewness(ERwc$DO_slope),2)

################################################
# make the density plots
## density plot for ERwater in this study
p0 <- ggplot() + 
  geom_density(data=ERwc, aes(x=DO_slope,fill='wc'),color='blue',adjust = 6)+
  geom_vline(aes(xintercept=median(ERwc$DO_slope)),color="blue",  size=1)+
  geom_vline(data=ERwc2, aes(xintercept=ERwc,color='lit'),linetype="dashed")+
  scale_x_cut(breaks=c(-0.13), which=c(1), scales=c(0.25, 1),space = 0.2)+ theme_bw()+ 
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

ggsave(file.path('.',"hist_density_plot_ERwater.png"),
       plot=p0, width = 4, height = 3, dpi = 300,device = "png") #grid.arrange(p1,p2, nrow=1)


## density plot for ERwater with legend
p1 <- ggplot() + 
  geom_density(data=ERwc, aes(x=DO_slope,fill='wc'),color='blue',adjust = 4)+
  geom_vline(aes(xintercept=median(ERwc$DO_slope)),color="blue",  size=1)+
  geom_vline(data=ERwc2, aes(xintercept=ERwc,color='lit'),linetype="dashed")+
  scale_x_cut(breaks=c(-0.13), which=c(1), scales=c(0.25, 1),space = 0.2)+ theme_bw()+ 
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
        #legend.key = element_rect(fill = "white", color = "black", linewidth = 0.5),
        legend.box.just = "left")

# Extract the colour legend - leg1
leg1 <- gtable_filter(ggplot_gtable(ggplot_build(p1)), "guide-box") 

ggsave(file.path('.',"hist_density_plot_ERwater_legend.png"),
       plot=p1, width = 4, height = 3, dpi = 300,device = "png") #grid.arrange(p1,p2, nrow=1)

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
  geom_rect(aes(xmin=-4.63,xmax=-0.02,ymin=0,ymax=0.08,colour="lit",fill='lit'))+ #ER lit
  geom_rect(aes(xmin=-0.11,xmax=0,ymin=0,ymax=0.08,colour="wc",fill='wc'),alpha=0.5)+ #ER WC
  scale_colour_manual("",breaks = c("tot", "wc", "lit"),labels = c(expression("ER"[tot]*""), expression("ER"[wc]*" range (this study)"), expression("ER"[wc]*" range (Lit) ")),
                      values = c("black", "blue", "#F9847B"),aesthetics = c("colour"))+
  scale_fill_manual("",breaks = c("tot", "wc", "lit"),labels = c(expression("ER"[tot]*""), expression("ER"[wc]*" range (this study)"), expression("ER"[wc]*" range (Lit) ")),
                      values = c("grey", "skyblue", "#F9847B"),guide = guide_legend(override.aes = list(alpha = .5)))+
  theme(
    legend.position = c(.4, .9),
    legend.justification = c( "top"),
    legend.margin = margin(6, 6, 6, 6),
    legend.text = element_text(size=12,hjust = 0, margin = margin(l = 0, r = 5, unit = "pt")),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),
    legend.key = element_rect(fill = "white", color = "black", linewidth = 0.2),
    legend.box.just = "right"
  )
ggsave(file.path('.',"hist_density_plot_ERtot_legend.png"),
       plot=p2, width = 4, height = 3, dpi = 300,device = "png") #grid.arrange(p1,p2, nrow=1)


# density plot for ERtotal with legend
p3 <- ggplot(ERriv, aes(x=ERvolumetric,color='tot',fill="tot")) + 
  geom_density()+ 
  geom_vline(aes(xintercept=median(ERvolumetric)), color='black', size=0.8)+ 
  # xlab(expression("ER"[tot]*"")) +
  # ylab('Density') + scale_fill_grey() + 
  theme_classic()+
  labs(x = expression("ER"[tot]*" (mg O"[2]*" L"^-1*" d"^-1*")"), y = 'Density', color = "Legend")+
  # geom_rect(aes(xmin=-4.63,xmax=-0.02,ymin=0,ymax=0.08,colour="lit",fill='lit'))+ #ER lit
  annotate("rect", xmin=-4.63,xmax=-0.02,ymin=0,ymax=0.08, alpha=0.6, fill="#F9847B") +
  geom_rect(aes(xmin=-0.11,xmax=0,ymin=0,ymax=0.08,colour="wc",fill='wc'),alpha=0.5)+ #ER WC
  #annotate("rect", xmin=-0.11,xmax=0,ymin=0,ymax=0.08,  fill="blue") +
  scale_colour_manual("",breaks = c("tot", "wc", "lit"),labels = c(expression("ER"[tot]*""), expression("ER"[wc]*" range (this study)"), expression("ER"[wc]*" range (Lit) ")),
                      values = c("black", "blue", "#F9847B"),aesthetics = c("colour"))+
  scale_fill_manual("",breaks = c("tot", "wc", "lit"),labels = c(expression("ER"[tot]*""), expression("ER"[wc]*" range (this study)"), expression("ER"[wc]*" range (Lit) ")),
                    values = c("grey", "lightblue", "#F9847B"))+
  theme(legend.position ="none")

ggsave(file.path('.',"hist_density_plot_bottom_ERtot_nolegend.png"),
       plot=p3, width = 4, height = 3, dpi = 300,device = "png") #grid.arrange(p1,p2, nrow=1)


