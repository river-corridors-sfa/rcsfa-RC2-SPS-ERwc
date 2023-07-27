# RC2 spatial study - Multiple linear regression for manual_chamber_data
# X Lin Nov 22 2022
################################################
# Read in data 
################################################
rm(list=ls(all=TRUE))
#setwd('C:/Users/linx882/XLin/automation of respiration calculations/RC_2/codes/Spatial_Study')
setwd('C:/Users/linx882/OneDrive - PNNL/XLin/automation of respiration calculations/spatial_study_21_manual_chamber_data_analysis-main')
library(MASS)
library(relaimpo)
library(visreg)
library(ggstatsplot)
library(caret)
library(leaps)
library(car)
library(patchwork)
# read in metadata
data<-read.csv(file.path('data','spatial_data.csv'))
names(data)[c(1:7)]<-c('Site_ID','Sample_Name','DIC','NPOC','TN','TSS','DO_slope')
data$DO_slope[data$DO_slope>0]<-0
sdata<- read.csv(file.path('data','SPS_Total_and_Normalized_Transformations_01-03-23.csv'))
names(sdata)<-c('Sample_Name','Transformations','Peaks','Normalized_Transformations')
data <- merge(data,sdata,by=c("Sample_Name"))

chemdata <- read.csv(file.path('data','v2_SFA_SpatialStudy_2021_SampleData','v2_SPS_NPOC_TN_DIC_TSS_Ions_Summary.csv'),skip=2)
chemdata <-chemdata[grep('SPS',chemdata$Sample_Name),]
names(chemdata)
chemdata <-chemdata[,c(2,4,18:20)]; 
names(chemdata)<-c('Sample_Name','DIC','NPOC','TN','TSS')
chemdata[chemdata=='-9999'] =NA
chemdata[c('DIC','NPOC','TN','TSS')] <- sapply(chemdata[c('DIC','NPOC','TN','TSS')],as.numeric)

data <-merge(data[,-grep(paste(c('DIC','NPOC','TN','TSS'), collapse = "|"),names(data))],chemdata,by=c("Sample_Name"))


cdata <- na.omit(data)

vars <- c('DIC','NPOC', 'TN','TSS','T_mean','TOT_BASIN_AREA','StreamOrde','Normalized_Transformations')
# correlation matrix
png(file.path('results',paste0('exploratory_variables_correlation_matrix',".png")),
    width = 10, height = 10, units = 'in', res = 600)
#par(mfrow=c(2,2)) 
chart.Correlation(data[vars], histogram=TRUE, pch=19)
dev.off()


################################################
# Stepwise Regression
# fit <- lm(DO_slope ~ DIC + NPOC + TN + TSS+T_mean+TOT_BASIN_AREA+StreamOrde, data = na.omit(cdata))
# step <- stepAIC(fit, direction="both")
# step$anova # display results


#define intercept-only model
intercept_only <- lm(DO_slope ~ 1, data=cdata)

#define model with all predictors
all <- lm(DO_slope ~ DIC + NPOC + TN+TSS+T_mean+TOT_BASIN_AREA+StreamOrde+Normalized_Transformations, data = cdata)

###################################
#perform forward stepwise regression
forward <- step(intercept_only, direction='forward', scope=formula(all), trace=1)
forward$anova
forward$coefficients
png(file.path('results',paste0('forward_best_coeff',".png")),
    width = 6, height = 5, units = 'in', res = 600)
#par(mfrow=c(1,2)) 
ggcoefstats(forward)
dev.off()
#############
#  lm fitting using selected variables from forward stepwise selection
bfit<- lm(DO_slope ~TN +TOT_BASIN_AREA+ T_mean, data = data)
#bfit<- lm(DO_slope ~TN +TOT_BASIN_AREA+ T_mean+StreamOrde+Transformations, data = data)
summary(bfit)

png(file.path('results',paste0('bestfit_plots3',".png")),
    width = 6, height = 6, units = 'in', res = 600)
par(mfrow=c(2,2)) 
plot(bfit)
dev.off()


png(file.path('results',paste0('avPlots_all_variables_sc3',".png")),
    width = 3, height = 6, units = 'in', res = 600)
avPlots(bfit,id=FALSE,layout=c(3,1)) 
dev.off()

png(file.path('results',paste0('crPlots_all_variables_sc3',".png")),
    width = 3, height = 6, units = 'in', res = 600)
crPlots(bfit,ylab='Partial-Residual (DO_slope)',main='Partial-Residual Plots', smooth=FALSE,id=FALSE,layout=c(3,1))
dev.off()

png(file.path('results',paste0('avPlots_3v_new',".png")),
    width = 3, height = 7, units = 'in', res = 600)
par(mfrow=c(3,1),mgp=c(2,1,0),mar=c(3.4,4.1,1,1.5))
#
avPlots(bfit, ~ TN,id=FALSE,main='',xlab='', #xlab=expression("Residuals - TDN (mg L"^-1*" )"),
        ylab=expression(paste("Residuals - ","ER"[wc]*" (mg O"[2]*" L"^-1*" day"^-1*")")))
mtext(expression("Residuals - TDN (mg L"^-1*" )"),side=1, line=2.5,cex =0.7)
avPlots(bfit, ~ TOT_BASIN_AREA,id=FALSE,main='',xlab='', #xlab=expression("Residuals - Drainage Area (km"^2*")"),
        ylab=expression(paste("Residuals - ","ER"[wc]*" (mg O"[2]*" L"^-1*" day"^-1*")")))
mtext(expression("Residuals - Drainage Area (km"^2*")"),side=1, line=2.5,cex =0.7)
avPlots(bfit, ~ T_mean,id=FALSE,main='', xlab='',#xlab=expression("Residuals - Temperature (°C)"),
        ylab=expression(paste("Residuals - ","ER"[wc]*" (mg O"[2]*" L"^-1*" day"^-1*")")))
mtext(expression("Residuals - Temperature (°C)"),side=1, line=2.5,cex =0.7)
#par(mfrow=c(1,3))
dev.off()
# 
# png(file.path(paste0('avPlots_StreamOrde',".png")),
#     width = 6, height = 6, units = 'in', res = 600)
# avPlots(bfit, ~ StreamOrde,id=FALSE)
# dev.off()


# 
png(file.path('results',paste0('crPlots_3v_new',".png")),
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


###################################
#perform forward stepwise regression
backward <- step(all, direction='backward', scope=formula(all), trace=0)
backward$anova
backward$coefficients
png(file.path('results',paste0('backward_best_coeff',".png")),
    width = 8, height = 6, units = 'in', res = 600)
#par(mfrow=c(1,2)) 
ggcoefstats(backward)
dev.off()




################################################
# remove TN due to the missing value
#define intercept-only model
intercept_only <- lm(DO_slope ~ 1, data=na.omit(data))

#define model with all predictors
all <- lm(DO_slope ~ DIC + NPOC + TSS+T_mean+TOT_BASIN_AREA+StreamOrde+Transformations, data = na.omit(data))

###################################
#perform forward stepwise regression
forward <- step(intercept_only, direction='forward', scope=formula(all), trace=0)
forward$anova
forward$coefficients
png(file.path('results',paste0('forward_best_coeff',".png")),
    width = 6, height = 5, units = 'in', res = 600)
#par(mfrow=c(1,2)) 
ggcoefstats(forward)
dev.off()

#perform forward stepwise regression
backward <- step(all, direction='backward', scope=formula(all), trace=0)
backward$anova
backward$coefficients
png(file.path('results',paste0('backward_best_coeff',".png")),
    width = 8, height = 6, units = 'in', res = 600)
#par(mfrow=c(1,2)) 
ggcoefstats(backward)
dev.off()

#############
bfit<- lm(DO_slope ~T_mean+DIC, data = cdata)
summary(bfit)
visreg(bfit, gg=TRUE)
png(file.path('results',paste0('bestfit_coeff',".png")),
    width = 6, height = 6, units = 'in', res = 600)
#par(mfrow=c(1,2)) 
ggcoefstats(bfit)
dev.off()


################################################
library(scales)
library(moments)
library(plotrix)
library(ggplot2)
library(gridExtra)
library(gtable)
library(grid)
library(ggbreak)

ERriv <- read.csv(file.path('data','mean_ERvolumetric_best_streamPULSEsites.csv'))
ERriv$ERvolumetric[ERriv$ERvolumetric>0]<-0
ERwc<-data

ERwc2 <- read.csv(file.path('data','ERwc_combined_lit_valuesV2.csv'))


sk1= round(skewness(ERriv$ERvolumetric),2)
sk2= round(skewness(ERwc$DO_slope),2)



colors <- c(expression("median ER"[wc]*"") = "blue", expression("ER"[wc]*" lit") = "black")
p1 <- ggplot(ERwc, aes(x=DO_slope)) + 
  geom_density(color="darkblue", fill="lightblue")+ 
  geom_vline(aes(xintercept=median(DO_slope)),color="blue",  size=1)+
  geom_vline(data=ERwc2, aes(xintercept=ERwc),color='red',linetype="dashed")
p1 + scale_x_cut(breaks=c(-0.13), which=c(1), scales=c(0.25, 1),space = 0.2)+ theme_bw()+
  xlab(expression("ER"[wc]*"")) +
  ylab('Density') + scale_fill_grey() + theme_classic()


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

ggsave(file.path('results',"hist_density_plot_top_gg_c0.png"),
       plot=p0, width = 4, height = 3, dpi = 300,device = "png") #grid.arrange(p1,p2, nrow=1)

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

ggsave(file.path('results',"hist_density_plot_top_gg_c1.png"),
       plot=p1, width = 4, height = 3, dpi = 300,device = "png") #grid.arrange(p1,p2, nrow=1)

# plotNew <- p0 + 
#   annotation_custom(grob = leg1, xmin = -0.075, xmax = -0.045, ymin = 20, ymax = 30)


p2 <- ggplot(ERriv, aes(x=ERvolumetric,color='tot',fill="tot")) + 
  geom_density()+ 
  geom_vline(aes(xintercept=median(ERvolumetric)), color='black', size=1)+ 
  # xlab(expression("ER"[tot]*"")) +
  # ylab('Density') + scale_fill_grey() + 
  theme_classic()+
  labs(x = expression("ER"[tot]*" (mg O"[2]*" L"^-1*" d"^-1*")"), y = 'Density', color = "Legend")+
  geom_rect(aes(xmin=-4.63,xmax=-0.02,ymin=0,ymax=0.08,colour="lit",fill='lit'))+ #ER lit
  geom_rect(aes(xmin=-0.11,xmax=0,ymin=0,ymax=0.08,colour="wc",fill='wc'))+ #ER WC
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
ggsave(file.path('results',"hist_density_plot_bottom_gg_legend2.png"),
       plot=p2, width = 4, height = 3, dpi = 300,device = "png") #grid.arrange(p1,p2, nrow=1)




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

ggsave(file.path('results',"hist_density_plot_bottom_gg2_NOLEGEND_c2.png"),
       plot=p3, width = 4, height = 3, dpi = 300,device = "png") #grid.arrange(p1,p2, nrow=1)



bigplot1 <- arrangeGrob(p1, p2,nrow=2)
ggsave(file.path('results',"hist_density_plot_gg1.png"),
       plot=bigplot1, width = 6, height = 8, dpi = 300,device = "png") #grid.arrange(p1,p2, nrow=1)

bigplot2 <- arrangeGrob(p1, p2,nrow=2)
ggsave(file.path('results',"hist_density_plot_gg2.png"),
       plot=bigplot2, width = 6, height = 4, dpi = 300,device = "png") #grid.arrange(p1,p2, nrow=1)
