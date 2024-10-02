library(tidyverse)
library(moments)
library(ggbreak)
library(gtable)


current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd("../..")
getwd()

################################################
# read in ERwc data 

data = read.csv(file.path('./Data/Multiple_linear_regression/spatial_data.csv')) %>% 
  select(-c(Water_Column_Respiration, Temperature_Mean))

erwc_sps = read.csv(file.path('./Data/Multiple_linear_regression/ERwc_Mean.csv')) %>% 
  select(-c(X)) %>% 
  mutate(Mean_ERwc = round(Mean_ERwc, 3)) %>% 
  mutate(Mean_Temp = round(Mean_Temp, 3)) %>% 
  left_join(data, by = "Site_ID")

# read in  ERtotal data
ERriv = read.csv(file.path('./Data/Appling_ERtot_analysis','mean_ERtot_bestSiteIDs.csv'))
#ERriv$ERvolumetric[ERriv$ERvolumetric>0]<-0

# ERwater data from literture
#ERwc2 <- read.csv(file.path('./Data/Multiple_linear_regression','ERwc_combined_lit_valuesV3.csv'))
erwc_lit <- read.csv(file.path('./Data/Water_column_respiration_published','Water_column_respiration_published_values.csv'))

## calculate the skewness for ERtotal and ER water in this study
sk1= round(skewness(ERriv$Total_Ecosystem_Respiration_Volumetric),2)
sk2= round(skewness(erwc_sps$Mean_ERwc),2)

################################################
# make the density plots
## density plot for ERwater in this study
p0 <- ggplot() + 
  geom_density(data=erwc_sps, aes(x=Mean_ERwc,fill='wc'),color='blue',adjust = 6)+
  geom_vline(aes(xintercept=median(erwc_sps$Mean_ERwc)),color="blue",  size=1)+
  geom_vline(data=erwc_lit, aes(xintercept= Water_Column_Respiration_Literature,color='lit'),linetype="dashed")+
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

ggsave(file.path('./Plots',"hist_density_plot_ERwater_ML.png"),
       plot=p0, width = 4, height = 3, dpi = 300,device = "png") #grid.arrange(p1,p2, nrow=1)


## density plot for ERwater with legend
# why is adjust different? what does this do?
p1 <- ggplot() + 
  geom_density(data=erwc_sps, aes(x=Mean_ERwc,fill='wc'),color='blue',adjust = 4)+
  geom_vline(aes(xintercept=median(erwc_sps$Mean_ERwc)),color="blue",  size=1)+
  geom_vline(data=erwc_lit, aes(xintercept=Water_Column_Respiration_Literature,color='lit'),linetype="dashed")+
  #scale_x_cut(breaks=c(-0.12), which=c(1), scales=c(0.25, 1),space = 0.2)+ theme_bw()+ 
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

ggsave(file.path('./Plots',"hist_density_plot_ERwater_legend_ML.png"),
       plot=p1, width = 4, height = 3, dpi = 300,device = "png") #grid.arrange(p1,p2, nrow=1)

# plotNew <- p0 + 
#   annotation_custom(grob = leg1, xmin = -0.075, xmax = -0.045, ymin = 20, ymax = 30)

# density plot for ERtotal 
# figure out ymins/maxes for ERwc_sps and ERwc_lit
p2 <- ggplot(ERriv, aes(x=Total_Ecosystem_Respiration_Volumetric,color='tot',fill="tot")) + 
  geom_density()+ 
  geom_vline(aes(xintercept=median(Total_Ecosystem_Respiration_Volumetric)), color='black', size=1)+ 
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


ggsave(file.path('./Plots',"hist_density_plot_ERtot_legend.png"),
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
  annotate("rect", xmin=-4.63,xmax=-0.02,ymin=0,ymax=0.03, alpha=0.6, fill="#F9847B") +
  geom_rect(aes(xmin=-0.11,xmax=0,ymin=0,ymax=0.03,colour="wc",fill='wc'),alpha=0.5)+ #ER WC
  #annotate("rect", xmin=-0.11,xmax=0,ymin=0,ymax=0.08,  fill="blue") +
  scale_colour_manual("",breaks = c("tot", "wc", "lit"),labels = c(expression("ER"[tot]*""), expression("ER"[wc]*" range (this study)"), expression("ER"[wc]*" range (Lit) ")),
                      values = c("black", "blue", "#F9847B"),aesthetics = c("colour"))+
  scale_fill_manual("",breaks = c("tot", "wc", "lit"),labels = c(expression("ER"[tot]*""), expression("ER"[wc]*" range (this study)"), expression("ER"[wc]*" range (Lit) ")),
                    values = c("grey", "lightblue", "#F9847B"))+
  xlim(-20, max(ERriv$ERtot))+
  theme(legend.position ="none")

ggsave(file.path('./Plots',"hist_density_plot_bottom_ERtot_nolegend.png"),
       plot=p3, width = 4, height = 3, dpi = 300,device = "png") #grid.arrange(p1,p2, nrow=1)

##############################################################
# use lme4 to estimate the means among the sites
rdata <- read.csv(file.path('./Data/Multiple_linear_regression/Minidot_Summary_Statistics_v2.csv'),skip=48)
rdata <- rdata[c('Site_ID','Dissolved_Oxygen_1_Slope','Dissolved_Oxygen_2_Slope','Dissolved_Oxygen_3_Slope')]
rdata2 <- rdata %>% 
  pivot_longer(
    cols = names(rdata)[2:4], 
    names_to = "replicates",
    values_to = "ERwc"
  )
rdata2$ERwc<- rdata2$ERwc*60*24

## within-sample variation
wsummary<-rdata2%>%group_by(Site_ID) %>%
  summarise(mean = mean(ERwc),
            sd = sd(ERwc),
            var = var(ERwc))
mean(wsummary$sd)
mean(wsummary$var)
## lme4 fitting

lfit <- lmer( ERwc ~  (1 | Site_ID), data=rdata2)  
summary(lfit)


###############################################################
# plot to visualize the value of ERwater in each Site
#odata=data[c('Site_ID','ERwc')]
odata=rdata2
odata <- odata %>%
  group_by(Site_ID) %>%
  mutate(mean = mean(ERwc),median =median(ERwc) )

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
  geom_dotplot(binaxis='y', stackdir='center',binwidth = 0.002,dotsize = 0.8) +
  ylab(expression("ER"[wc]*" (mg O"[2]*" L"^-1*"day"^-1*")")) + xlab("Site ID") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = c(0, 0))+
  scale_x_discrete(labels=sites)
DotPlotFin <- DotPlot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(DotPlotFin)
ggsave(file.path("./Plots",paste0('ERwc_dotplot_mean_rank',".png")), 
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
  geom_dotplot(binaxis='y', stackdir='center',binwidth = 0.002,dotsize = 0.8) +
  ylab(expression("ER"[wc]*" (mg O"[2]*" L"^-1*"day"^-1*")")) + xlab("Site ID") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = c(0, 0))+
  scale_x_discrete(labels=sites)
DotPlotFin <- DotPlot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(DotPlotFin)
ggsave(file.path("./Plots",paste0('ERwc_dotplot_median_rank',".png")), 
       plot=DotPlotFin, width = 6, height = 3, dpi = 300,device = "png") #grid.arrange(p1,p2, nrow=1)