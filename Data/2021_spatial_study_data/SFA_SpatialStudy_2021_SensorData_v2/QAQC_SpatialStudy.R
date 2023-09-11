# RC2 spatial study
# time series plots and linear regression for BarotrollAtm, MantaRiver and Minidot data
# summary statistics for sensor data
# violin plots for sensor data 
# X Lin Mar 21 2022
########################################################################## 
# Read in metadata 
########################################################################## 
rm(list=ls(all=TRUE))
source('./helper.R')
# read in metadata
mindir<-getwd()
mfile<- 'SPS_Sensor_Field_Metadata.csv' # name of metadata file
metadata<-read_metadata(mindir,mfile)

##########################################################################  
## make time series plot for BarotrollAtm data
## fit linear regression 
## save plots and summary statistics
indir = file.path(getwd(),'BarotrollAtm/Data')
plot_dir= file.path(getwd(),'BarotrollAtm/tsplots')
if(!dir.exists(plot_dir)==T){dir.create(plot_dir)}
summary_dir = file.path(getwd(),'BarotrollAtm')
#
fsdata = BA_tsplot_summary(indir,plot_dir,summary_dir)

##make violin plot for BarotrollAtm data 
plot_dir= file.path(getwd(),'BarotrollAtmData/vplots')
if(!dir.exists(plot_dir)==T){dir.create(plot_dir)}
BA_vplot(indir,plot_dir,summary_dir)

##########################################################################  
## make time series plot for MantaRiver data
## fit linear regression 
## save plots and summary statistics
indir = file.path(getwd(),'MantaRiver/Data')
plot_dir= file.path(getwd(),'MantaRiver/tsplots')
if(!dir.exists(plot_dir)==T){dir.create(plot_dir)}
summary_dir = file.path(getwd(),'MantaRiver')
#
fsdata = MR_tsplot_summary(indir,plot_dir,summary_dir)

##make violin plot for MantaRiver data 
plot_dir= file.path(getwd(),'MantaRiver/vplots')
if(!dir.exists(plot_dir)==T){dir.create(plot_dir)}
MR_vplot(indir,plot_dir,summary_dir)

##########################################################################  
## make time series plot for Minidot data
## fit linear regression 
## save plots and summary statistics
indir = file.path(getwd(),'MinidotManualChamber/Data')
plot_dir= file.path(getwd(),'MinidotManualChamber/tsplots')
if(!dir.exists(plot_dir)==T){dir.create(plot_dir)}
summary_dir =  file.path(getwd(),'MinidotManualChamber')
#
fsdata = MD_tsplot_summary(indir,plot_dir,summary_dir)

##make violin plot for Minidot data
plot_dir=file.path(getwd(),'MinidotManualChamber/vplots')
if(!dir.exists(plot_dir)==T){dir.create(plot_dir)}
MD_vplot(indir,plot_dir,summary_dir)


