# RC2 spatial study
# Functions to help visualize BarotrollAtm, MantaRiver and Minidot data
# X Lin Mar 21th 2022
########################################################################## 
#rm(list=ls(all=TRUE))

libs <- c('ncdf4','raster','maptools','rgdal','GISTools','readxl',
          'tidyverse','tidyr', 'ggplot2','sp','lubridate','feather','vioplot')
for(il in 1:length(libs)){
  if (!is.element(libs[il], installed.packages()[,1])){
    install.packages(libs[il], dep = TRUE)}
}

Sys.setenv(TZ='GMT')
library(raster)
library(ncdf4)
library(maptools)
library(rgdal)
library(GISTools)
library(tidyverse)
library(tidyr)
library(sp)
library(ggplot2)
library(vioplot)
library(lubridate)
library(stringr)
library(readxl)
##########################################################################
## function to read metadata 
read_metadata<-function(indir,file){
  ## indir: : directory where metadata are saved
  ## file: filename of metadata
  metadata <- read.csv(file.path(indir,file))
  names(metadata)[grep("Date", names(metadata))[1]]<-'Date'
  names(metadata)[grep("Time", names(metadata))]<-gsub('_.24_hr_hh_mm.','',names(metadata)[grep("Time", names(metadata))])
  names(metadata)[grep("Time", names(metadata))]<-gsub('_.24_hr_hh.mm.','',names(metadata)[grep("Time", names(metadata))])
  if (any(sub('.*\\/', '', metadata$Date)==2021)==TRUE) {
    metadata['Date'] <-as.character(as.Date(metadata$Date,format="%m/%d/%Y"))
  }else{
    metadata['Date'] <-as.character(as.Date(metadata$Date,format="%m/%d/%y"))
  }
  #metadata['Date'] <-as.character(as.Date(metadata$Date,format="%m/%d/%Y"))
  scols= grep(paste(c('Manta_Time','Manta_Start','Manta_End'), collapse = "|"),names(metadata))
  for (s in 1:length(scols)){
    sdate <- metadata[,scols[s]]
    tid1 <- which(str_count(sdate, pattern = "/")<2); tid2 <- which(str_count(sdate, pattern = "/")==2)
    sdate[tid1]<- paste(metadata$Date[tid1],sdate[tid1])
    metadata[tid1,scols[s]]<-as.character(as.POSIXct(sdate[tid1],"%Y-%m-%d %H:%M",tz='GMT'))
    metadata[tid2,scols[s]]<-as.character(as.POSIXct(sdate[tid2],"%m/%d/%Y %H:%M",tz='GMT'))
  }
  #grepl(paste(c('00065_cd',"00060_cd"), collapse = "|"),names(hydro_data))
  #metadata['Date'] <-as.character(as.Date(metadata$Date,format="%m/%d/%y"))
  # metadata$Location[which(metadata$Location=="Cle Elum")]="Cle_Elum"
  # metadata$Location[which(metadata$Location=="Temporal Sites")]="Temporal_Sites"
  return(metadata)
}

##########################################################################
## calculated r2
r2_general <-function(preds,actual){ 
  return(1- sum((preds - actual) ^ 2)/sum((actual - mean(actual))^2))
}

## find outliers
find_anomaly<-function(data,col){
  idx = which(colnames(data)==col)
  tdata<-data.frame(time=data$DateTime,var = data[,idx])
  d_ts<- as_tbl_time(tdata, time)
  tso<-d_ts%>%
    #time_decompose(var, method = "stl",trend = "10 minutes") %>%
    anomalize(var, method = "iqr", alpha = 0.03, max_anoms = 0.5) #iqr
  outliers=which(tso$anomaly=='Yes')
  return(outliers)
}

##########################################################################
## function to fit regression line to test data, and return summary data 
BAfitplot<-function(data,col,ylab,outliers=0){
  ## sdata:  test data of specific ID with three replicates (e.g. S13R)
  ## col: column name
  ## ylab: yaxis label for the selected column 'col'
  ## outliers: index of outliers
  colors = c('black','blue','red')
  ltys = c(1,2,3)
  #slope <- c();rmse <-c();r2<-c();drange<-c()
  cid = which(colnames(data)==col) #match(col, names(data))
  #data = sdata[,cid]
  data['x'] = as.numeric(data$DateTime) #c(1:nrow(data))
  fit_l=lm(data[,cid]~x,data) #grep("^T$", colnames(data))
  plot(data$DateTime,data[,cid],
       xaxt="n",
       xlab='DateTime',
       ylab=ylab, 
       #main = paste0("Site ID ",unique(sdata$ID)),
       ylim=c(min(data[col]),max(data[col])),
       xlim=c(min(data$DateTime),max(data$DateTime)),
       type='b',#pch= c(16), 
       col='black')
  if (length(outliers)>0){
    points(data$DateTime[outliers],data[outliers,cid],col='red',pch=8,cex=2)
  }
  axis.POSIXct(1,at=seq(from=min(data$DateTime), to=max(data$DateTime), by = 1800),
               labels=format(seq(from=min(data$DateTime), to=max(data$DateTime), by = 1800),"%m/%d %H:%M"),las=1)
  lines(data$x,fit_l$fitted.values,col='red', lty=2, lwd=2)
  slope<-format(summary(fit_l)$coefficients[2], digits=2)
  r2<-format(summary(fit_l)$r.squared, digits=2)
  rmse<- format(sqrt(mean((data[,cid] - fit_l$fitted.values)^2)), digits=3)
  nrmse<- round(sqrt(mean((data[,cid] - fit_l$fitted.values)^2))/mean(data[,cid],na.rm=TRUE), digits=3)
  dmean<- format(mean(data[,cid]), digits=3)
  drange <-format(max(data[,cid])-min(data[,cid]), digits=3)
  legend(x = "top",
         inset=c(1.1,-0.16), xpd=TRUE,
         lty=2, bty='n',horiz=TRUE,
         #legend=paste0('ID ',unique(data$ID),"\nSlope=",slope,"\nRMSE=",rmse,'\n'),
         legend=paste0('Site-ID ',unique(data$ID),"  Slope=",slope,"  RMSE=",rmse),
         col='red', #pch=c(16, 17),
         text.col = 'red',
         cex=0.8,pt.cex=1) #bty="n",,ncol=3
  return(list(r=slope,r2=r2,rmse=rmse,nrmse=nrmse,mean=dmean,range= drange))
}

##########################################################################
## time series plot and linear regression for BarotrollAtm data
BA_tsplot_summary<-function(indir,plot_dir,summary_dir){
  ## indir: : directory where test data are saved
  ## ids: IDs of available test data  (e.g. S13R)
  ## ids_patterns: ID patterns used to select test data files
  ## plot_dir: directory to save the generated plots
  ## summary_dir: directory to save lm regression summary data 
  test_files = list.files(file.path(indir),pattern = c(".csv"))
  #test_files =grep(paste(test_dates, collapse = "|"),test_files,value = TRUE)
  ids =gsub("_[^_]*$","",test_files,perl=T) #unique(sub("\\_.*", "", test_files))
  sidx<-sort(parse_number(ids), index.return=TRUE)$ix
  ids =ids[sidx] #ids_patterns = paste0(ids,'_')
  mdata =list()
  for (i in 1:length(ids)){
    sfiles <-list.files(file.path(indir),pattern = ids[i])
    date <- str_extract(sfiles, "[0-9]{4}-[0-9]{2}-[0-9]{2}")
    sdata<-data.frame() #sdata <-list() #
    for (j in 1:length(sfiles)){
      lines <- readLines(file.path(indir,sfiles[j]))
      sidx<-grep("DateTime", lines)[2]
      data <- read.delim(file.path(indir,sfiles[j]),skip = sidx-1,sep = ",",header = TRUE)
      #data <- read.csv(file.path(indir,sfiles[j]),header = TRUE)
      names(data) <- sub("\\..*", "",names(data))
      data$DateTime  <- as.POSIXct(trimws(data$DateTime, "l"),"%Y-%m-%d %H:%M:%S",tz='GMT')
      # int <- interval(min(data$Time)+minutes(8),max(data$Time)-minutes(5))
      # data <- data[data$Time%within%int,]
      sdata<-rbind(sdata,data)  #sdata[[j]] <-data ##
    }
    if (nrow(sdata)>0){
      #if (length(sdata)>0){
      cols<- c('Pressure','Air_Temperature')
      ylabs = c("Pressure(mBar)","Temperature(deg C)")
      png(file.path(plot_dir, #
                    paste0('TsPlot_',ids[i],'_',date,'.png')), #'Testdata_',
          width = 6, height = 4, units = 'in', res = 300)
      par(mfrow=c(2,1),mgp=c(2,1,0),mar=c(2.1,3.1,1,1.5))
      ### plot with linear fit
      rdata=data.frame(matrix(0, nrow = 6,ncol=length(cols), byrow = TRUE))
      names(rdata) = cols
      for (k in 1:length(cols)){
        #idx = grep(cols[k],names(sdata))
        idx = which(colnames(data)==cols[k])
        if(length(idx)==1){
          #names(sdata)[cid] <-fnames[k]
          outliers = find_anomaly(data=sdata,col=cols[k])
          if (length(outliers)>0){
            cat(paste0('Variable: ',cols[k],'\n'))
            print(sdata[outliers,c(1,idx)])
            cat(paste0('\n'))
            rd=BAfitplot(sdata,col =cols[k],ylab=ylabs[k],outliers=outliers)
          }else{
            cat(paste0('Variable: ',cols[k],'\n'))
            print('None was found!')
            cat(paste0('\n'))
            rd=BAfitplot(sdata,col =cols[k],ylab=ylabs[k])
          }
          #
          rdata[,k]<-matrix(Reduce(c,rd))
        }
      }
      dev.off()
      ## save lm fit summary data
      #rds = mapply(c, rd1, rd2,rd3, SIMPLIFY=FALSE)
      #d <- data.frame(matrix(Reduce(c,rd), nrow = 4, byrow = TRUE))
      #cnames<- apply(expand.grid(c(1,2,3), cols), 1, function(x) paste(x[2], x[1], sep="_")) 
      #names(d)<- cnames
      rdata['SID']<- ids[i];rdata['Test_date']<-date;
      rdata['Metrics']<- c('lm_slope','R_square','RMSE','NRMSE','Mean','Range')
      rdata<- rdata[c('SID','Test_date','Metrics',cols)]
      mdata[[i]] <-rdata
    }
  }
  df <- Reduce(rbind, mdata)
  names(df)<- c('Site_ID','DateTime',	'Metrics',	'Pressure','Air_Temperature')
  #outdir = 'C:/Users/linx882/XLin/automation of respiration calculations'
  write.csv(df,file.path(summary_dir,paste0('RC2_BarotrollAtm_lm_summary','.csv')),
            na = "NA",row.names = FALSE)
  return(df)
}


##########################################################################
## violin plot to check the distribution of variables for BarotrollAtm data
BA_vplot<-function(indir,plot_dir,summary_dir){
  ## indir: : directory where test data are saved
  ## plot_dir: directory to save the generated plots
  ## summary_dir: directory to save lm regression summary data 
  test_files = list.files(file.path(indir),pattern = c(".csv"))
  #test_files =grep(paste(test_dates, collapse = "|"),test_files,value = TRUE)
  ids =gsub("_[^_]*$","",test_files,perl=T) #unique(sub("\\_.*", "", test_files))
  sidx<-sort(parse_number(ids), index.return=TRUE)$ix
  ids =ids[sidx] #ids_patterns = paste0(ids,'_')
  mdata =list()
  for (i in 1:length(ids)){
    sfiles <-list.files(file.path(indir),pattern = ids[i])
    date <- str_extract(sfiles, "[0-9]{4}-[0-9]{2}-[0-9]{2}")
    sdata<-data.frame() #sdata <-list() #
    for (j in 1:length(sfiles)){
      #data <- read.csv(file.path(indir,sfiles[j]),header = TRUE)
      lines <- readLines(file.path(indir,sfiles[j]))
      sidx<-max(grep("DateTime", lines))#[which(grep("Eureka_Manta", lines)<50)])
      #data <- read.delim(file.path(indir,date,location,file),skip = sidx-1,sep = ",",header = TRUE)
      data <- read.csv(file.path(indir,sfiles[j]),skip = sidx-1,header = TRUE)
      names(data) <- sub("\\..*", "",names(data))
      data$DateTime  <- as.POSIXct(trimws(data$DateTime, "l"),"%Y-%m-%d %H:%M:%S",tz='GMT')
      # int <- interval(min(data$Time)+minutes(8),max(data$Time)-minutes(5))
      # data <- data[data$Time%within%int,]
      sdata<-rbind(sdata,data)  #sdata[[j]] <-data ##
    }
    if (nrow(sdata)>0){
      #if (length(sdata)>0){
      cols<- c('Pressure','Air_Temperature')
      ylabs = c("Pressure(mBar)","Temperature(deg C)")
      png(file.path(plot_dir, #
                    paste0('Vplot_',ids[i],'_',date,'.png')),
          width = 4, height = 4, units = 'in', res = 300)
      par(mfrow=c(2,1),mgp=c(2,1,0),mar=c(2.1,3.1,1,1.5))
      for (k in 1:length(cols)){
        idx = which(colnames(data)==cols[k])
        if(length(idx)==1){
          vioplot(sdata[,idx],names=ids[i],
                  ylab = ylabs[k],
                  lineCol = "red",     # Color of the line
                  colMed = "red",col="lightgreen")
          #stripchart(sdata[,idx],method = "jitter",pch = 19,vertical = TRUE,col = "gold",add=TRUE)
        }
      }
      dev.off()
    }
  }
}



##########################################################################
## function to fit regression line to test data, and return summary data 
MRfitplot<-function(data,col,ylab,outliers=0){
  ## sdata:  test data of specific ID with three replicates (e.g. S13R)
  ## col: column name
  ## ylab: yaxis label for the selected column 'col'
  #sid =unique(sdata$Series_ID)
  colors = c('black','blue','red')
  ltys = c(1,2,3)
  #slope <- c();rmse <-c();r2<-c();drange<-c()
  cid = which(colnames(data)==col) #match(col, names(data))
  #data = sdata[,cid]
  data = data[order(data$DateTime),]
  data['x'] =as.numeric(data$DateTime)# c(1:nrow(data))
  fit_l=lm(data[,cid]~x,data) #grep("^T$", colnames(data))
  plot(data$DateTime,data[,cid],
       xaxt="n",
       xlab='DateTime',
       ylab=ylab, 
       ylim=c(min(data[col]),max(data[col])),
       xlim=c(min(data$DateTime),max(data$DateTime)),
       type='b',#pch= c(16), 
       col='black')
  if (length(outliers)>0){
    points(data$DateTime[outliers],data[outliers,cid],col='red',pch=8,cex=2)
  }
  axis.POSIXct(1,at=seq(from=min(data$DateTime), to=max(data$DateTime), by = 1800),
               labels=format(seq(from=min(data$DateTime), to=max(data$DateTime), by = 1800),"%m/%d %H:%M"),las=1)
  lines(data$x,fit_l$fitted.values,col='red', lty=2, lwd=2)
  slope<-format(summary(fit_l)$coefficients[2], digits=2)
  r2<-format(summary(fit_l)$r.squared, digits=2)
  rmse<- format(sqrt(mean((data[,cid] - fit_l$fitted.values)^2)), digits=3)
  nrmse<- round(sqrt(mean((data[,cid] - fit_l$fitted.values)^2))/mean(data[,cid],na.rm=TRUE), digits=3)
  dmean<- format(mean(data[,cid]), digits=3)
  drange <-format(max(data[,cid])-min(data[,cid]), digits=3)
  legend(x = "top",
         inset=c(1.1,-0.12), xpd=TRUE,
         lty=2, bty='n',horiz=TRUE,
         #legend=paste0('ID ',unique(data$ID),"\nSlope=",slope,"\nRMSE=",rmse,'\n'),
         legend=paste0('Site ID ',unique(data$ID),'_',unique(data$Location),"  Slope=",slope,"  RMSE=",rmse),
         col='red', #pch=c(16, 17),
         text.col = 'red',
         cex=0.8,pt.cex=1) #bty="n",,ncol=3
  return(list(r=slope,r2=r2,rmse=rmse,nrmse=nrmse,mean=dmean,range= drange))
}
##########################################################################
## time series plot and linear regression for MantaRiver data
MR_tsplot_summary<-function(indir,plot_dir,summary_dir){
  ## indir: : directory where test data are saved
  ## ids: IDs of available test data  (e.g. S13R)
  ## ids_patterns: ID patterns used to select test data files
  ## plot_dir: directory to save the generated plots
  ## summary_dir: directory to save lm regression summary data 
  ############################################
  test_files = list.files(file.path(indir),pattern = c(".csv"))
  #test_files =grep(paste(test_dates, collapse = "|"),test_files,value = TRUE)
  ids =gsub("_[^_]*$","",test_files,perl=T) #unique(sub("\\_.*", "", test_files))
  sidx<-sort(parse_number(ids), index.return=TRUE)$ix
  ids =ids[sidx] #ids_patterns = paste0(ids,'_')
  mdata =list()
  for (i in 1:length(ids)){
    sfiles <-list.files(file.path(indir),pattern = ids[i])
    date <- str_extract(sfiles, "[0-9]{4}-[0-9]{2}-[0-9]{2}")
    sdata<-data.frame() #sdata <-list() #
    for (j in 1:length(sfiles)){
      #data <- read.csv(file.path(indir,sfiles[j]),header = TRUE)
      lines <- readLines(file.path(indir,sfiles[j]))
      sidx<-max(grep("DateTime", lines))#[which(grep("Eureka_Manta", lines)<50)])
      #data <- read.delim(file.path(indir,date,location,file),skip = sidx-1,sep = ",",header = TRUE)
      data <- read.csv(file.path(indir,sfiles[j]),skip = sidx-1,header = TRUE)
      names(data) <- sub("\\..*", "",names(data))
      data$DateTime  <- as.POSIXct(trimws(data$DateTime, "l"),"%Y-%m-%d %H:%M:%S",tz='GMT')
      # int <- interval(min(data$Time)+minutes(8),max(data$Time)-minutes(5))
      # data <- data[data$Time%within%int,]
      sdata<-rbind(sdata,data)  #sdata[[j]] <-data ##
    }
    if (nrow(sdata)>0){
      #if (length(sdata)>0){
      cols<- c("Temperature",'Depth',"Specific_Conductance","Turbidity",'pH',"Battery")
      ylabs = c("Temperature(deg C)","Depth(m)","Specific_Conductance(uS/cm)","Turbidity(FNU)","pH",'Battery(V)')
      png(file.path(plot_dir, #
                    paste0('TsPlot_',ids[i],'_',date,'.png')), #'Testdata_',
          width = 8, height = 6, units = 'in', res = 300)
      par(mfrow=c(3,2),mgp=c(2,1,0),mar=c(2.1,3.1,2,2.5))
      ### plot with linear fit
      rdata=data.frame(matrix(0, nrow = 6,ncol=length(cols), byrow = TRUE))
      names(rdata) = cols
      for (k in 1:length(cols)){
        #idx = grep(cols[k],names(sdata))
        idx = which(colnames(data)==cols[k])
        if(length(idx)==1){
          #names(sdata)[cid] <-fnames[k]
          outliers = find_anomaly(data=sdata,col=cols[k])
          if (length(outliers)>0){
            cat(paste0('Variable: ',cols[k],'\n'))
            print(sdata[outliers,c(1,idx)])
            cat(paste0('\n'))
            rd=MRfitplot(sdata,col =cols[k],ylab=ylabs[k],outliers=outliers)
          }else{
            cat(paste0('Variable: ',cols[k],'\n'))
            print('None was found!')
            cat(paste0('\n'))
            rd=MRfitplot(sdata,col =cols[k],ylab=ylabs[k])
          }
          #
          rdata[,k]<-matrix(Reduce(c,rd))
        }
      }
      dev.off()
      rdata['SID']<- unique(sdata$Site_ID);rdata['Location']<- unique(sdata$Location);
      rdata['Test_date']<-date;
      rdata['Metrics']<- c('lm_slope','R_square','RMSE','NRMSE','Mean','Range')
      rdata<- rdata[c('SID','Test_date','Metrics',cols)]
      mdata[[i]] <-rdata
    }
  }
  df <- Reduce(rbind, mdata)
  #outdir = 'C:/Users/linx882/XLin/automation of respiration calculations'
  write.csv(df,file.path(summary_dir,paste0('RC2_MantaRiver_lm_summary','.csv')),
            na = "NA",row.names = FALSE)
  return(df)
}

##########################################################################
## violin plot to check the distribution of variables for MantaRiver data
MR_vplot<-function(indir, plot_dir,summary_dir){
  ## indir: : directory where test data are saved
  ## ids: IDs of available test data  (e.g. S13R)
  ## ids_patterns: ID patterns used to select test data files
  ## plot_dir: directory to save the generated plots
  ## summary_dir: directory to save lm regression summary data 
  test_files = list.files(file.path(indir),pattern = c(".csv"))
  #test_files =grep(paste(test_dates, collapse = "|"),test_files,value = TRUE)
  ids =gsub("_[^_]*$","",test_files,perl=T) #unique(sub("\\_.*", "", test_files))
  sidx<-sort(parse_number(ids), index.return=TRUE)$ix
  ids =ids[sidx] #ids_patterns = paste0(ids,'_')
  mdata =list()
  for (i in 1:length(ids)){
    sfiles <-list.files(file.path(indir),pattern = ids[i])
    date <- str_extract(sfiles, "[0-9]{4}-[0-9]{2}-[0-9]{2}")
    sdata<-data.frame() #sdata <-list() #
    for (j in 1:length(sfiles)){
      lines <- readLines(file.path(indir,sfiles[j]))
      sidx<-max(grep("DateTime", lines))#[which(grep("Eureka_Manta", lines)<50)])
      #data <- read.delim(file.path(indir,date,location,file),skip = sidx-1,sep = ",",header = TRUE)
      data <- read.csv(file.path(indir,sfiles[j]),skip = sidx-1,header = TRUE)
      names(data) <- sub("\\..*", "",names(data))
      data$DateTime  <- as.POSIXct(trimws(data$DateTime, "l"),"%Y-%m-%d %H:%M:%S",tz='GMT')
      # int <- interval(min(data$Time)+minutes(8),max(data$Time)-minutes(5))
      # data <- data[data$Time%within%int,]
      sdata<-rbind(sdata,data)  #sdata[[j]] <-data ##
    }
    if (nrow(sdata)>0){
      #if (length(sdata)>0){
      cols<- c("Temperature",'Depth',"Specific_Conductance","Turbidity",'pH',"Battery")
      ylabs = c("Temperature(deg C)","Depth(m)","Specific_Conductance(uS/cm)","Turbidity(FNU)","pH",'Battery(V)')
      png(file.path(plot_dir, #
                    paste0('VPlot_','MantaRiver_',ids[i],'_',date,'.png')), #'Testdata_',
          width = 6, height = 4, units = 'in', res = 300)
      par(mfrow=c(3,2),mgp=c(2,1,0),mar=c(2.1,3.1,1,2.5))
      for (k in 1:length(cols)){
        idx = which(colnames(data)==cols[k])
        if(length(idx)==1){
          vioplot(sdata[,idx],names=ids[i],
                  ylab = ylabs[k],
                  lineCol = "red",     # Color of the line
                  colMed = "red",col="lightgreen")
          #stripchart(sdata[,idx],method = "jitter",pch = 19,vertical = TRUE,col = "gold",add=TRUE)
        }
      }
      dev.off()
    }
  }
}

##########################################################################
## function to fit regression line to test data, and return summary data for Minidot data
MDfitplot<-function(sdata,col,ylab){
  ## sdata:  test data of specific ID with three replicates (e.g. S13R)
  ## col: column name
  ## ylab: yaxis label for the selected column 'col'
  sid =unique(sdata$Minidot_SN)
  colors = c('black','blue','red')
  ltys = c(1,2,3)
  slope <- c();rmse <-c();nrmse <-c();r2<-c();dmean<-c();drange<-c()
  for (i in 1:length(sid)){
    data = sdata[which(sdata$Minidot_SN==sid[i]),]
    cid = which(colnames(data)==col) #match(col, names(data))
    data['x'] = as.numeric(data$DateTime)#c(1:nrow(data))
    fit_l=lm(data[,cid]~x,data) #grep("^T$", colnames(data))
    if (i==1){
      plot(data$DateTime,data[,cid],
           xaxt="n",
           xlab='DateTime',
           ylab=ylab, 
           main = paste0("Site ID ",unique(sdata$ID)),
           ylim=c(min(sdata[col]),max(sdata[col])),
           xlim=c(min(sdata$DateTime),max(sdata$DateTime)),
           type='b',#pch= c(16), 
           col=colors[i])
      axis.POSIXct(1,at=seq(from=min(sdata$DateTime), to=max(sdata$DateTime), by = 1800),
                   labels=format(seq(from=min(sdata$DateTime), to=max(sdata$DateTime), by = 1800),"%m/%d %H:%M"),las=1)
      lines(data$x,fit_l$fitted.values,col=colors[i], lty=ltys[i], lwd=2)
    }else{
      lines(data$DateTime,data[,cid],type='b',col=colors[i])
      lines(data$x,fit_l$fitted.values,col=colors[i], lty=ltys[i], lwd=2)
    }
    slope[i]<-format(summary(fit_l)$coefficients[2], digits=2)
    r2[i]<-format(summary(fit_l)$r.squared, digits=2)
    rmse[i]<- format(sqrt(mean((data[,cid] - fit_l$fitted.values)^2)), digits=3)
    nrmse[i]<- round(sqrt(mean((data[,cid] - fit_l$fitted.values)^2))/mean(data[,cid],na.rm=TRUE), digits=3)
    dmean[i] <-format(mean(data[,cid]), digits=3)
    drange[i] <-format(max(data[,cid])-min(data[,cid]), digits=3)
  }
  legend(x = "left",
         inset=c(1,0), xpd=TRUE,
         lty=ltys, bty='n',#horiz=TRUE,
         ##horiz=TRUE, bty='n',
         legend=paste0('Series ID ',sid,"\nslope=",slope,"\nRMSE=",rmse,'\n'), 
         col=colors, #pch=c(16, 17),
         text.col = colors,
         cex=0.8,pt.cex=1) #bty="n",,ncol=3
  return(list(r=slope,r2=r2,rmse=rmse,nrmse=nrmse,mean=dmean,range= drange))
}

##########################################################################
## time series plot and linear regression for Minidot data
MD_tsplot_summary<-function(indir,plot_dir,summary_dir){
  ## indir: : directory where test data are saved
  ## plot_dir: directory to save the generated plots
  ## summary_dir: directory to save lm regression summary data
  test_files = list.files(file.path(indir),pattern = c(".csv"))
  ids =str_split_fixed(test_files, "_", 3)[,2] #unique(sub("\\_.*", "", test_files));#ids = ids[ids!='RC3']
  #ids_patterns =c(ids,'RC3_01','RC3_02','RC3_03')
  ids_patterns = paste0(ids,'_')
  mdata =list()
  for (i in 1:length(ids)){
    sfiles <-list.files(file.path(indir),pattern = ids_patterns[i])
    date <- str_extract(sfiles, "[0-9]{4}-[0-9]{2}-[0-9]{2}")[1]
    sdata<-data.frame() #sdata <-list() #
    for (j in 1:length(sfiles)){
      lines <- readLines(file.path(indir,sfiles[j]))
      sidx<-max(grep("DateTime", lines))#[which(grep("Eureka_Manta", lines)<50)])
      #data <- read.delim(file.path(indir,date,location,file),skip = sidx-1,sep = ",",header = TRUE)
      data <- read.csv(file.path(indir,sfiles[j]),skip = sidx-1,header = TRUE)
      names(data) <- sub("\\..*", "",names(data))
      data <-na.omit(data)
      data$DateTime  <- as.POSIXct(trimws(data$DateTime, "l"),"%Y-%m-%d %H:%M:%S",tz='GMT')
      # int <- interval(min(data$Time)+minutes(8),max(data$Time)-minutes(5))
      # data <- data[data$Time%within%int,]
      sdata<-rbind(sdata,data)  #sdata[[j]] <-data ##
    }
    if (nrow(sdata)>0){
      #if (length(sdata)>0){
      cols = c('Temperature',"Dissolved_Oxygen","Dissolved_Oxygen_Saturation")
      lcols = c('Temperature',"Dissolved_Oxygen","Dissolved_Oxygen_Saturation")
      ylabs = c('Temperature (deg C)','Dissolved Oxygen (mg/l)',"Dissolved Oxygen Saturation(%)")
      png(file.path(plot_dir, #
                    paste0('TsPlot_',ids[i],'_',unique(sdata$Location),'_',date,'.png')), #'Testdata_',
          width = 5, height = 6, units = 'in', res = 300)
      par(mfrow=c(3,1),mgp=c(2,1,0),mar=c(2.1,3.1,2,8.5))
      ### plot with linear fit
      rd1=MDfitplot(sdata,col= cols[1],ylab=ylabs[1])
      rd2=MDfitplot(sdata,col= cols[2],ylab=ylabs[2])
      rd3=MDfitplot(sdata,col= cols[3],ylab=ylabs[3])
      dev.off()
      ## save lm fit summary data
      # rd1=rfit_summary(sdata,cols[1],ylab=ylabs[1])
      # rd2=rfit_summary(sdata,cols[2],ylab=ylabs[2])
      # rd3=rfit_summary(sdata,cols[3],ylab=ylabs[3])
      rds = mapply(c, rd1, rd2,rd3, SIMPLIFY=FALSE)
      d <- data.frame(matrix(Reduce(c,rds), nrow = 1, byrow = TRUE))
      cnames<- apply(expand.grid(c(1,2,3), lcols), 1, function(x) paste(x[2], x[1], sep="_")) 
      mnames<-c('Slope','Rsquare','RMSE','NRMSE','Mean','Range')
      cnames<- apply(expand.grid(cnames, mnames), 1, function(x) paste(x[1], x[2], sep="_"))
      names(d)<- cnames
      d['Site_ID']<- unique(sdata$Site_ID);d['Date']<-date;
      #d['Metrics']<- c('lm_slope','R_square','RMSE','Mean','Range')
      d<- d[c('Site_ID','Date',
              grep('Temperature',names(d),value=TRUE),
              grep(paste(c('Dissolved_Oxygen_1','Dissolved_Oxygen_2','Dissolved_Oxygen_3'), collapse = "|"),names(d),value=TRUE),
              grep('Dissolved_Oxygen_Saturation',names(d),value=TRUE))
      ]
      mdata[[i]] <-d
    }
  }
  df <- Reduce(rbind, mdata)
  #outdir = 'C:/Users/linx882/XLin/automation of respiration calculations'
  write.csv(df,file.path(summary_dir,paste0('RC2_Minidot_lm_summary','.csv')),
            na = "NA",row.names = FALSE)
  return(df)
}

##########################################################################
## violin plot to check the distribution of variables for Minidot data
MD_vplot<-function(indir,plot_dir,summary_dir){
  ## indir: : directory where test data are saved
  ## plot_dir: directory to save the generated plots
  ## summary_dir: directory to save lm regression summary data 
  test_files = list.files(file.path(indir),pattern = c(".csv"))
  ids = unique(sub("\\_.*", "", test_files));ids = ids[ids!='RC3']
  ids_patterns =c(ids,'RC3_01','RC3_02','RC3_03')
  ids_patterns = paste0(ids_patterns,'_')
  mdata =list()
  for (i in 1:length(ids)){
    sfiles <-list.files(file.path(indir),pattern = ids_patterns[i])
    date <- str_extract(sfiles, "[0-9]{4}-[0-9]{2}-[0-9]{2}")[1]
    sdata<-data.frame() #sdata <-list() #
    for (j in 1:length(sfiles)){
      lines <- readLines(file.path(indir,sfiles[j]))
      sidx<-max(grep("DateTime", lines))#[which(grep("Eureka_Manta", lines)<50)])
      #data <- read.delim(file.path(indir,date,location,file),skip = sidx-1,sep = ",",header = TRUE)
      data <- read.csv(file.path(indir,sfiles[j]),skip = sidx-1,header = TRUE)
      names(data) <- sub("\\..*", "",names(data))
      data$DateTime  <- as.POSIXct(trimws(data$DateTime, "l"),"%Y-%m-%d %H:%M:%S",tz='GMT')
      # int <- interval(min(data$Time)+minutes(8),max(data$Time)-minutes(5))
      # data <- data[data$Time%within%int,]
      sdata<-rbind(sdata,data)  #sdata[[j]] <-data ##
    }
    if (nrow(sdata)>0){
      #if (length(sdata)>0){
      cols = c('Temperature','Dissolved_Oxygen',"Dissolved_Oxygen_Saturation")
      ylabs = c('Temperature (deg C)','Dissolved Oxygen (mg/l)',"Dissolved Oxygen Saturation(%)")
      sdata$Series_ID<-as.factor(sdata$Series_ID)
      SID = unique(sdata$Series_ID)
      
      png(file.path(plot_dir, #
                    paste0('VPlot_',ids[i],'_',unique(sdata$Location),'_',date,'.png')), #'Testdata_',
          width = 5, height = 6, units = 'in', res = 300)
      par(mfrow=c(3,1),mgp=c(2,1,0),mar=c(3.1,3.1,1,8.5))
      for (k in 1:length(cols)){
        idx = which(colnames(sdata)==cols[k])
        if(length(idx)==1){
          vioplot(sdata[,idx]~sdata$Series_ID,names=SID,
                  ylab = ylabs[k], xlab = paste0(ids[i],'_',unique(sdata$Location)),
                  lineCol = "red",     # Color of the line
                  colMed = "red",col="lightgreen")
          #stripchart(sdata[,idx],method = "jitter",pch = 19,vertical = TRUE,col = "gold",add=TRUE)
        }
      }
      dev.off()
    }
  }
}

