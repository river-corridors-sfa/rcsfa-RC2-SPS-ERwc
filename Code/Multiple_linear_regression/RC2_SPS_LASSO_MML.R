# Code to run LASSO on SPS ERwc data 


# Libraries ---------------------------------------------------------------
library(tidyverse)
library(dplyr)
library(corrplot)
#library(ggpubr)
#library(ggpmisc)
# library(factoextra)
# library(stringr)
library(glmnet)
# library(magick)


# Working Directory -------------------------------------------------------
current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd("../..")
getwd()

rm(list=ls());graphics.off()

# Functions ---------------------------------------------------------------

cube_root <- function(x) sign(x) * (abs(x))^(1/3)

# function for pearson corr matrix
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = (cor(x, y, method = c("pearson")))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) {cex.cor <- 0.8/strwidth(txt)} else {cex = cex.cor}
  text(0.5, 0.5, txt, cex = cex.cor * (1 + r)/1)
  
  # if(missing(cex.cor)) {cex <- 1.2/strwidth(txt)} else {cex = cex.cor}
  # text(0.5, 0.5, txt, cex = cex * sin(sqrt(abs(r))))
  
  test <- cor.test(x,y)
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " "))
  #text(0.5, 0.5, txt, cex = cex * r)
  text(.5, .8, Signif, cex=cex, col=2)
  
}


# Read in Data ------------------------------------------------------------

mean_erwc = read.csv("./Data/Multiple_linear_regression/ERwc_Mean.csv") %>% 
  select(-X)

## are these the only variables I want?
geo = read.csv("./Data/Multiple_linear_regression/v2_RCSFA_Extracted_Geospatial_Data_2023-06-21 (1).csv") %>% 
  select(c(site, streamorde, totdasqkm)) %>% 
  dplyr::rename(Site_ID = site)

npoc_tn = read.csv("./Data/Multiple_linear_regression/v2_SPS_NPOC_TN.csv", skip = 2) %>% 
  filter(grepl("SPS", Sample_Name)) %>% 
  select(c(Sample_Name, X00681_NPOC_mg_per_L_as_C, X00602_TN_mg_per_L_as_N)) %>% 
  rename(NPOC = X00681_NPOC_mg_per_L_as_C) %>% 
  rename(TN = X00602_TN_mg_per_L_as_N) %>% 
  mutate(TN = if_else(grepl("Below", TN), as.numeric(.035), as.numeric(TN))) %>% 
  mutate(NPOC = as.numeric(NPOC)) %>% 
  separate(Sample_Name, c("Parent", "Rep"), sep = "-") %>% 
  group_by(Parent) %>% 
  #mutate(cv_npoc = sd(NPOC)/mean(NPOC)) %>% 
  mutate(cv_tn = sd(TN)/mean(TN))

npoc_tn$flag <- NA

#make final data frame
npoc_tn_final <- as.data.frame(matrix(NA, ncol = 7, nrow = 1))

colnames(npoc_tn_final) = c("tn.temp", "Parent", "Rep",  "NPOC",  "cv_tn", "cv_tn_rem", "flag")

#Number of unique parent IDs with replicates to loop through
unique.samples = unique(npoc_tn$Parent)

for (i in 1:length(unique.samples)) {
  
  ## Subset replicates
  data_subset = subset(npoc_tn, npoc_tn$Parent == unique.samples[i])
  
  ## Pull out TN percent values
  tn.temp = as.numeric(data_subset$TN)
  
  ## Calculate standard deviation, average, and coefficient of variation of rates
  tn.temp.sd <- sd(tn.temp)
  tn.temp.mean <- mean(tn.temp)
  CV = abs(tn.temp.sd/tn.temp.mean)
  
  #looping to get 3 best samples
  for (sample.reduction in 1:3)  {
    
    if (length(tn.temp) > 2) {
      
      dist.temp = as.matrix(abs(dist(tn.temp)))
      dist.comp = numeric() 
      
      for(tn.now in 1:ncol(dist.temp)) {
        
        dist.comp = rbind(dist.comp,c(tn.now,sum(dist.temp[,tn.now])))
        
      }
      
      dist.comp[,2] = as.numeric(dist.comp[,2])
      tn.temp = tn.temp[-which.max(dist.comp[,2])]
      
      tn.temp.sd <- sd(tn.temp)
      tn.temp.mean <- mean(tn.temp)
      tn.temp.cv <- abs(tn.temp.sd/tn.temp.mean)
      CV = tn.temp.cv
      tn.temp.range <- max(tn.temp) - min(tn.temp)
      range = tn.temp.range
      
    } 
  }
  
  if (length(tn.temp) >= 2) {
    
    tn.combined <- as.data.frame(tn.temp)
    
    tn.removed <- merge(tn.combined, data_subset, by.x = "tn.temp", by.y = "TN", all.x = TRUE)
    
    tn.removed <- tn.removed[!duplicated(tn.removed$Rep), ]
    
    tn.removed$cv_tn_rem = as.numeric(abs(sd(tn.temp)/mean(tn.temp)))
    

    
  }
  
  npoc_tn_final = rbind(tn.removed, npoc_tn_final)
  
  rm('tn.temp')
}

mean_npoc_tn_clean = npoc_tn_final %>% 
  rename(TN = tn.temp) %>% 
  group_by(Parent) %>% 
  mutate(Mean_NPOC = mean(NPOC)) %>% 
  mutate(Mean_TN = mean(TN)) %>% 
  distinct(Parent, .keep_all = T) %>% 
  select(c(Parent, Mean_NPOC, Mean_TN)) %>% 
  mutate(Parent = str_replace(Parent, "_OCN", ""))

## Try with published NO3 data

ions = read.csv("./Data/Multiple_linear_regression/v2_SPS_Ions.csv", skip = 2) %>% 
  filter(grepl("SPS", Sample_Name)) %>% 
  select(c(Sample_Name, X71851_NO3_mg_per_L_as_NO3)) 

cv_ions = ions %>% 
  separate(Sample_Name, c("Sample_Name", "Rep"), sep = "-") %>% 
  rename(NO3 = X71851_NO3_mg_per_L_as_NO3) %>% 
  mutate(NO3 = as.numeric(NO3)) %>% 
  filter(NO3 != -9999) %>% 
  group_by(Sample_Name) %>% 
  mutate(cv_no3 = sd(NO3)/mean(NO3)) 

#make final data frame
no3_final <- as.data.frame(matrix(NA, ncol = 5, nrow = 1))

colnames(no3_final) = c("no3.temp", "Sample_Name", "Rep",  "cv_no3", "cv_no3_rem")

#Number of unique parent IDs with replicates to loop through
unique.samples = unique(cv_ions$Sample_Name)

for (i in 1:length(unique.samples)) {
  
  ## Subset replicates
  data_subset = subset(cv_ions, cv_ions$Sample_Name == unique.samples[i])
  
  ## Pull out TN percent values
  no3.temp = as.numeric(data_subset$NO3)
  
  ## Calculate standard deviation, average, and coefficient of variation of rates
  no3.temp.sd <- sd(no3.temp)
  no3.temp.mean <- mean(no3.temp)
  CV = abs(no3.temp.sd/no3.temp.mean)
  
  #looping to get 3 best samples
  for (sample.reduction in 1:3)  {
    
    if (length(no3.temp) > 2) {
      
      dist.temp = as.matrix(abs(dist(no3.temp)))
      dist.comp = numeric() 
      
      for(no3.now in 1:ncol(dist.temp)) {
        
        dist.comp = rbind(dist.comp,c(no3.now,sum(dist.temp[,no3.now])))
        
      }
      
      dist.comp[,2] = as.numeric(dist.comp[,2])
      no3.temp = no3.temp[-which.max(dist.comp[,2])]
      
      no3.temp.sd <- sd(no3.temp)
      no3.temp.mean <- mean(no3.temp)
      no3.temp.cv <- abs(no3.temp.sd/no3.temp.mean)
      CV = no3.temp.cv
      no3.temp.range <- max(no3.temp) - min(no3.temp)
      range = no3.temp.range
      
    } 
  }
  
  if (length(no3.temp) >= 1) {
    
    no3.combined <- as.data.frame(no3.temp)
    
    no3.removed <- merge(no3.combined, data_subset, by.x = "no3.temp", by.y = "NO3", all.x = TRUE)
    
    no3.removed <- no3.removed[!duplicated(no3.removed$Rep), ]
    
    no3.removed$cv_no3_rem = as.numeric(abs(sd(no3.temp)/mean(no3.temp)))
    
    
    
  }
  
  no3_final = rbind(no3.removed, no3_final)
  
  rm('no3.temp')
}

mean_no3_clean = no3_final %>% 
  rename(NO3 = no3.temp) %>% 
  group_by(Sample_Name) %>% 
  distinct(Rep, .keep_all = T) %>% 
  mutate(Mean_NO3 = mean(NO3)) %>% 
  ungroup() %>% 
  distinct(Sample_Name, .keep_all = T) %>% 
  mutate(Parent = str_replace(Sample_Name, "_ION", "")) %>% 
  select(c(Parent, Mean_NO3))

# Clean NO3 that we don't trust (below 0.5 mg/L standard) - give all these a value of 0.25 mg/L  

lod_no3 = cv_ions %>% 
  group_by(Sample_Name) %>% 
  mutate(NO3_lod = ifelse(NO3 < 0.5, 0.25, NO3)) %>% 
  mutate(cv_no3_lod = sd(NO3_lod)/mean(NO3_lod)) %>% 
  filter(NO3_lod != 0.52) %>% # remove outlier form SPS_0028 - only one with LODs
  mutate(Mean_NO3_lod = mean(NO3_lod)) %>% 
  ungroup() %>% 
  distinct(Sample_Name, .keep_all = T) %>% 
  mutate(Parent = str_replace(Sample_Name, "_ION", "")) %>% 
  select(c(Parent, Mean_NO3_lod))

#Check LOD
tss = read.csv("./Data/Multiple_linear_regression/SPS_TSS.csv", skip = 2) %>% 
  filter(grepl("SPS", Sample_Name)) %>% 
  select(c(Sample_Name, X00530_TSS_mg_per_L)) %>% 
  rename(TSS = X00530_TSS_mg_per_L) %>% 
  mutate(TSS = if_else(grepl("Below", TSS), as.numeric(.743), as.numeric(TSS))) %>% 
  separate(Sample_Name, c("Parent","Rep"), sep = "-") %>% 
  mutate(Parent = str_replace(Parent, "_TSS", "")) %>% 
  select(c(Parent, TSS))

dic = read.csv("C:/Users/laan208/OneDrive - PNNL/Shared Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/RC2/Boye_Files/SPS/SPS_Water_DIC_Boye_2024-09-10.csv", skip = 2) %>% 
  filter(grepl("SPS", Sample_Name)) %>% 
  select(c(Sample_Name, X00691_DIC_mg_per_L_as_C)) %>% 
  rename(DIC = X00691_DIC_mg_per_L_as_C)

mean_dic = dic %>% 
  mutate(DIC = as.numeric(DIC)) %>% 
  separate(Sample_Name, c("Parent", "Rep"), sep = "-") %>% 
  group_by(Parent) %>% 
  mutate(mean_DIC = mean(DIC)) %>% 
  select(c(Parent, mean_DIC)) %>% 
  distinct(Parent, .keep_all = T) %>% 
  mutate(Parent = str_replace(Parent, "_DIC", ""))

sample = left_join(mean_dic, mean_npoc_tn_clean) %>% 
  left_join(tss) %>% 
  left_join(mean_no3_clean) %>% 
  left_join(lod_no3) %>% 
  rename(Sample_Name = Parent)


# need to get DIC from VGC, check Ions from Sophia, and pull together in summary file code, add ultrameter water chemistry to this? add manta data to this? 
# sample = read.csv("./Data/Multiple_linear_regression/Summary_Not_Cleaned.csv") %>% 
#   mutate(Sample_Name = str_remove(Parent, "_MEAN")) %>%   select(-c(X, Parent, X00000_NH4_mg_per_L_as_NH4_mean, X01130_Li_mg_per_L_mean, X00653_PO4_mg_per_L_as_PO4_mean, X71856_NO2_mg_per_L_as_NO2_mean)) %>% 
#   rename_with(~ str_sub(., 8), starts_with("X")) %>% 
#   mutate(across(everything(), ~if_else(. == -9999, NA, .))) %>% 
#   select(c(Sample_Name, TSS_mg_per_L, NPOC_mg_per_L_as_C_mean, DIC_mean))
#   
mapping = read.csv("./Data/Multiple_linear_regression/v2_SPS_Sensor_Field_Metadata.csv") %>% 
  select(c(Site_ID, Sample_Name))

om <- read.csv(file.path("./Data/OM_transformation_analysis/SPS_Total_and_Normalized_Transformations_01-03-23.csv")) 

# Merge Data --------------------------------------------------------------

all_data = left_join(mean_erwc, mapping, by = "Site_ID") %>% 
  left_join(geo, by = "Site_ID") %>% 
  full_join(sample, by = "Sample_Name") %>% 
  left_join(om, by = "Sample_Name") %>% 
  select(-c(Sample_Name))

## Shorten Names 

new_names = c(ERwc = "Mean_ERwc", Temp = "Mean_Temp", StrOrd = "streamorde", TotDr = "totdasqkm", Trans = "Total_Number_of_Transformations", Peaks = "Number_of_Peaks", NormTrans = "Normalized_Transformations"#, TSS = "TSS_mg_per_L", DIC = "DIC_mean", NPOC = "NPOC_mg_per_L_as_C_mean"#, TN = "TN_mg_per_L_as_N_mean"#, Br = #"Br_mg_per_L_mean", Ca = #"Ca_mg_per_L_mean", Cl = #"Cl_mg_per_L_mean", Fl = #"F_mg_per_L_mean", Mg = #"Mg_mg_per_L_mean", NO3 = #"NO3_mg_per_L_as_NO3_mean", K = #"K_mg_per_L_mean", Na = #"Na_mg_per_L_mean", SO4 = #"SO4_mg_per_L_as_SO4_mean"
                )

new_data <- all_data %>% 
  rename(!!!new_names) %>% 
  column_to_rownames("Site_ID")

## Look at histograms

long_data = new_data %>% 
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

ggplot() + 
  geom_histogram(long_data, mapping = aes(x = value)) + 
  facet_wrap(~ variable, scales = "free") +
  theme_minimal()
  
# Transform data ----------------------------------------------------------

#decide how you want to do this, eg, cube or log transform

#should everything be transformed (eg Peaks, stream order, Transformations?)

# log gives a lot of NAs even with + 1
log_data = new_data %>% 
  mutate_if(is.numeric, log10)

log_data_plus_one = new_data %>% 
  mutate_if(is.numeric, ~ log10(. + 1))

# cube root allows you to keep sign
# remove NA values from analysis (might change later)
cube_data = new_data %>% 
 mutate(across(where(is.numeric), cube_root)) %>% 
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x)) %>% 
  drop_na()

long_cube_data = cube_data %>% 
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

ggplot() + 
  geom_histogram(long_cube_data, mapping = aes(x = value)) + 
  facet_wrap(~ variable, scales = "free") +
  theme_minimal()


# Check Co-Linearity ------------------------------------------------------
## Pearson without cube transformation
scale_data = as.data.frame(scale(new_data))%>% 
  rename_with(where(is.numeric), .fn = ~ paste0("scale_", .x)) %>% 
  drop_na()

scale_pearson <- cor(scale_data, method = "pearson")

png(file = paste0("./Figures/", as.character(Sys.Date()),"_Scale_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)


corrplot(scale_pearson,type = "upper", method = "number", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25,  title = "Effect Samples Pearson Correlation")

dev.off()


## Scale data before it goes into correlation matrix

scale_cube_data = as.data.frame(scale(cube_data))%>% 
  rename_with(where(is.numeric), .fn = ~ paste0("scale_", .x))

scale_cube_pearson <- cor(scale_cube_data, method = "pearson")

png(file = paste0("./Figures/", as.character(Sys.Date()),"_Cube_Scale_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)


corrplot(scale_cube_pearson,type = "upper", method = "number", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25,  title = "Effect Samples Pearson Correlation")

dev.off()

## Keep just TN ####
pearson_df <- as.data.frame(scale_cube_pearson)

row_names_pearson <- rownames(pearson_df)

pearson_df$Variable <- row_names_pearson

# Melt the dataframe for plotting
pearson_melted <- reshape2::melt(pearson_df, id.vars = "Variable") %>% 
  filter(value != 1) %>% 
  mutate(value = abs(value)) %>% # do this so it removes in order, and doesn't leave out high negative correlations
  filter(!grepl("ERwc", Variable)) %>% # remove ERwc, don't want it to be removed %>% 
  filter(!grepl("Mean_NO3|Mean_NO3_lod",Variable) & !grepl("Mean_NO3|Mean_NO3_lod", variable))

# pull out erwc correlations only
erwc_melted <- pearson_melted %>% 
  filter(grepl("ERwc", variable)) 

choose_melted <- pearson_melted %>% 
  filter(!grepl("ERwc", variable)) %>%
  #distinct(value, .keep_all = TRUE) %>% 
  left_join(erwc_melted, by = "Variable") %>% 
  rename(Variable_1 = Variable) %>% 
  rename(Variable_2 = variable.x) %>% 
  rename(Correlation = value.x) %>% 
  rename(Variable_1_ERwc_Correlation = value.y) %>% 
  select(-c(variable.y)) %>% 
  left_join(erwc_melted, by = c("Variable_2" = "Variable")) %>% 
  rename(Variable_2_ERwc_Correlation = value) %>% 
  select(-c(variable))

loop_melt = choose_melted %>% 
  arrange(desc(Correlation))

# Pearson correlation coefficient to remove above
correlation = 0.7

## Start loop to remove highly correlated (> 0.5)
erwc_filter = function(loop_melt) {
  
  rows_to_keep = rep(TRUE, nrow(loop_melt))
  
  for (i in seq_len(nrow(loop_melt))) {
    
    if (!rows_to_keep[i]) next
    
    row = loop_melt[i, ]
    
    if (row$Correlation < correlation) next
    
    if(row$Variable_1_ERwc_Correlation >= row$Variable_2_ERwc_Correlation) {
      
      var_to_keep = row$Variable_1
      var_to_remove = row$Variable_2
      
    } else {
      
      var_to_keep = row$Variable_2
      var_to_remove = row$Variable_1
      
    }
    
    loop_melt$Variable_to_Keep[i] = var_to_keep
    loop_melt$Variable_to_Remove[i] = var_to_remove
    
    for (j in seq(i + 1, nrow(loop_melt))) {
      
      if(loop_melt$Variable_1[j] == var_to_remove || loop_melt$Variable_2[j] == var_to_remove) {
        
        rows_to_keep[j] = FALSE
        
      }
      
    }
    
    
  }
  
  return(loop_melt[rows_to_keep, ])
  
}

filtered_data = erwc_filter(loop_melt) 

# pull out variables to remove
removed_variables = filtered_data %>% 
  distinct(Variable_to_Remove)

# pull out all variables 
all_variables = erwc_melted %>% 
  select(c(Variable))

# remove variables from all variables to get variables to keep for LASSO 
kept_variables = erwc_melted[!(erwc_melted$Variable %in% removed_variables$Variable_to_Remove), ]


col_to_keep = unique(kept_variables$Variable)
col_to_keep = c(col_to_keep, "scale_cube_ERwc")

scale_cube_variables = scale_cube_data[, col_to_keep, drop = FALSE]


# Start LASSO -------------------------------------------------------------

## LASSO with Correlation Matrix Selected Variables 
set.seed(42)
## Set response variable (Cube_Effect_Size) and scale
yvar <- data.matrix(scale_cube_variables$scale_cube_ERwc)
mean(yvar)
sd(yvar)

## Set predictor variables and scale
exclude_col = "scale_cube_ERwc"

x_cube_variables = as.data.frame(scale_cube_variables[, !(names(scale_cube_variables) %in% exclude_col)])
#mean(x_cube_variables$scale_cube_Temp)
#sd(x_cube_variables$scale_cube_Temp)

xvars <- data.matrix(x_cube_variables)

lasso = cv.glmnet(xvars, yvar, alpha = 1, nfolds = 5,
                  ,standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                  #,standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                  # , standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

best_lambda <- lasso$lambda.min
best_lambda

plot(lasso)

best_lasso_model <- glmnet(xvars, yvar, alpha = 1, lambda = best_lambda, family = "gaussian",
                           , standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                           #  , standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                           #, standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

#lasso_coefs = coef(best_lasso_model)
coef(best_lasso_model)

yvar_predict <- predict(best_lasso_model, s = best_lambda, newx = xvars)

sst <- sum((yvar - mean(yvar))^2)
sse <- sum((yvar_predict - yvar)^2)

rsq = 1 - sse/sst

rsq

## Using Published NO3 ####

# Melt the dataframe for plotting
pearson_melted <- reshape2::melt(pearson_df, id.vars = "Variable") %>% 
  filter(value != 1) %>% 
  mutate(value = abs(value)) %>% # do this so it removes in order, and doesn't leave out high negative correlations
  filter(!grepl("ERwc", Variable)) %>% # remove ERwc, don't want it to be removed %>% 
  filter(!grepl("Mean_TN|Mean_NO3_lod",Variable) & !grepl("Mean_TN|Mean_NO3_lod", variable))

# pull out erwc correlations only
erwc_melted <- pearson_melted %>% 
  filter(grepl("ERwc", variable)) 

choose_melted <- pearson_melted %>% 
  filter(!grepl("ERwc", variable)) %>%
  #distinct(value, .keep_all = TRUE) %>% 
  left_join(erwc_melted, by = "Variable") %>% 
  rename(Variable_1 = Variable) %>% 
  rename(Variable_2 = variable.x) %>% 
  rename(Correlation = value.x) %>% 
  rename(Variable_1_ERwc_Correlation = value.y) %>% 
  select(-c(variable.y)) %>% 
  left_join(erwc_melted, by = c("Variable_2" = "Variable")) %>% 
  rename(Variable_2_ERwc_Correlation = value) %>% 
  select(-c(variable))

loop_melt = choose_melted %>% 
  arrange(desc(Correlation))

# Pearson correlation coefficient to remove above
correlation = 0.7

## Start loop to remove highly correlated (> 0.5)
erwc_filter = function(loop_melt) {
  
  rows_to_keep = rep(TRUE, nrow(loop_melt))
  
  for (i in seq_len(nrow(loop_melt))) {
    
    if (!rows_to_keep[i]) next
    
    row = loop_melt[i, ]
    
    if (row$Correlation < correlation) next
    
    if(row$Variable_1_ERwc_Correlation >= row$Variable_2_ERwc_Correlation) {
      
      var_to_keep = row$Variable_1
      var_to_remove = row$Variable_2
      
    } else {
      
      var_to_keep = row$Variable_2
      var_to_remove = row$Variable_1
      
    }
    
    loop_melt$Variable_to_Keep[i] = var_to_keep
    loop_melt$Variable_to_Remove[i] = var_to_remove
    
    for (j in seq(i + 1, nrow(loop_melt))) {
      
      if(loop_melt$Variable_1[j] == var_to_remove || loop_melt$Variable_2[j] == var_to_remove) {
        
        rows_to_keep[j] = FALSE
        
      }
      
    }
    
    
  }
  
  return(loop_melt[rows_to_keep, ])
  
}

filtered_data = erwc_filter(loop_melt) 

# pull out variables to remove
removed_variables = filtered_data %>% 
  distinct(Variable_to_Remove)

# pull out all variables 
all_variables = erwc_melted %>% 
  select(c(Variable))

# remove variables from all variables to get variables to keep for LASSO 
kept_variables = erwc_melted[!(erwc_melted$Variable %in% removed_variables$Variable_to_Remove), ]


col_to_keep = unique(kept_variables$Variable)
col_to_keep = c(col_to_keep, "scale_cube_ERwc")

scale_cube_variables = scale_cube_data[, col_to_keep, drop = FALSE]


# Start LASSO -------------------------------------------------------------

## LASSO with Correlation Matrix Selected Variables 
set.seed(42)
## Set response variable (Cube_Effect_Size) and scale
yvar <- data.matrix(scale_cube_variables$scale_cube_ERwc)
mean(yvar)
sd(yvar)

## Set predictor variables and scale
exclude_col = "scale_cube_ERwc"

x_cube_variables = as.data.frame(scale_cube_variables[, !(names(scale_cube_variables) %in% exclude_col)])
#mean(x_cube_variables$scale_cube_Temp)
#sd(x_cube_variables$scale_cube_Temp)

xvars <- data.matrix(x_cube_variables)

lasso = cv.glmnet(xvars, yvar, alpha = 1, nfolds = 5,
                  ,standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                  #,standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                  # , standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

best_lambda <- lasso$lambda.min
best_lambda

plot(lasso)

best_lasso_model <- glmnet(xvars, yvar, alpha = 1, lambda = best_lambda, family = "gaussian",
                           , standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                           #  , standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                           #, standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

#lasso_coefs = coef(best_lasso_model)
coef(best_lasso_model)

yvar_predict <- predict(best_lasso_model, s = best_lambda, newx = xvars)

sst <- sum((yvar - mean(yvar))^2)
sse <- sum((yvar_predict - yvar)^2)

rsq = 1 - sse/sst

rsq

## Using 0.25 NO3 ####

# Melt the dataframe for plotting
pearson_melted <- reshape2::melt(pearson_df, id.vars = "Variable") %>% 
  filter(value != 1) %>% 
  mutate(value = abs(value)) %>% # do this so it removes in order, and doesn't leave out high negative correlations
  filter(!grepl("ERwc", Variable)) %>% # remove ERwc, don't want it to be removed %>% 
  filter(!grepl("Mean_TN",Variable) & !grepl("Mean_TN", variable)) %>% 
  filter(Variable != "scale_cube_Mean_NO3") %>% 
  filter(variable != "scale_cube_Mean_NO3")

# pull out erwc correlations only
erwc_melted <- pearson_melted %>% 
  filter(grepl("ERwc", variable)) 

choose_melted <- pearson_melted %>% 
  filter(!grepl("ERwc", variable)) %>%
  #distinct(value, .keep_all = TRUE) %>% 
  left_join(erwc_melted, by = "Variable") %>% 
  rename(Variable_1 = Variable) %>% 
  rename(Variable_2 = variable.x) %>% 
  rename(Correlation = value.x) %>% 
  rename(Variable_1_ERwc_Correlation = value.y) %>% 
  select(-c(variable.y)) %>% 
  left_join(erwc_melted, by = c("Variable_2" = "Variable")) %>% 
  rename(Variable_2_ERwc_Correlation = value) %>% 
  select(-c(variable))

loop_melt = choose_melted %>% 
  arrange(desc(Correlation))

# Pearson correlation coefficient to remove above
correlation = 0.7

## Start loop to remove highly correlated (> 0.5)
erwc_filter = function(loop_melt) {
  
  rows_to_keep = rep(TRUE, nrow(loop_melt))
  
  for (i in seq_len(nrow(loop_melt))) {
    
    if (!rows_to_keep[i]) next
    
    row = loop_melt[i, ]
    
    if (row$Correlation < correlation) next
    
    if(row$Variable_1_ERwc_Correlation >= row$Variable_2_ERwc_Correlation) {
      
      var_to_keep = row$Variable_1
      var_to_remove = row$Variable_2
      
    } else {
      
      var_to_keep = row$Variable_2
      var_to_remove = row$Variable_1
      
    }
    
    loop_melt$Variable_to_Keep[i] = var_to_keep
    loop_melt$Variable_to_Remove[i] = var_to_remove
    
    for (j in seq(i + 1, nrow(loop_melt))) {
      
      if(loop_melt$Variable_1[j] == var_to_remove || loop_melt$Variable_2[j] == var_to_remove) {
        
        rows_to_keep[j] = FALSE
        
      }
      
    }
    
    
  }
  
  return(loop_melt[rows_to_keep, ])
  
}

filtered_data = erwc_filter(loop_melt) 

# pull out variables to remove
removed_variables = filtered_data %>% 
  distinct(Variable_to_Remove)

# pull out all variables 
all_variables = erwc_melted %>% 
  select(c(Variable))

# remove variables from all variables to get variables to keep for LASSO 
kept_variables = erwc_melted[!(erwc_melted$Variable %in% removed_variables$Variable_to_Remove), ]


col_to_keep = unique(kept_variables$Variable)
col_to_keep = c(col_to_keep, "scale_cube_ERwc")

scale_cube_variables = scale_cube_data[, col_to_keep, drop = FALSE]


# Start LASSO -------------------------------------------------------------

## LASSO with Correlation Matrix Selected Variables 
set.seed(42)
## Set response variable (Cube_Effect_Size) and scale
yvar <- data.matrix(scale_cube_variables$scale_cube_ERwc)
mean(yvar)
sd(yvar)

## Set predictor variables and scale
exclude_col = "scale_cube_ERwc"

x_cube_variables = as.data.frame(scale_cube_variables[, !(names(scale_cube_variables) %in% exclude_col)])
#mean(x_cube_variables$scale_cube_Temp)
#sd(x_cube_variables$scale_cube_Temp)

xvars <- data.matrix(x_cube_variables)

lasso = cv.glmnet(xvars, yvar, alpha = 1, nfolds = 5,
                  ,standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                  #,standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                  # , standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

best_lambda <- lasso$lambda.min
best_lambda

plot(lasso)

best_lasso_model <- glmnet(xvars, yvar, alpha = 1, lambda = best_lambda, family = "gaussian",
                           , standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                           #  , standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                           #, standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

#lasso_coefs = coef(best_lasso_model)
coef(best_lasso_model)

yvar_predict <- predict(best_lasso_model, s = best_lambda, newx = xvars)

sst <- sum((yvar - mean(yvar))^2)
sse <- sum((yvar_predict - yvar)^2)

rsq = 1 - sse/sst

rsq



