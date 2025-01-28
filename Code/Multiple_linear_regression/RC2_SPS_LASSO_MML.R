
# Code to run LASSO on SPS ERwc data for Table 2

# Makes Figures 3b, 4, S2

# Author: Maggi Laan (maggi.laan@gmail.com)

# Libraries ------------------------------------------------------------
library(tidyverse)
library(corrplot)
library(ggpubr)
library(ggpmisc)
library(glmnet)
library(readxl)
library(reshape2)
library(viridis)


# Working Directory ----------------------------------------------------
current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd("../..")
getwd()

rm(list=ls());graphics.off()

# Functions ------------------------------------------------------------

# Cube Root Transformation
cube_root <- function(x) sign(x) * (abs(x))^(1/3)

# Pearson corr matrix
pear.panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = (cor(x, y, method = c("pearson")))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) {cex.cor <- 0.8/strwidth(txt)} 
  text(0.5, 0.5, txt, cex = cex.cor * (1 + abs(r))/2)
  
  # if(missing(cex.cor)) {cex <- 1.2/strwidth(txt)} else {cex = cex.cor}
  # text(0.5, 0.5, txt, cex = cex * sin(sqrt(abs(r))))
  
  test <- cor.test(x,y)
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " "))
  #text(0.5, 0.5, txt, cex = cex * r)
  text(.5, .8, Signif, cex=cex.cor, col=2)
  
}

# Scatter Plots 
panel.smooth <- function(x, y) {
  points(x, y, pch = 19, col = rgb(0.1, 0.2, 0.5, alpha = 0.3))
  abline(lm(y ~ x), col = 'blue', lty = 2)
}

# Histograms
panel.hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1))
  
  h <- hist(x, plot = FALSE, breaks = "FD")
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  
  rect(breaks[-nB], 0, breaks[-1], y, col="grey", border="white", ...)
}

# Read in Data ---------------------------------------------------------

mean_erwc = read.csv("./Data/Multiple_linear_regression/ERwc_Mean.csv") %>% 
  select(-X)

## Pull out stream order and total drainage area
geo = read.csv("https://github.com/river-corridors-sfa/Geospatial_variables/raw/refs/heads/main/v2_RCSFA_Extracted_Geospatial_Data_2023-06-21.csv") %>% 
  select(c(site, streamorde, totdasqkm)) %>% 
  dplyr::rename(Site_ID = site)

# DIC Data
data = read.csv("./Data/Multiple_linear_regression/v3_SFA_SpatialStudy_2021_SampleData/v3_SPS_Water_Sample_Data_Summary.csv", skip = 2) %>% 
  filter(grepl("SPS", Sample_Name)) %>% 
  select(c(Sample_Name, Mean_00691_DIC_mg_per_L_as_C))
  
## DOC/TDN 
npoc_tn = read.csv("./Data/Multiple_linear_regression/v3_SFA_SpatialStudy_2021_SampleData/v3_SPS_Water_NPOC_TN.csv", skip = 2) %>% 
  filter(grepl("SPS", Sample_Name)) %>% 
  select(c(Sample_Name, X00681_NPOC_mg_per_L_as_C, X00602_TN_mg_per_L_as_N, Methods_Deviation)) %>% 
  rename(NPOC = X00681_NPOC_mg_per_L_as_C) %>% 
  rename(TN = X00602_TN_mg_per_L_as_N) %>% 
  mutate(TN = if_else(grepl("Below", TN), as.numeric(.035), as.numeric(TN)), missing = NA_real_) %>% #set samples below standard or LOD to half of LOD (0.035)
  mutate(NPOC = as.numeric(NPOC)) %>% 
  filter(!grepl("OUTLIER", Methods_Deviation))%>% 
  separate(Sample_Name, c("Parent", "Rep"), sep = "-") %>% 
  group_by(Parent) %>% 
  mutate(Mean_NPOC = mean(NPOC)) %>% 
  mutate(Mean_TN = mean(TN)) %>% 
  ungroup() %>% 
  distinct(Parent, .keep_all = T) %>% 
  mutate(Sample_Name = str_replace(Parent, "_OCN", "_Water")) %>% 
  select(c(Sample_Name, Mean_NPOC, Mean_TN))


## NO3, SO4, Cl
ions = read.csv("./Data/Multiple_linear_regression/v3_SFA_SpatialStudy_2021_SampleData/v3_SPS_Water_Ions.csv", skip = 2) %>% 
  filter(grepl("SPS", Sample_Name)) %>% 
  mutate(NO3_mg_per_L = if_else(grepl("Nitrate", X71851_NO3_mg_per_L_as_NO3), as.numeric(0.035), as.numeric(X71851_NO3_mg_per_L_as_NO3))) %>% #set samples below standard or LOD to half of LOD (0.035)
  mutate(Sample_Name = str_replace(Sample_Name, "ION", "Water")) %>% 
  separate(Sample_Name, c("Sample_Name", "Rep"), sep = "-") %>% 
  mutate(Cl_mg_per_L = as.numeric(X00940_Cl_mg_per_L)) %>% 
  mutate(SO4_mg_per_L = as.numeric(X00945_SO4_mg_per_L_as_SO4)) %>% 
  select(c(Sample_Name, NO3_mg_per_L, Cl_mg_per_L, SO4_mg_per_L)) 
  
##TSS
#Check LOD
tss = read.csv("./Data/Multiple_linear_regression/v3_SFA_SpatialStudy_2021_SampleData/v2_SPS_Water_TSS.csv", skip = 2) %>% 
  filter(grepl("SPS", Sample_Name)) %>% 
  select(c(Sample_Name, X00530_TSS_mg_per_L)) %>% 
  rename(TSS = X00530_TSS_mg_per_L) %>% 
  mutate(TSS = if_else(grepl("Below", TSS), as.numeric(.743), as.numeric(TSS))) %>% 
  mutate(TSS = if_else(TSS <= 0.24, 0.12, TSS)) %>% 
  separate(Sample_Name, c("Parent","Rep"), sep = "-") %>% 
  mutate(Sample_Name = str_replace(Parent, "_TSS", "_Water")) %>% 
  select(c(Sample_Name, TSS))

## Organic Matter
# check this is still good
om <- read.csv(file.path("./Data/OM_transformation_analysis/SPS_Total_and_Normalized_Transformations_01-03-23.csv")) 

# Join all sample data
sample = left_join(data, npoc_tn) %>% 
  left_join(ions) %>% 
  left_join(tss) %>% 
  rename(Mean_DIC = Mean_00691_DIC_mg_per_L_as_C) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "_Water", "")) %>% 
  mutate(Mean_DIC = as.numeric(Mean_DIC)) %>% 
  left_join(om) 

# add ultrameter water chemistry to this? add manta data to this? 

## Get Site ID and Sample Names to merge with geospatial data
mapping = read.csv("./Data/Multiple_linear_regression/v2_SPS_Sensor_Field_Metadata.csv") %>% 
  select(c(Site_ID, Sample_Name))


# Merge Data -----------------------------------------------------------

all_data = left_join(mean_erwc, mapping, by = "Site_ID") %>% 
  left_join(geo, by = "Site_ID") %>% 
  full_join(sample, by = "Sample_Name") %>% 
  select(-c(Sample_Name))


## Clean Data ####

## Shorten Names 
new_names = c(ERwc = "Mean_ERwc", Temp = "Mean_Temp", StrOrd = "streamorde", TotDr = "totdasqkm", Transformations = "Total_Number_of_Transformations", Peaks = "Number_of_Peaks", NormTrans = "Normalized_Transformations")

## Remove values > 0.5, which are biologically unrealistic. In this dataset, these are likely from diffusion processes as [DO] starts ~5. ERwc < 0.5 is difficult to distinguish from 0, so values are kept

new_data <- all_data %>% 
  rename(!!!new_names) %>% 
  column_to_rownames("Site_ID") %>%
  filter(ERwc < 0.5) %>% # removes S38 and S83
  drop_na() # S68 doesn't have geospatial data


# Data Visualization ------------------------------------------------------

## Look at histograms of data ####

new_data %>% 
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>% 
ggplot() + 
  geom_histogram(mapping = aes(x = value)) + 
  facet_wrap(~ variable, scales = "free") +
  theme_minimal()

## Pearson correlation before transformations ####

png(file = paste0("./Figures/", as.character(Sys.Date()),"_Pairs_Pearson_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)

pairs(new_data,
      lower.panel = panel.smooth, 
      upper.panel = pear.panel.cor, 
      diag.panel = panel.hist,
      labels = colnames(new_data),
      cex.labels = 0.8) 

dev.off()

# Cube Root Transform data ---------------------------------------------

# cube root allows you to keep sign

cube_data = new_data %>% 
 mutate(across(where(is.numeric), cube_root)) %>% 
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x)) 

## Look at histograms of cube data
cube_data %>% 
  rownames_to_column("Sample_Name") %>% 
  pivot_longer(!Sample_Name, names_to = "variable", values_to = "value") %>% ggplot() + 
  geom_histogram( mapping = aes(x = value)) + 
  facet_wrap(~ variable, scales = "free") +
  theme_minimal()


# Check Co-Linearity ---------------------------------------------------

## Pearson Correlation Matrix of Transformed Data ####

## Scale data before it goes into correlation matrix (doesn't change pearson values - scaling here for consistency before LASSO)

scale_cube_data = as.data.frame(scale(cube_data))%>% 
  rename_with(where(is.numeric), .fn = ~ paste0("scale_", .x))

# check that mean is 0 and sd is 1
round(mean(scale_cube_data$scale_cube_ERwc), 4)
sd(scale_cube_data$scale_cube_ERwc)

# Matrix of Pearson correlation values
scale_cube_pearson <- cor(scale_cube_data, method = "pearson")

png(file = paste0("./Figures/", as.character(Sys.Date()),"_Pairs_Scale_Cube_Pearson_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)

pear_data = scale_cube_data %>% 
  rename_all(~ sub("scale_cube_", "", .)) %>% 
  rename("Temperature" = "Temp") %>% 
  rename("Stream Order" = "StrOrd") %>% 
  rename("Total Drainage" = "TotDr") %>% 
  rename("DIC" = "Mean_DIC") %>% 
  rename("DOC" = "Mean_NPOC") %>% 
  rename("TDN" = "Mean_TN") %>% 
  rename("NO3" = "NO3_mg_per_L") %>% 
  rename("Cl" = "Cl_mg_per_L") %>% 
  rename("SO4" = "SO4_mg_per_L") %>% 
  rename("Norm. Transformations" = "NormTrans")
  

pairs(pear_data,
      lower.panel = panel.smooth, 
      upper.panel = pear.panel.cor, 
      diag.panel = panel.hist,
      labels = colnames(pear_data),
      cex.labels = 0.5) 

dev.off()

# Make pearson correlation matrix into df
pearson_df = as.data.frame(scale_cube_pearson)  %>% 
  rownames_to_column("Variable")

## Melt the Pearson dataframe for removing variables ####
pearson_melted <- reshape2::melt(pearson_df, id.vars = "Variable") %>% 
  filter(value != 1) %>% # remove self-correlations
  mutate(value = abs(value)) %>% # do this so it removes in order, and doesn't leave out high negative correlations
  filter(!grepl("ERwc", Variable)) %>% # remove ERwc from first column
  #filter(!grepl("Mean_TN",Variable) & !grepl("Mean_TN", variable)) %>% # Run twice, once keeping TN and once removing
  filter(!grepl("StrOrd", Variable) & !grepl("StrOrd", variable)) #try removing stream order to keep drainage area 

# Pull out ERwc correlations only
erwc_melted <- pearson_melted %>% 
  filter(grepl("ERwc", variable)) 

# Remove ERwc from this DF
loop_melt <- pearson_melted %>% 
  filter(!grepl("ERwc", variable)) %>%
  left_join(erwc_melted, by = "Variable") %>% 
  rename(Variable_1 = Variable) %>% 
  rename(Variable_2 = variable.x) %>% 
  rename(Correlation = value.x) %>% 
  rename(Variable_1_ERwc_Correlation = value.y) %>% 
  select(-c(variable.y)) %>% 
  left_join(erwc_melted, by = c("Variable_2" = "Variable")) %>% 
  rename(Variable_2_ERwc_Correlation = value) %>% 
  select(-c(variable)) %>% 
  arrange(desc(Correlation))

# Pearson correlation coefficient to remove above
correlation = 0.7

## Start loop to remove highly correlated (> 0.7)
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
kept_variables = erwc_melted %>% 
  filter(!(Variable %in% removed_variables$Variable_to_Remove))

col_to_keep = unique(kept_variables$Variable)
col_to_keep = c(col_to_keep, "scale_cube_ERwc")

scale_cube_variables = scale_cube_data[, col_to_keep, drop = FALSE]


# Start LASSO ----------------------------------------------------------

## Loop through LASSO to get average over a lot of seeds ####

num_seeds = 100
seeds = sample(1:500, num_seeds)

## Set response variable (ERwc)
yvar <- data.matrix(scale_cube_variables$scale_cube_ERwc)

# list for storing LASSO iterations
norm_coeffs = list()
lasso_coefs_pull = list()
r2_scores = numeric(num_seeds)

## Set predictor variables
exclude_col = "scale_cube_ERwc"

x_cube_variables = as.data.frame(scale_cube_variables[, !(names(scale_cube_variables) %in% exclude_col)])
#mean(x_cube_variables$scale_cube_Temp)
#sd(x_cube_variables$scale_cube_Temp)

xvars <- data.matrix(x_cube_variables)


for (i in 1:num_seeds) {

  seed = seeds[i]
  set.seed(seed)
  
lasso = cv.glmnet(xvars, yvar, alpha = 1, nfolds = 5,
                  standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                  #,standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                  # , standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

best_lambda <- lasso$lambda.min
#best_lambda
#plot(lasso)

best_lasso_model <- glmnet(xvars, yvar, alpha = 1, lambda = best_lambda, family = "gaussian",
                           standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                           #  , standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                           #, standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

lasso_coefs = as.matrix(coef(best_lasso_model, s = best_lambda))

lasso_coefs_pull[[as.character(seed)]] = lasso_coefs[-1, , drop = FALSE]

norm_coeffs_scale = lasso_coefs/max(abs(lasso_coefs[-1]))

norm_coeffs[[as.character(seed)]] = norm_coeffs_scale[-1, , drop = FALSE]

y_pred = predict(best_lasso_model, newx = xvars, s = best_lambda)

sst = sum((yvar - mean(yvar))^2)
sse = sum((y_pred - yvar)^2)
r2_scores[i] = 1 - (sse / sst)

}

lasso_coef_mat = as.data.frame(do.call(cbind, lasso_coefs_pull)) 

colnames(lasso_coef_mat) = make.names(colnames(lasso_coef_mat), unique = T)

# Make DF of all LASSO results with mean and std. dev (un-normalized)
lasso_coef_means = lasso_coef_mat %>% 
  mutate(RowNames = rownames(lasso_coef_mat)) %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(contains("s1"))), 
         sd = sd(c_across(contains("s1")))) %>% 
  relocate(mean, .before = s1) %>% 
  relocate(sd, .before = s1) %>% 
  relocate(RowNames, .before = mean)
  

# Bind all normalized LASSO results from 100 iterations
norm_coeffs_matrix = do.call(cbind, norm_coeffs)

mean_coeffs = as.data.frame(norm_coeffs_matrix, row.names = rownames(norm_coeffs_matrix))

colnames(mean_coeffs) = make.names(colnames(mean_coeffs), unique = T)

# Make DF of all LASSO results with mean and std. dev  
mean_coeffs_df = mean_coeffs %>% 
  mutate(RowNames = rownames(mean_coeffs)) %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(contains("s1"))), 
         sd = sd(c_across(contains("s1")))) %>% 
  relocate(mean, .before = s1) %>% 
  relocate(sd, .before = s1) %>% 
  relocate(RowNames, .before = mean)

results_r2 = as.data.frame(r2_scores) 
mean(results_r2$r2_scores)
sd(results_r2$r2_scores)


## With scale, cube, pearson > 0.7 removals
  # TN, Temp, TSS

## With scale, cube, pearson > 0.7 removals, TN removed
  # NPOC, Temp, NO3, TSS, Peaks

# Scatter Plots of Data -----------------------------------------------

## Figure 2b ####
cube_totdr_plot = ggplot(cube_data, aes(y = cube_ERwc, x = cube_TotDr)) +
  geom_point(aes(color = StrOrd), size =4) + theme_bw() +
  scale_color_viridis(name = "Stream Order", discrete = T)+
  stat_cor(data = cube_data, label.x = 2.5, label.y = -1.7, size = 6, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = cube_data, label.x = 2.5, label.y = -1.9, size = 6, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = cube_data, se = FALSE, linetype = "dashed", linewidth = 2)+ 
  xlab(expression("Total Drainage Area (km"^2*")"^(1/3))) +
  ylab(expression("ER"[wc]*" (mg O"[2]*" L"^-1*" d"^-1*")"^(1/3))) +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) 

cube_totdr_plot

ggsave(file.path('./Figures',"Fibure3b_Cube_TotDr_Scatter_Plots.pdf"), plot=cube_totdr_plot, width = 10, height = 8, dpi = 300,device = "pdf") 

## Figure 4 - Cube root scatter plots ####

cube_npoc_plot = ggplot(cube_data, aes(y = cube_ERwc, x = cube_Mean_NPOC)) +
  geom_point(shape = 1, size = 3) + theme_bw() + 
  stat_cor(data = cube_data, label.x = 0.815, label.y = -1.7, size = 4, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = cube_data, label.x = 0.815, label.y = -1.9, size = 4, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = cube_data, se = FALSE, linetype = "dashed") + 
  xlab(expression("DOC (mg L"^-1*")"^(1/3))) +
  ylab("")
  #ylab(expression("ER"[wc]*" (mg O"[2]*" L"^-1*" d"^-1*")"^(1/3)))

cube_tss_plot = ggplot(cube_data, aes(y = cube_ERwc, x = cube_TSS)) +
  geom_point(shape = 1, size = 3) + theme_bw() + 
  stat_cor(data = cube_data, label.x = 0.525, label.y = -1.7, size = 4, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = cube_data, label.x = 0.525, label.y = -1.9, size = 4, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = cube_data, se = FALSE, linetype = "dashed")+ 
  xlab(expression("TSS (mg L"^-1*")"^(1/3))) +
  ylab("")
  #ylab(expression("ER"[wc]*" (mg O"[2]*" L"^-1*" d"^-1*")"^(1/3)))

cube_no3_plot = ggplot(cube_data, aes(y = cube_ERwc, x = cube_NO3_mg_per_L)) +
  geom_point(shape = 1, size = 3) + theme_bw() + 
  stat_cor(data = cube_data, label.x = 0.30, label.y = -1.7, size = 4, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = cube_data, label.x = 0.30, label.y = -1.9, size = 4, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = cube_data, se = FALSE, linetype = "dashed")+ 
  xlab(expression("NO"[3]*" (mg L"^-1*")"^(1/3))) +
  ylab("")
  #ylab(expression("ER"[wc]*" (mg O"[2]*" L"^-1*" d"^-1*")"^(1/3)))

cube_tn_plot = ggplot(cube_data, aes(y = cube_ERwc, x = cube_Mean_TN)) +
  geom_point(shape = 1, size = 3) + theme_bw() + 
  stat_cor(data = cube_data, label.x = 0.35, label.y = -1.7, size = 4, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = cube_data, label.x = 0.35, label.y = -1.9, size = 4, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = cube_data, se = FALSE, linetype = "dashed")+ 
  xlab(expression("TDN (mg L"^-1*")"^(1/3))) +
  ylab("")
 # ylab(expression("ER"[wc]*" (mg O"[2]*" L"^-1*" d"^-1*")"^(1/3)))

cube_temp_plot = ggplot(cube_data, aes(y = cube_ERwc, x = cube_Temp)) +
  geom_point(shape = 1, size = 3) + theme_bw() + 
  stat_cor(data = cube_data, label.x = 2.025, label.y = -1.7, size = 4, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = cube_data, label.x = 2.025, label.y = -1.9, size = 4, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = cube_data, se = FALSE, linetype = "dashed")+ 
  xlab(expression("Temperature (°C)"^(1/3))) +
  ylab("")
  #ylab(expression("ER"[wc]*" (mg O"[2]*" L"^-1*" d"^-1*")"^(1/3)))

cube_lasso_comb = ggarrange(cube_tn_plot, cube_npoc_plot, cube_temp_plot, cube_tss_plot, cube_no3_plot, nrow = 2, ncol = 3, labels = c("(a)", "(b)", "(c)", "(d)", "(e)"), label.x = c(0.85, 0.85, 0.85, 0.85, 0.85), label.y = c(0.95, 0.95, 0.95, 0.95, 0.95))

cube_lasso_comb

# Annotate Figure by adding common "Effect Size" y-axis
cube_lasso_ann = annotate_figure(cube_lasso_comb, left = text_grob(expression("ER"[wc]*" (mg O"[2]*" L"^-1*" d"^-1*")"^(1/3)), rot = 90, size = 15))

ggsave(file.path('./Figures',"Figure4_Cube_Lasso_Combined_Scatter_Plots.png"), plot=cube_lasso_ann, width = 12, height = 8, dpi = 300,device = "png") 


## Untransformed Scatter Plots ####

npoc_plot = ggplot(new_data, aes(y = ERwc, x = Mean_NPOC)) +
  geom_point() + theme_bw() + 
  stat_cor(data = new_data, label.x = 2.25, label.y = 1.5, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = new_data, label.x = 2.25, label.y = 0.9, size = 3, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = new_data, se = FALSE) + 
  ggtitle("NPOC")

tss_plot = ggplot(new_data, aes(y = ERwc, x = TSS)) +
  geom_point() + theme_bw() + 
  stat_cor(data = new_data, label.x = 45, label.y = 1.5, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = new_data, label.x = 45, label.y = 0.9, size = 3, digits = 2, aes(label = paste(..p.label..)))+ #not sig.
  stat_poly_line(data = new_data, se = FALSE)+ 
  ggtitle("TSS") 

tn_plot = ggplot(new_data, aes(y = ERwc, x = Mean_TN)) +
  geom_point() + theme_bw() + 
  stat_cor(data = new_data, label.x = 1.25, label.y = 1.5, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = new_data, label.x = 1.25, label.y = 0.9, size = 3, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = new_data, se = FALSE)+ 
  ggtitle("TN")

transformations_plot = ggplot(new_data, aes(y = ERwc, x = Transformations)) +
  geom_point() + theme_bw() + 
  stat_cor(data = new_data, label.x = 57500, label.y = 1.5, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = new_data, label.x = 57500, label.y = 0.9, size = 3, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = new_data, se = FALSE)+ 
  ggtitle("Transformations")

norm_plot = ggplot(new_data, aes(y = ERwc, x = NormTrans)) +
  geom_point() + theme_bw() + 
  stat_cor(data = new_data, label.x = 7.75, label.y = 1.5, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = new_data, label.x = 7.75, label.y = 0.9, size = 3, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = new_data, se = FALSE)+ 
  ggtitle("Normalized Transformations")

new_data = new_data %>% 
  mutate(size = ifelse(StrOrd <= 3, "small", ifelse(StrOrd == 7, "large", "mid"))) # try assigning stream order "sizes"

totdr_plot = ggplot(new_data, aes(y = ERwc, x = TotDr)) +
  geom_point(aes(color = size)) + theme_bw() + 
  stat_cor(data = new_data, label.x = 10000, label.y = 1.5, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = new_data, label.x = 10000, label.y = 0.9, size = 3, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = new_data, se = FALSE)+ 
  xlab("Total Drainage")

no3_plot = ggplot(new_data, aes(y = ERwc, x = NO3_mg_per_L)) +
  geom_point() + theme_bw() + 
  stat_cor(data = new_data, label.x = 5.5, label.y = 1.5, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = new_data, label.x = 5.5, label.y = 0.9, size = 3, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = new_data, se = FALSE)+ 
  ggtitle("NO3")

temp_plot = ggplot(new_data, aes(y = ERwc, x = Temp)) +
  geom_point() + theme_bw() + 
  stat_cor(data = new_data, label.x = 17, label.y = 1.5, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = new_data, label.x = 17, label.y = 0.9, size = 3, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = new_data, se = FALSE)+ 
  ggtitle("Temperature")

strord_plot = ggplot(new_data, aes(y = ERwc, x = StrOrd)) +
  geom_point() + theme_bw() + 
  stat_cor(data = new_data, label.x = 5.50, label.y = 1.5, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = new_data, label.x = 5.50, label.y = 0.9, size = 3, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = new_data, se = FALSE)+ 
  xlab("Stream Order")

cl_plot = ggplot(new_data, aes(y = ERwc, x = Cl_mg_per_L)) +
  geom_point() + theme_bw() + 
  stat_cor(data = new_data, label.x =5, label.y = 1.5, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = new_data, label.x = 5, label.y = 0.9, size = 3, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = new_data, se = FALSE)+ 
  ggtitle("Cl")

dic_plot = ggplot(new_data, aes(y = ERwc, x = Mean_DIC)) +
  geom_point() + theme_bw() + 
  stat_cor(data = new_data, label.x =27.5, label.y = 1.5, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = new_data, label.x = 27.5, label.y = 0.9, size = 3, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = new_data, se = FALSE)+ 
  ggtitle("DIC")

so4_plot = ggplot(new_data, aes(y = ERwc, x = SO4_mg_per_L)) +
  geom_point() + theme_bw() + 
  stat_cor(data = new_data, label.x =8, label.y = 1.5, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = new_data, label.x = 8, label.y = 0.9, size = 3, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = new_data, se = FALSE)+ 
  ggtitle("SO4")

peaks_plot = ggplot(new_data, aes(y = ERwc, x = Peaks)) +
  geom_point() + theme_bw() + 
  stat_cor(data = new_data, label.x =7600, label.y = 1.5, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = new_data, label.x = 7600, label.y = 0.9, size = 3, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = new_data, se = FALSE)+ 
  ggtitle("Peaks")

combined = ggarrange(npoc_plot, tss_plot, tn_plot, transformations_plot, norm_plot, no3_plot, totdr_plot, temp_plot, strord_plot, dic_plot, cl_plot, so4_plot, peaks_plot, nrow = 4, ncol = 4)

combined

lasso_comb = ggarrange(npoc_plot, temp_plot, no3_plot, tss_plot, nrow = 2, ncol = 2)

lasso_comb

ggsave(file.path('./Figures',"lasso_combined_scatter_plots.png"), plot=lasso_comb, width = 12, height = 12, dpi = 300,device = "png") 
