
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
setwd("./..")
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

mean_erwc = read.csv("./Data/ERwc_Mean.csv")

## Pull out stream order and total drainage area
geo = read.csv("https://github.com/river-corridors-sfa/Geospatial_variables/raw/refs/heads/main/Archived_versions/v2_RCSFA_Extracted_Geospatial_Data_2023-06-21.csv") %>% 
  select(c(site, streamorde, totdasqkm)) %>% 
  dplyr::rename(Site_ID = site)

## Download data and add to Published_Data folder
##  https://data.ess-dive.lbl.gov/view/doi%3A10.15485%2F1898914, v3_SFA_SpatialStudy_2021_SampleData.zip

# DIC 
data = read.csv("./Data/Published_Data/v3_SFA_SpatialStudy_2021_SampleData/v3_SPS_Water_Sample_Data_Summary.csv", skip = 2) %>% 
  filter(grepl("SPS", Sample_Name)) %>% 
  select(c(Sample_Name, Mean_00691_DIC_mg_per_L_as_C))
  
## DOC/TDN 
npoc_tn = read.csv("./Data/Published_Data/v3_SFA_SpatialStudy_2021_SampleData/v3_SPS_Water_NPOC_TN.csv", skip = 2) %>% 
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


## NO3, SO4, Cl data link
ions = read.csv("./Data/Published_Data/v3_SFA_SpatialStudy_2021_SampleData/v3_SPS_Water_Ions.csv", skip = 2) %>% 
  filter(grepl("SPS", Sample_Name)) %>% 
  mutate(NO3_mg_per_L = if_else(grepl("Nitrate", X71851_NO3_mg_per_L_as_NO3), as.numeric(0.035), as.numeric(X71851_NO3_mg_per_L_as_NO3))) %>% #set samples below standard or LOD to half of LOD (0.035)
  mutate(Sample_Name = str_replace(Sample_Name, "ION", "Water")) %>% 
  separate(Sample_Name, c("Sample_Name", "Rep"), sep = "-") %>% 
  mutate(Cl_mg_per_L = as.numeric(X00940_Cl_mg_per_L)) %>% 
  mutate(SO4_mg_per_L = as.numeric(X00945_SO4_mg_per_L_as_SO4)) %>% 
  select(c(Sample_Name, NO3_mg_per_L, Cl_mg_per_L, SO4_mg_per_L)) 
  
## TSS 
tss = read.csv("./Data/Published_Data/v3_SFA_SpatialStudy_2021_SampleData/v2_SPS_Water_TSS.csv", skip = 2) %>% 
  filter(grepl("SPS", Sample_Name)) %>% 
  select(c(Sample_Name, X00530_TSS_mg_per_L)) %>% 
  rename(TSS = X00530_TSS_mg_per_L) %>% 
  mutate(TSS = if_else(grepl("Below", TSS), as.numeric(.743), as.numeric(TSS))) %>% 
  mutate(TSS = if_else(TSS <= 0.24, 0.12, TSS)) %>% 
  separate(Sample_Name, c("Parent","Rep"), sep = "-") %>% 
  mutate(Sample_Name = str_replace(Parent, "_TSS", "_Water")) %>% 
  select(c(Sample_Name, TSS))

## Organic Matter
om <- read.csv(file.path("./Data/OM_transformation_analysis/SPS_Total_and_Normalized_Transformations_01-03-23.csv")) 

# Join all sample data
sample = left_join(data, npoc_tn) %>% 
  left_join(ions) %>% 
  left_join(tss) %>% 
  rename(Mean_DIC = Mean_00691_DIC_mg_per_L_as_C) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "_Water", "")) %>% 
  mutate(Mean_DIC = as.numeric(Mean_DIC)) %>% 
  left_join(om) 


## Get Site ID and Sample Names to merge with geospatial data: https://data.ess-dive.lbl.gov/view/doi%3A10.15485%2F1892052, v3_SFA_SpatialStudy_2021_SensorData.zip
mapping = read.csv("./Data/Published_Data/v3_SFA_SpatialStudy_2021_SensorData/v2_SPS_Sensor_Field_Metadata.csv") %>% 
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

scale_cube_variables = scale_cube_data %>% 
  select(-c(scale_cube_StrOrd))#[, col_to_keep, drop = FALSE]

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

xvars <- data.matrix(x_cube_variables)


for (i in 1:num_seeds) {

  seed = seeds[i]
  set.seed(seed)
  
# cross validation
lasso = cv.glmnet(xvars, yvar, alpha = 1, nfolds = 5,
                  standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                  #,standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                  # , standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

best_lambda <- lasso$lambda.min

#best_lambda = lasso$lambda.1se
# best_lambda
# plot(lasso)

best_lasso_model <- glmnet(xvars, yvar, alpha = 1, lambda = best_lambda, family = "gaussian",
                           standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                           #  , standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                           #, standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

coef(best_lasso_model)

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
  select(where(~!any(is.nan(.)))) %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(contains("s1"))), 
         sd = sd(c_across(contains("s1")))) %>% 
  relocate(mean, .before = s1) %>% 
  relocate(sd, .before = s1) %>% 
  relocate(RowNames, .before = mean)

results_r2 = as.data.frame(r2_scores) %>% 
  filter(r2_scores != 0)
mean(results_r2$r2_scores)
sd(results_r2$r2_scores)

 ## With scale, cube, pearson > 0.7 removals
  # TN, Temp, TSS

## With scale, cube, pearson > 0.7 removals, TN removed
  # NPOC, Temp, NO3, TSS, Peaks

# Scatter Plots of Data -----------------------------------------------

totdr_plot = ggplot(new_data, aes(y = ERwc, x = TotDr)) +
  geom_point(aes(color = as.factor(StrOrd), size =4)) + theme_bw() +
  scale_color_viridis(name = "Stream Order", discrete = T)+
  stat_cor(data = new_data, label.x = 2.5, label.y = -1.7, size = 6, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = new_data, label.x = 2.5, label.y = -2.1, size = 6, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = new_data, se = FALSE, linetype = "dashed", linewidth = 2)+ 
  xlab(expression("Total Drainage Area (km"^2*")"^(1/3))) +
  ylab(expression("ER"[wc]*" (mg O"[2]*" L"^-1*" d"^-1*")"^(1/3))) +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) 

## Figure 2b ####
cube_totdr_plot = ggplot(cube_data, aes(y = cube_ERwc, x = cube_TotDr)) +
  geom_point(aes(color = as.factor(cube_StrOrd^3)), size =4) + theme_bw() +
  scale_color_viridis(name = "Stream Order", discrete = T)+
  stat_cor(data = cube_data, label.x = 2.5, label.y = -1.7, size = 6, digits = 2, aes(label = paste(..r.label..)))+
  #stat_cor(data = cube_data, label.x = 2.5, label.y = -1.9, size = 6, digits = 2, aes(label = paste(..p.label..)))+
  #stat_poly_line(data = cube_data, se = FALSE, linetype = "dashed", linewidth = 2)+ 
  xlab(expression("Total Drainage Area (km"^2*")"^(1/3))) +
  ylab(expression("ER"[wc]*" (g O"[2]*" m"^-3*" d"^-1*")"^(1/3))) +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) 

cube_totdr_plot

ggsave(file.path('./Figures',"Figure3b_Cube_TotDr_Scatter_Plots.pdf"), plot=cube_totdr_plot, width = 10, height = 8, dpi = 300,device = "pdf") 


## Try boxplots of stream order?

ggplot(cube_data, aes(y = cube_ERwc, x = (cube_StrOrd^3), fill = (cube_StrOrd^3), group = cube_StrOrd)) +
  geom_boxplot() +
  theme_bw() +
  #geom_point(aes(color = as.factor(cube_StrOrd^3)), size =4) + theme_bw() +
  scale_fill_viridis(name = "Stream Order")+
  stat_cor(data = cube_data, label.x = 2.5, label.y = -1.7, size = 6, digits = 2, aes(label = paste(..r.label..)))+
  #stat_cor(data = cube_data, label.x = 2.5, label.y = -1.9, size = 6, digits = 2, aes(label = paste(..p.label..)))+
  #stat_poly_line(data = cube_data, se = FALSE, linetype = "dashed", linewidth = 2)+ 
  xlab(expression("Total Drainage Area (km"^2*")"^(1/3))) +
  ylab(expression("ER"[wc]*" (g O"[2]*" m"^-3*" d"^-1*")"^(1/3))) +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) 


## Figure 4 - Cube root scatter plots ####

cube_npoc_plot = ggplot(cube_data, aes(y = cube_ERwc, x = cube_Mean_NPOC)) +
  geom_point(shape = 1, size = 3) + theme_bw() + 
  stat_cor(data = cube_data, label.x = 0.815, label.y = -1.7, size = 5, digits = 2, aes(label = paste(..r.label..)))+
  #stat_cor(data = cube_data, label.x = 0.815, label.y = -1.9, size = 4, digits = 2, aes(label = paste(..p.label..)))+
  #stat_poly_line(data = cube_data, se = FALSE, linetype = "dashed") + 
  xlab(expression("DOC (mg L"^-1*")"^(1/3))) +
  ylab("")
  #ylab(expression("ER"[wc]*" (mg O"[2]*" L"^-1*" d"^-1*")"^(1/3)))

cube_tss_plot = ggplot(cube_data, aes(y = cube_ERwc, x = cube_TSS)) +
  geom_point(shape = 1, size = 3) + theme_bw() + 
  stat_cor(data = cube_data, label.x = 0.525, label.y = -1.7, size = 5, digits = 2, aes(label = paste(..r.label..)))+
  #stat_cor(data = cube_data, label.x = 0.525, label.y = -1.9, size = 4, digits = 2, aes(label = paste(..p.label..)))+
  #stat_poly_line(data = cube_data, se = FALSE, linetype = "dashed")+ 
  xlab(expression("TSS (mg L"^-1*")"^(1/3))) +
  ylab("")
  #ylab(expression("ER"[wc]*" (mg O"[2]*" L"^-1*" d"^-1*")"^(1/3)))

cube_no3_plot = ggplot(cube_data, aes(y = cube_ERwc, x = cube_NO3_mg_per_L)) +
  geom_point(shape = 1, size = 3) + theme_bw() + 
  stat_cor(data = cube_data, label.x = 0.30, label.y = -1.7, size = 5, digits = 2, aes(label = paste(..r.label..)))+
  #stat_cor(data = cube_data, label.x = 0.30, label.y = -1.9, size = 4, digits = 2, aes(label = paste(..p.label..)))+
  #stat_poly_line(data = cube_data, se = FALSE, linetype = "dashed")+ 
  xlab(expression("NO"[3]*" (mg L"^-1*")"^(1/3))) +
  ylab("")
  #ylab(expression("ER"[wc]*" (mg O"[2]*" L"^-1*" d"^-1*")"^(1/3)))

cube_tn_plot = ggplot(cube_data, aes(y = cube_ERwc, x = cube_Mean_TN)) +
  geom_point(shape = 1, size = 3) + theme_bw() + 
  stat_cor(data = cube_data, label.x = 0.35, label.y = -1.7, size = 5, digits = 2, aes(label = paste(..r.label..)))+
  #stat_cor(data = cube_data, label.x = 0.35, label.y = -1.9, size = 4, digits = 2, aes(label = paste(..p.label..)))+
 # stat_poly_line(data = cube_data, se = FALSE, linetype = "dashed")+ 
  xlab(expression("TDN (mg L"^-1*")"^(1/3))) +
  ylab("")
 # ylab(expression("ER"[wc]*" (mg O"[2]*" L"^-1*" d"^-1*")"^(1/3)))

cube_temp_plot = ggplot(cube_data, aes(y = cube_ERwc, x = cube_Temp)) +
  geom_point(shape = 1, size = 3) + theme_bw() + 
  stat_cor(data = cube_data, label.x = 2.025, label.y = -1.7, size = 5, digits = 2, aes(label = paste(..r.label..)))+
  #stat_cor(data = cube_data, label.x = 2.025, label.y = -1.9, size = 4, digits = 2, aes(label = paste(..p.label..)))+
  #stat_poly_line(data = cube_data, se = FALSE, linetype = "dashed")+ 
  xlab(expression("Temperature (°C)"^(1/3))) +
  ylab("")
  #ylab(expression("ER"[wc]*" (mg O"[2]*" L"^-1*" d"^-1*")"^(1/3)))

# cube_lasso_comb = ggarrange(cube_tn_plot, cube_npoc_plot, cube_temp_plot, cube_tss_plot, cube_no3_plot, nrow = 2, ncol = 3, labels = c("(a)", "(b)", "(c)", "(d)", "(e)"), label.x = c(0.85, 0.85, 0.85, 0.85, 0.85), label.y = c(0.95, 0.95, 0.95, 0.95, 0.95))

cube_lasso_comb = ggarrange(cube_tn_plot, cube_temp_plot, cube_npoc_plot, cube_tss_plot, nrow = 2, ncol = 2, labels = c("(a)", "(b)", "(c)", "(d)"), label.x = c(0.85, 0.85, 0.85, 0.85), label.y = c(0.95, 0.95, 0.95, 0.95))

cube_lasso_comb

# Annotate Figure by adding common "Effect Size" y-axis
cube_lasso_ann = annotate_figure(cube_lasso_comb, left = text_grob(expression("ER"[wc]*" (g O"[2]*" m"^-3*" d"^-1*")"^(1/3)), rot = 90, size = 15))

ggsave(file.path('./Figures',"Figure4_Cube_Lasso_Combined_Scatter_Plots.png"), plot=cube_lasso_ann, width = 9, height = 6.5, dpi = 300,device = "png") 

