# Code to run LASSO on SPS ERwc data 


# Libraries ---------------------------------------------------------------
library(tidyverse)
#library(dplyr)
library(corrplot)
library(ggpubr)
library(ggpmisc)
# library(factoextra)
# library(stringr)
library(glmnet)
# library(magick)
library(readxl)
#library(hal9001)
#library(GGally)


# Working Directory -------------------------------------------------------
current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd("../..")
getwd()

rm(list=ls());graphics.off()

# Functions ---------------------------------------------------------------

cube_root <- function(x) sign(x) * (abs(x))^(1/3)

# Read in Data ------------------------------------------------------------

mean_erwc = read.csv("./Data/Multiple_linear_regression/ERwc_Mean.csv") %>% 
  select(-X)

## are these the only variables I want?
geo = read.csv("./Data/Multiple_linear_regression/v2_RCSFA_Extracted_Geospatial_Data_2023-06-21 (1).csv") %>% 
  select(c(site, streamorde, totdasqkm)) %>% 
  dplyr::rename(Site_ID = site)

data = read.csv("./Data/Multiple_linear_regression/v3_SFA_SpatialStudy_2021_SampleData/v3_SPS_Water_Sample_Data_Summary.csv", skip = 2) %>% 
  filter(grepl("SPS", Sample_Name)) %>% 
  select(c(Sample_Name, Mean_00691_DIC_mg_per_L_as_C))
  

npoc_tn = read.csv("./Data/Multiple_linear_regression/v3_SFA_SpatialStudy_2021_SampleData/v3_SPS_Water_NPOC_TN.csv", skip = 2) %>% 
  filter(grepl("SPS", Sample_Name)) %>% 
  select(c(Sample_Name, X00681_NPOC_mg_per_L_as_C, X00602_TN_mg_per_L_as_N, Methods_Deviation)) %>% 
  rename(NPOC = X00681_NPOC_mg_per_L_as_C) %>% 
  rename(TN = X00602_TN_mg_per_L_as_N) %>% 
  mutate(TN = if_else(grepl("Below", TN), as.numeric(.035), as.numeric(TN))) %>% 
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


## Try with published NO3 data

ions = read.csv("./Data/Multiple_linear_regression/v3_SFA_SpatialStudy_2021_SampleData/v3_SPS_Water_Ions.csv", skip = 2) %>% 
  filter(grepl("SPS", Sample_Name)) %>% 
  mutate(NO3_mg_per_L = ifelse(grepl("Nitrate", X71851_NO3_mg_per_L_as_NO3), 0.035, as.numeric(X71851_NO3_mg_per_L_as_NO3))) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "ION", "Water")) %>% 
  separate(Sample_Name, c("Sample_Name", "Rep"), sep = "-") %>% 
  mutate(Cl_mg_per_L = as.numeric(X00940_Cl_mg_per_L)) %>% 
  mutate(SO4_mg_per_L = as.numeric(X00945_SO4_mg_per_L_as_SO4)) %>% 
  select(c(Sample_Name, NO3_mg_per_L, Cl_mg_per_L, SO4_mg_per_L)) 
  

#Check LOD
tss = read.csv("./Data/Multiple_linear_regression/v3_SFA_SpatialStudy_2021_SampleData/v2_SPS_Water_TSS.csv", skip = 2) %>% 
  filter(grepl("SPS", Sample_Name)) %>% 
  select(c(Sample_Name, X00530_TSS_mg_per_L)) %>% 
  rename(TSS = X00530_TSS_mg_per_L) %>% 
  mutate(TSS = if_else(grepl("Below", TSS), as.numeric(.743), as.numeric(TSS))) %>% 
  separate(Sample_Name, c("Parent","Rep"), sep = "-") %>% 
  mutate(Sample_Name = str_replace(Parent, "_TSS", "_Water")) %>% 
  select(c(Sample_Name, TSS))

# Use means in DP
# dic = read.csv("./Data/Multiple_linear_regression/v3_SFA_SpatialStudy_2021_SampleData/v3_SPS_Water_DIC.csv", skip = 2) %>% 
#   filter(grepl("SPS", Sample_Name)) %>% 
#   select(c(Sample_Name, X00691_DIC_mg_per_L_as_C)) %>% 
#   rename(DIC = X00691_DIC_mg_per_L_as_C)


sample = left_join(data, npoc_tn) %>% 
  left_join(ions) %>% 
  left_join(tss) %>% 
  rename(Mean_DIC = Mean_00691_DIC_mg_per_L_as_C) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "_Water", "")) %>% 
  mutate(Mean_DIC = as.numeric(Mean_DIC))


# add ultrameter water chemistry to this? add manta data to this? 

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

new_names = c(ERwc = "Mean_ERwc", Temp = "Mean_Temp", StrOrd = "streamorde", TotDr = "totdasqkm", Transformations = "Total_Number_of_Transformations", Peaks = "Number_of_Peaks", NormTrans = "Normalized_Transformations"#, TSS = "TSS_mg_per_L", DIC = "DIC_mean", NPOC = "NPOC_mg_per_L_as_C_mean"#, TN = "TN_mg_per_L_as_N_mean"#, Br = #"Br_mg_per_L_mean", Ca = #"Ca_mg_per_L_mean", Cl = #"Cl_mg_per_L_mean", Fl = #"F_mg_per_L_mean", Mg = #"Mg_mg_per_L_mean", NO3 = #"NO3_mg_per_L_as_NO3_mean", K = #"K_mg_per_L_mean", Na = #"Na_mg_per_L_mean", SO4 = #"SO4_mg_per_L_as_SO4_mean"
                )

new_data <- all_data %>% 
  rename(!!!new_names) %>% 
  column_to_rownames("Site_ID") %>% 
  filter(ERwc < 0.5) %>% 
  drop_na()

## Look at histograms

long_data = new_data %>% 
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

ggplot() + 
  geom_histogram(long_data, mapping = aes(x = value)) + 
  facet_wrap(~ variable, scales = "free") +
  theme_minimal()

##Spearman correlation before transformations

spearman <- cor(new_data, method = "spearman", use = "complete.obs")

png(file = paste0("./Figures/", as.character(Sys.Date()),"_Scale_Spearman_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)

corrplot(spearman,type = "upper", method = "number", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25,  title = "Spearman Correlation")

dev.off()

spear.panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
  
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  
  r = (cor(x, y, method = c("spearman")))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  
  if(missing(cex.cor)) {cex.cor <- 0.8/strwidth(txt)}
  text(0.5, 0.5, txt, cex = cex.cor * (1 + abs(r))/2)
  
  # if(missing(cex.cor)) {cex <- 1.2/strwidth(txt)} else {cex = cex.cor}
  # text(0.5, 0.5, txt, cex = cex * sin(sqrt(abs(r))))
  
  test <- cor.test(x,y, method = "spearman")
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " "))
  #text(0.5, 0.5, txt, cex = cex * r)
  text(.5, .8, Signif, cex=cex.cor, col=2)
  
}

panel.smooth <- function(x, y) {
  points(x, y, pch = 19, col = rgb(0.1, 0.2, 0.5, alpha = 0.3))
  abline(lm(y ~ x), col = 'blue', lty = 2)
}

panel.hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1))
  
  h <- hist(x, plot = FALSE, breaks = "FD")
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  
  rect(breaks[-nB], 0, breaks[-1], y, col="grey", border="white", ...)
}

png(file = paste0("./Figures/", as.character(Sys.Date()),"_Pairs_Spearman_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)

pairs(new_data,
      lower.panel = panel.smooth, 
      upper.panel = spear.panel.cor, 
      diag.panel = panel.hist,
      labels = colnames(new_data),
      cex.labels = 0.8) 

dev.off()

##Pearson correlation before transformations
# function for pearson corr matrix

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

png(file = paste0("./Figures/", as.character(Sys.Date()),"_Pairs_Pearson_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)

pairs(new_data,
      lower.panel = panel.smooth, 
      upper.panel = pear.panel.cor, 
      diag.panel = panel.hist,
      labels = colnames(new_data),
      cex.labels = 0.8) 

dev.off()


       
# Transform data ----------------------------------------------------------

#decide how you want to do this, eg, cube or log transform

#should everything be transformed (eg Peaks, stream order, Transformations?)

# log gives a lot of NAs even with + 1
log_data = new_data %>% 
  mutate_if(is.numeric, log10)

long_log_data = log_data %>% 
  rownames_to_column("Sample_Name") %>% 
  pivot_longer(!Sample_Name, names_to = "variable", values_to = "value")

ggplot() + 
  geom_histogram(long_log_data, mapping = aes(x = value)) + 
  facet_wrap(~ variable, scales = "free") +
  theme_minimal()

log_data_plus_one = new_data %>% 
  mutate_if(is.numeric, ~ log10(. + 1))

# cube root allows you to keep sign
# remove NA values from analysis (might change later)
cube_data = new_data %>% 
 mutate(across(where(is.numeric), cube_root)) %>% 
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x)) %>% 
  drop_na()

long_cube_data = cube_data %>% 
  rownames_to_column("Sample_Name") %>% 
  pivot_longer(!Sample_Name, names_to = "variable", values_to = "value")

ggplot() + 
  geom_histogram(long_cube_data, mapping = aes(x = value)) + 
  facet_wrap(~ variable, scales = "free") +
  theme_minimal()


# Check Co-Linearity ------------------------------------------------------
## Pearson without cube transformation
scale_data = as.data.frame(scale(new_data)) %>% 
  rename_with(where(is.numeric), .fn = ~ paste0("scale_", .x)) %>% 
  drop_na()

mean(scale_data$scale_ERwc)
sd(scale_data$scale_ERwc)

scale_pearson <- cor(scale_data, method = "pearson")

png(file = paste0("./Figures/", as.character(Sys.Date()),"_Pairs_Scale_Pearson_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)

pairs(scale_data,
      lower.panel = panel.smooth, 
      upper.panel = pear.panel.cor, 
      diag.panel = panel.hist,
      labels = colnames(new_data),
      cex.labels = 0.8) 

dev.off()


## Scale data before it goes into correlation matrix

scale_cube_data = as.data.frame(scale(cube_data))%>% 
  rename_with(where(is.numeric), .fn = ~ paste0("scale_", .x))

mean(scale_cube_data$scale_cube_ERwc)
sd(scale_cube_data$scale_cube_ERwc)

scale_cube_pearson <- cor(scale_cube_data, method = "pearson")

png(file = paste0("./Figures/", as.character(Sys.Date()),"_Pairs_Scale_Cube_Pearson_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)

pairs(scale_cube_data,
      lower.panel = panel.smooth, 
      upper.panel = pear.panel.cor, 
      diag.panel = panel.hist,
      labels = colnames(scale_cube_data),
      cex.labels = 0.5) 

dev.off()


png(file = paste0("./Figures/", as.character(Sys.Date()),"_Scale_Cube_Pairs_Spearman_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)

pairs(scale_cube_data,
      lower.panel = panel.smooth, 
      upper.panel = spear.panel.cor, 
      diag.panel = panel.hist,
      labels = colnames(scale_cube_data),
      cex.labels = 0.5) 

dev.off()

## Keep NO3, remove TN ####
pearson_df <- as.data.frame(scale_cube_pearson)

row_names_pearson <- rownames(pearson_df)

pearson_df$Variable <- row_names_pearson

spearman_cor = cor(scale_cube_data, method = c("spearman"))
spearman_df = as.data.frame(spearman_cor)
row_names_spearman = rownames(spearman_df)
spearman_df$Variable = row_names_spearman

# Melt the dataframe for plotting
pearson_melted <- reshape2::melt(pearson_df, id.vars = "Variable") %>% 
  filter(value != 1) %>% 
  mutate(value = abs(value)) %>% # do this so it removes in order, and doesn't leave out high negative correlations
  filter(!grepl("ERwc", Variable))# %>% # remove ERwc, don't want it to be removed %>% 
  #filter(!grepl("Mean_TN",Variable) & !grepl("Mean_TN", variable))

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
kept_variables = erwc_melted[!(erwc_melted$Variable %in% removed_variables$Variable_to_Remove), ]


col_to_keep = unique(kept_variables$Variable)
col_to_keep = c(col_to_keep, "scale_cube_ERwc")

scale_cube_variables = scale_cube_data[, col_to_keep, drop = FALSE]


# Start LASSO -------------------------------------------------------------
  ## LASSO with Correlation Matrix Selected Variables 

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
                  standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                  #,standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                  # , standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

best_lambda <- lasso$lambda.min
best_lambda

plot(lasso)

best_lasso_model <- glmnet(xvars, yvar, alpha = 1, lambda = best_lambda, family = "gaussian",
                           standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                           #  , standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                           #, standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

lasso_coefs = coef(best_lasso_model)
coef(best_lasso_model)

yvar_predict <- predict(best_lasso_model, s = best_lambda, newx = xvars)

sst <- sum((yvar - mean(yvar))^2)
sse <- sum((yvar_predict - yvar)^2)

rsq = 1 - sse/sst

rsq #0.428

#check residuals

res = yvar - yvar_predict

plot(res ~ yvar_predict)

## Loop through LASSO to get average over a lot of seeds ####

num_seeds = 100
seeds = sample(1:500, num_seeds)

## Set response variable (Cube_Effect_Size) and scale
yvar <- data.matrix(scale_cube_variables$scale_cube_ERwc)
mean(yvar)
sd(yvar)

norm_coeffs = list()
r2_scores = numeric(num_seeds)

## Set predictor variables and scale
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

norm_coeffs_scale = lasso_coefs/max(abs(lasso_coefs[-1]))

norm_coeffs[[as.character(seed)]] = norm_coeffs_scale[-1, , drop = FALSE]

y_pred = predict(best_lasso_model, newx = xvars, s = best_lambda)

sst = sum((yvar - mean(yvar))^2)
sse = sum((y_pred - yvar)^2)
r2_scores[i] = 1 - (sse / sst)

aic_value = AIC(best_lasso_model$glmnet.fit)

}

norm_coeffs_matrix = do.call(cbind, norm_coeffs)

mean_coeffs = as.data.frame(norm_coeffs_matrix, row.names = rownames(norm_coeffs_matrix))

colnames(mean_coeffs) = make.names(colnames(mean_coeffs), unique = T)
  
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
  # TN, Temp, TSS, Peaks

## With scale, cube, pearson > 0.7 removals, TN removed
  # NPOC, Temp, NO3, TSS, Peaks

npoc_plot = ggplot(new_data, aes(y = ERwc, x = Mean_NPOC)) +
  geom_point() + theme_bw() + 
  stat_cor(data = new_data, label.x = 2.25, label.y = 1.5, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = new_data, label.x = 2.25, label.y = 0.9, size = 3, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = new_data, se = FALSE) + 
  ggtitle("NPOC")

tss_plot = ggplot(new_data, aes(y = ERwc, x = TSS)) +
  geom_point() + theme_bw() + 
  #geom_smooth()
  stat_cor(data = new_data, label.x = 45, label.y = 1.5, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = new_data, label.x = 45, label.y = 0.9, size = 3, digits = 2, aes(label = paste(..p.label..)))+
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

totdr_plot = ggplot(new_data, aes(y = ERwc, x = TotDr)) +
  geom_point() + theme_bw() + 
  stat_cor(data = new_data, label.x = 10000, label.y = 1.5, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = new_data, label.x = 10000, label.y = 0.9, size = 3, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = new_data, se = FALSE)+ 
  ggtitle("Total Drainage")

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
  ggtitle("Stream Order")

cl_plot = ggplot(new_data, aes(y = ERwc, x = Cl_mg_per_L)) +
  geom_point() + theme_bw() + 
  stat_cor(data = new_data, label.x =5, label.y = 1.5, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = new_data, label.x = 5, label.y = 0.9, size = 3, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = new_data, se = FALSE)+ 
  ggtitle("Cl")

dic_plot = ggplot(new_data, aes(y = ERwc, x = mean_DIC)) +
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

## Cube root scatter plots

cube_npoc_plot = ggplot(cube_data, aes(y = cube_ERwc, x = cube_Mean_NPOC)) +
  geom_point() + theme_bw() + 
  stat_cor(data = cube_data, label.x = 1.3, label.y = 0.75, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = cube_data, label.x = 1.3, label.y = 0.55, size = 3, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = cube_data, se = FALSE) + 
  ggtitle("Cube NPOC")

cube_tss_plot = ggplot(cube_data, aes(y = cube_ERwc, x = cube_TSS)) +
  geom_point() + theme_bw() + 
  #geom_smooth()
  stat_cor(data = cube_data, label.x = 3.5, label.y = 0.75, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = cube_data, label.x = 3.5, label.y = 0.55, size = 3, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = cube_data, se = FALSE)+ 
  ggtitle("Cube TSS")

cube_no3_plot = ggplot(cube_data, aes(y = cube_ERwc, x = cube_NO3_mg_per_L)) +
  geom_point() + theme_bw() + 
  stat_cor(data = cube_data, label.x = 1.65, label.y = 0.75, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = cube_data, label.x = 1.65, label.y = 0.55, size = 3, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = cube_data, se = FALSE)+ 
  ggtitle("Cube NO3")

cube_temp_plot = ggplot(cube_data, aes(y = cube_ERwc, x = cube_Temp)) +
  geom_point() + theme_bw() + 
  stat_cor(data = cube_data, label.x = 2.65, label.y = 0.75, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = cube_data, label.x = 2.65, label.y = 0.55, size = 3, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = cube_data, se = FALSE)+ 
  ggtitle("Cube Temperature")

cube_lasso_comb = ggarrange(cube_npoc_plot, cube_temp_plot, cube_no3_plot, cube_tss_plot, nrow = 2, ncol = 2)

cube_lasso_comb

ggsave(file.path('./Figures',"cube_lasso_combined_scatter_plots.png"), plot=cube_lasso_comb, width = 12, height = 12, dpi = 300,device = "png") 

