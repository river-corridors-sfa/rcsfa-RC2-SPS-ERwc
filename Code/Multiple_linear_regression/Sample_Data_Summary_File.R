## Read in SPS 2021 data and make into summary file

library(tidyverse)

current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd("../..")
getwd()

## CV function

cv <- function(x) {
  if(length(x) > 1){
    (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))*100
  } else {
    NA
  }
}

## TSS 

tss_path = "./Data/Multiple_linear_regression/SPS_TSS.csv"

tss = read_csv(tss_path, skip = 2) %>% 
  slice(-1:-11, -59) %>% 
  select(c("Sample_Name", "00530_TSS_mg_per_L", "Methods_Deviation")) 

# decide what to do about LOD
tss_mean = tss %>% 
  mutate(Parent = str_replace(Sample_Name, "TSS-1", "MEAN")) %>% 
  mutate(TSS_mg_per_L = as.numeric(str_extract(`00530_TSS_mg_per_L`, "\\d+\\.*\\d*"))) %>% 
  select(c("Parent", "TSS_mg_per_L")) 

## DIC - why only a few replicates? no method deviations

dic_path = "./Data/Multiple_linear_regression/v2_SPS_DIC.csv"

dic = read_csv(dic_path, skip = 2) %>% 
  slice(-1:-11, -107) %>% 
  select(c("Sample_Name", "00691_DIC_mg_per_L_as_C")) 

#decide what to do about LOD/CV
dic_mean = dic %>% 
  mutate(DIC_mg_L = as.numeric(`00691_DIC_mg_per_L_as_C`))%>% 
  separate(Sample_Name, c("Parent", "Rep"), sep = "-") %>% 
  group_by(Parent) %>% 
  summarize(DIC_mean = mean(DIC_mg_L, na.rm = TRUE),
            DIC_cv = cv(DIC_mg_L),
            DIC_count = n()) %>% 
  ungroup() %>% 
  mutate(Parent = str_replace(Parent, "DIC", "MEAN"))

## Ions - also missing some replicates here

ion_path = "./Data/Multiple_linear_regression/v2_SPS_Ions.csv"

ion = read_csv(ion_path, skip = 2) %>% 
  slice(-1:-11, -151) %>% 
  select(-c(Field_Name, Material)) 

ion_mean = ion %>% 
  mutate(across(c(`00000_NH4_mg_per_L_as_NH4`:`00945_SO4_mg_per_L_as_SO4`), ~as.numeric(str_extract(., "\\d+\\.*\\d*")))) %>% 
  separate(Sample_Name, c("Parent", "Rep"), sep = "-") %>% 
  mutate(across(c(`00000_NH4_mg_per_L_as_NH4`:`00945_SO4_mg_per_L_as_SO4`), ~if_else(. == 9999, NA, .))) %>% 
  group_by(Parent) %>% 
  summarize(across(c(`00000_NH4_mg_per_L_as_NH4`:`00945_SO4_mg_per_L_as_SO4`),
                   list(mean = ~ mean(., na.rm = TRUE),
                        cv = ~ cv(.),
                        count = ~ sum(!is.na(.))))) %>%   ungroup() %>% 
  mutate(Parent = str_replace(Parent, "ION", "MEAN"))


## NPOC/TN - CV flags here as well

cn_path = "./Data/Multiple_linear_regression/v2_SPS_NPOC_TN.csv"

cn = read_csv(cn_path, skip = 2) %>% 
  slice(-1:-11, -153) %>% 
  select(-c(Field_Name, Material)) 

cn_mean = cn %>% 
  mutate(across(c(`00681_NPOC_mg_per_L_as_C`:`00602_TN_mg_per_L_as_N`), ~as.numeric(str_extract(., "\\d+\\.*\\d*")))) %>% 
  separate(Sample_Name, c("Parent", "Rep"), sep = "-") %>% 
  mutate(across(c(`00681_NPOC_mg_per_L_as_C`:`00602_TN_mg_per_L_as_N`), ~if_else(. == 9999, NA, .))) %>% 
  group_by(Parent) %>% 
  summarize(across(c(`00681_NPOC_mg_per_L_as_C`:`00602_TN_mg_per_L_as_N`),
                   list(mean = ~ mean(., na.rm = TRUE),
                        cv = ~ cv(.),
                        count = ~ sum(!is.na(.)))))%>%   ungroup() %>% 
  mutate(Parent = str_replace(Parent, "OCN", "MEAN"))

## Summary ####

sum_file = full_join(tss_mean, cn_mean, by = "Parent") %>% 
  full_join(dic_mean, by = "Parent") %>% 
  full_join(ion_mean, by = "Parent") %>% 
  mutate(across(where(is.numeric), ~ ifelse(is.na(.) | is.nan(.), -9999, .)))%>%
  select(-matches("count|cv"))

write.csv(sum_file, "./Data/Multiple_linear_regression/Summary_Not_Cleaned.csv" )

