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
  slice(-1:-11) %>% 
  select(c("Sample_Name", "00530_TSS_mg_per_L", "Methods_Deviation")) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "TSS-1", "MEAN"))

## DIC - why only a few replicates? no method deviations

dic_path = "./Data/Multiple_linear_regression/v2_SPS_DIC.csv"

dic = read_csv(dic_path, skip = 2) %>% 
  slice(-1:-11, -107) %>% 
  select(c("Sample_Name", "00691_DIC_mg_per_L_as_C")) 

dic_mean = dic %>% 
  mutate(DIC_mg_L = as.numeric(`00691_DIC_mg_per_L_as_C`))%>% 
  separate(Sample_Name, c("Parent", "Rep"), sep = "-") %>% 
  group_by(Parent) %>% 
  summarize(DIC_mean = mean(DIC_mg_L, na.rm = TRUE),
            DIC_cv = cv(DIC_mg_L),
            DIC_count = n())

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
                        count = ~ sum(!is.na(.)))))

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
                        count = ~ sum(!is.na(.)))))

