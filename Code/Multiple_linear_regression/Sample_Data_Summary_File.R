## Read in SPS 2021 data and make into summary file

library(tidyverse)

current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd("../..")
getwd()


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
  summarize(mean_DIC = mean(DIC_mg_L, na.rm = TRUE),
            cv_DIC = (sd(DIC_mg_L, na.rm = T)/mean(DIC_mg_L, na.rm = T))*100,
            count_DIC = n())

## Ions - also missing some replicates here

ion_path = "./Data/Multiple_linear_regression/v2_SPS_Ions.csv"

ion = read_csv(ion_path, skip = 2) %>% 
  slice(-1:-11, -151) %>% 
  select(-c(Field_Name, Material)) 
