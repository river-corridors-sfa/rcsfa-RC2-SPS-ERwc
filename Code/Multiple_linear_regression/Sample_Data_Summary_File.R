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
  mutate(Sample_Name = str_replace(Sample_Name, "TSS", "MEAN"))

## DIC 

dic_path = "./Data/Multiple_linear_regression/v2_SPS_DIC.csv"

dic = read_csv(dic_path, skip = 2) %>% 
  slice(-1:-11, -107) %>% 
  select(c("Sample_Name", "00530_TSS_mg_per_L", "Methods_Deviation")) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "TSS", ""))
