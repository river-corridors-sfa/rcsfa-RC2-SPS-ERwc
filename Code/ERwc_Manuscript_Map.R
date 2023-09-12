# ==============================================================================
#
# Make map for ERswc manuscript 
#
# Status: Complete
#
# ==============================================================================
#
# Author: Brieanne Forbes 
# 4 August 2023
#
# ==============================================================================
library(tidyverse) #keep it tidy
library(raster) # work with rasters, NOTE: masks dplyr::select
library(janitor) # clean_names()
library(ggthemes) # theme_map()
library(ggsflabel) # add labels to sf objects
library(ggnewscale) # set multiple color scales
library(ggspatial) # add north arrow and scale bar
library(nhdplusTools) # get watershed boundary/flowlines
library(elevatr) # pull elevation maps
library(sf) # tidy spatial
library(spData)
library(cowplot)
library(rstudioapi)
library(viridis)

rm(list=ls(all=T))

# Setting wd to parent folder
current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
setwd("./..")

# ================================= User inputs ================================

metadata_file <- 'C:/Users/forb086/Downloads/v2_SFA_SpatialStudy_2021_SampleData/SPS_Sample_Field_Metadata.csv'

data_file <- './Data/Multiple_linear_regression/spatial_data.csv'

yrb_shp_dir <- './Data/Maps/YakimaRiverBasin_Boundary'

cluster_shp_dir <- './Data/Maps'

common_crs = 4326

# ============================== read in and merge =============================

metadata <- read_csv(metadata_file) %>%
  dplyr::select(Site_ID, Latitude, Longitude)

data <- read_csv(data_file) %>%
  dplyr::select(IDmapping.Site_ID, DO.slope.mgday)

merge <- data %>%
  left_join(metadata, by = c('IDmapping.Site_ID' = 'Site_ID')) %>%
  rename(Site_ID = IDmapping.Site_ID,
         ER_wc = DO.slope.mgday) %>%
  arrange(ER_wc)

# %>%
#   mutate(ER_wc = if_else(ER_wc > 0, 0, ER_wc))
  
# ============================ read in YRB shp file ============================

YRB_shp <- list.files(yrb_shp_dir, 'shp', full.names = T)

YRB_boundary <- read_sf(YRB_shp) %>%
  st_transform(common_crs)

# ============================ convert to sf object ============================

sites <- st_as_sf(merge, coords = c('Longitude','Latitude'), crs = common_crs)

# ======================== pull NHD data and elevation =========================

YRB_flowlines <- get_nhdplus(AOI = YRB_boundary$geometry, streamorder = 3)

elevation_raw <- get_elev_raster(YRB_boundary$geometry, z = 10)

elevation_crop <- mask(elevation_raw, YRB_boundary)

elevation <- as.data.frame(elevation_crop, xy = T) %>%
  as_tibble() %>%
  rename("long" = x,
         "lat" = y,
         "elevation" = 3) %>% #column index > name (changing resolution changes colname)
  filter(!is.na(elevation))

# ============================ read in cluster shp file ========================

# cluster_shp <- list.files(cluster_shp_dir, 'shp', full.names = T)
# 
# cluster <- read_sf(cluster_shp) %>%
#   st_transform(common_crs)


# ========================= create insert map ======================

data("us_states", package = "spData")
us_states_4326 = st_transform(us_states, crs = 4326)

wa <- us_states_4326 %>% filter(NAME == "Washington")

insert <- ggplot() +
  geom_sf(data = us_states_4326, fill = "white") + 
  geom_sf(data = wa, fill = "black",colour = "black")+
  geom_sf(data = YRB_boundary, colour = "red", fill = 'red') +
  labs(x = "", y = "")+
  theme_map()

# ========================= create map of ER water column (elevation) ======================

ER_wc_map <- ggplot()+
  geom_sf(data = YRB_boundary)+
  geom_raster(data = elevation, aes(long, lat, fill = elevation), show.legend = F, alpha = 0.4)+
  scale_fill_gradient(low = 'white', high = 'black')+
  geom_sf(data = YRB_flowlines, color = "royalblue", alpha = 0.6)+
  new_scale_fill()+
  geom_sf(data = sites, aes(color = ER_wc, size = ER_wc), show.legend = T) +
  scale_fill_viridis(option = 'B', begin = 0.2)+
  scale_color_viridis(option = 'B', begin = 0.2)+
  scale_size(range = c(1.5, 6), trans = 'reverse')+
  theme_map() +
  labs(x = "", y = "", color = "Water Column\nRespiration\n(mg O2 L-1 day-1)") +
  ggspatial::annotation_scale(
    location = "br",
    pad_x = unit(0.5, "in"),
    bar_cols = c("black", "white")) +
  ggspatial::annotation_north_arrow(
    location = "tl", which_north = "true",
    pad_x = unit(1.5, "in"),
    # pad_y = unit(0.5, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("black", "white"),
      line_col = "grey20"))

full <- ggdraw() +
  draw_plot(ER_wc_map) +
  draw_plot(insert, x = 0.4, y = 0.4, width = 0.3, height = 0.3)

ggsave('./Data/Maps/SPS_ER_Water_Column_Map.pdf',
       full,
       width = 8,
       height = 5
)

# ========================= create map of ER water column (cluster) ======================

# ER_wc_map_cluster <- ggplot()+
#   geom_sf(data = YRB_boundary)+
#   geom_sf(data = cluster, aes(fill = as.factor(ClusterNum), color = as.factor(ClusterNum)), show.legend = T)+
#   scale_fill_manual(values = alpha(c('forestgreen', 'steelblue1', 'green4', 'tan', 'palegreen', 'peachpuff'), 0.3))+
#   scale_color_manual(values = alpha(c('forestgreen', 'steelblue1', 'green4', 'tan', 'palegreen', 'peachpuff'), 0.2))+
#   geom_sf(data = YRB_flowlines, color = "royalblue", alpha = 0.6)+
#   new_scale_fill()+
#   new_scale_color()+
#   geom_sf(data = sites, aes(color = ER_wc, size = ER_wc), show.legend = T) +
#   scale_fill_viridis(option = 'B', begin = 0.3)+
#   scale_color_viridis(option = 'B', begin = 0.3)+
#   scale_size(range = c(2, 6.5), trans = 'reverse')+
#   theme_map() + 
#   labs(x = "", y = "", color = "Water Column\nRespiration\n(mg O2 L-1 day-1)") + 
#   ggspatial::annotation_scale(
#     location = "br",
#     pad_x = unit(0.5, "in"), 
#     bar_cols = c("black", "white")) +
#   ggspatial::annotation_north_arrow(
#     location = "tl", which_north = "true",
#     pad_x = unit(1.5, "in"), 
#     # pad_y = unit(0.5, "in"),
#     style = ggspatial::north_arrow_nautical(
#       fill = c("black", "white"),
#       line_col = "grey20"))
# 
# full_cluster <- ggdraw() +
#   draw_plot(ER_wc_map_cluster) +
#   draw_plot(insert, x = 0.4, y = 0.4, width = 0.3, height = 0.3)
# 
# ggsave('./Data/Map/SPS_ER_Water_Column_Map_Cluster.pdf',
#        full_cluster,
#        width = 8,
#        height = 5
# )



