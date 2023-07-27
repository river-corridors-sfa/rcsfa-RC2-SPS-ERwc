

# Load Packages:
library(sf)


# Define Input and output Paths:
input_path = "//pnl/projects/ColumbiaGIS/Tables/zonal/_R_/ClusterAnalysis_Files_Scripts/R_Inputs/"
output_path = "//pnl/projects/ColumbiaGIS/Tables/zonal/_R_/ClusterAnalysis_Files_Scripts/R_Outputs/"



# Define Data Used in Analysis -------------------------------------------
n_clusters = 6
Data = Data_norm_8
Centers = centers_8



# Run Cluster Analysis ---------------------------------------------------
x = Data[,2:ncol(Data)]
set.seed(1234)
kpres = kmeans(x, Centers, nstart = 10, iter.max = 150)



# Extract Results --------------------------------------------------------

# Sum of Squares:
R = matrix(NA, 1, 4 + n_clusters)

R[1,1] = kpres$totss
R[1,2] = kpres$tot.withinss
R[1,3] = kpres$betweenss
R[1,4] = (kpres$betweenss/kpres$totss)*100
R[1,5:ncol(R)] = kpres$withinss

R = as.data.frame(R)
colnames(R) = c("Total_SS", "Total_Within_SS", "Between_SS", "Goodness_(BSS/TSS*100)", paste0("WSS_Cluster", seq(1:n_clusters)))


# Clusters:
catchment_ID = Data[,1]
clusters = kpres$cluster
P = as.data.frame(cbind(catchment_ID, clusters))


# Join clusters to shape file:
shapefile = st_read(paste0(input_path, "CB_catchments_clean_UTM11.shp"))
cluster_data = P
shape_cluster = left_join(shapefile, cluster_data, by = c("GRIDCODE" = "catchment_ID"))





# Save Results in Outputs Directory --------------------------------------

# Create new directory and export results:
dir.create(paste0(output_path, "Kmeans_Results_", n_clusters, "Clusters_", Sys.Date()))
new_output_path = paste0(output_path, "Kmeans_Results_", n_clusters, "Clusters_", Sys.Date(), "/")


### NOTE: if running multiple analyses in a day, results will overwrite!!! 
### If this is not the intention, make sure to save previous results elsewhere

write.csv(R, paste0(new_output_path, "Kmeans_SSResults_", Sys.Date(), ".csv"), row.names = FALSE)
write.csv(P, paste0(new_output_path, "Kmeans_Clusters_", Sys.Date(), ".csv"), row.names = FALSE)
st_write(shape_cluster, paste0(new_output_path, "CB_catchments_clean_joined_with_Kmeans_Results_", Sys.Date(), ".shp"))



























