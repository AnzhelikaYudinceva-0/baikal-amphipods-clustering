library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(vegan)
library(cluster)
library(dplyr)
library(dbscan)
library(mclust)
library(fpc)
library(ggsci)
library(factoextra)

# Data preparation
#uncomment and modify the line below:
# setwd("/path/to/your/data/folder")

mydata <- read.csv("depth_data.csv", sep = ";")

mydata_ed <- (mydata[,-6])

data_for_analysis <- pivot_longer(mydata_ed, cols = c(minimum_frequent_depth, maximum_frequent_depth, minimum_depth, maximum_depth),
                                  names_to = "Depth_type", values_to = "Depth")

# Plotting a boxplot

Boxplot_Depth <- ggplot(data_for_analysis, aes(x = Depth_type, y = Depth, fill = Depth_type)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "depth_distribution",
       x = "Depth_type", y = "Depth",
       fill = "Depth_type") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("My_Boxplot_Depth.pdf", plot = Boxplot_Depth, width = 6, height = 6, device = "pdf")

# Plotting a scatterplot

Scatterplot_Depth <- ggplot(mydata_ed, aes(x = minimum_depth, y = maximum_depth)) +
  geom_point(aes(color = maximum_depth - minimum_depth), alpha = 0.7, size = 3) +
  scale_color_gradient(low = "blue", high = "red") + 
  labs(title = "Minimum depth vs Maximum_depth",
       x = "minimum_depth", y = "maximum_depth",
       color = "Difference_in_depths") + 
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


ggsave("My_Scatterplot_Depth.pdf", plot = Scatterplot_Depth, width = 6, height = 6, device = "pdf")


# Clustering

# data scaling

scaled_data <- scale(mydata[, c(2, 3, 4, 5)])

pca_result <- prcomp(scaled_data, center = FALSE, scale = TRUE)

pca_data <- pca_result$x[, 1:2]

# Clustering using DBSCAN and GMM methods

# DBSCAN без PCA
dbscan_model <- dbscan::dbscan(scaled_data, eps = 0.1, minPts = 4)
dbscan_clusters <- dbscan_model$cluster

# DBSCAN with PCA
dbscan_model_pca <- dbscan::dbscan(pca_data, eps = 0.1, minPts = 4)
dbscan_clusters_pca <- dbscan_model_pca$cluster

# Helper function to calculate silhouette index (excluding noise points)

calculate_silhouette_dbscan <- function(data, clusters) {
  non_noise <- clusters != 0
  data_non_noise <- data[non_noise, ]
  clusters_non_noise <- clusters[non_noise]
  
  if (length(unique(clusters_non_noise)) > 1) {
    sil <- silhouette(clusters_non_noise, dist(data_non_noise))
    return(mean(sil[, "sil_width"]))
  }
  return(NA)
}

# Helper function to create summary table with min/max depth values per cluster

create_summary_table <- function(depth_data, clusters) {
  do.call(rbind, lapply(unique(clusters), function(cluster_num) {
    cluster_data <- depth_data[clusters == cluster_num, ]
    if(nrow(cluster_data) > 0) {
      min_vals <- apply(cluster_data, 2, min)
      max_vals <- apply(cluster_data, 2, max)
    } else {
      min_vals <- rep(NA, ncol(depth_data))
      max_vals <- rep(NA, ncol(depth_data))
    }
    data.frame(
      cluster = cluster_num,
      depth_col = colnames(depth_data),
      min_value = round(min_vals, 2),
      max_value = round(max_vals, 2)
    )
  }))
}

# Calculate number of clusters (excluding noise cluster 0)

dbscan_nclusters <- length(unique(dbscan_clusters[dbscan_clusters != 0]))
dbscan_nclusters_pca <- length(unique(dbscan_clusters_pca[dbscan_clusters_pca != 0]))

# Count noise points

dbscan_noise_points <- sum(dbscan_clusters == 0)
dbscan_noise_points_pca <- sum(dbscan_clusters_pca == 0)

# Calculate silhouette indices

dbscan_silhouette_avg <- calculate_silhouette_dbscan(scaled_data, dbscan_clusters)
dbscan_silhouette_avg_pca <- calculate_silhouette_dbscan(pca_data, dbscan_clusters_pca)

# Visualize DBSCAN without PCA

DBSCAN_without_PCA <- ggplot(as.data.frame(scaled_data), 
       aes(x = minimum_frequent_depth, y = maximum_frequent_depth, color = factor(dbscan_clusters))) +
  geom_point() +
  labs(title = "DBSCAN without PCA", 
       x = "Minimum frequent depth (scaled)", 
       y = "Maximum frequent depth (scaled)") +
  theme_bw()

ggsave("DBSCAN_without_PCA.pdf", plot = DBSCAN_without_PCA, width = 6, height = 6, device = "pdf")

# Visualize DBSCAN with PCA

DBSCAN_with_PCA <- ggplot(as.data.frame(pca_data), 
       aes(x = PC1, y = PC2, color = factor(dbscan_clusters_pca))) +
  geom_point() +
  labs(title = "DBSCAN Clustering with PCA", x = "PC1", y = "PC2") +
  theme_bw()

ggsave("DBSCAN_with_PCA.pdf", plot = DBSCAN_with_PCA, width = 6, height = 6, device = "pdf")

# Print comparison results

cat("Number_of_clusters:\n")
cat("Without_PCA:", dbscan_nclusters, "\n")
cat("With_PCA:", dbscan_nclusters_pca, "\n\n")

cat("Noise_points:\n")
cat("Without_PCA:", dbscan_noise_points, "\n")
cat("With_PCA:", dbscan_noise_points_pca, "\n\n")

# Print silhouette results

if (!is.na(dbscan_silhouette_avg)) {
  cat("Average silhouette index for DBSCAN without PCA:", dbscan_silhouette_avg, "\n")
} else {
  cat("Cannot calculate silhouette index (only one cluster)\n")
}

if (!is.na(dbscan_silhouette_avg_pca)) {
  cat("Average silhouette index for DBSCAN with PCA:", dbscan_silhouette_avg_pca, "\n")
} else {
  cat("Cannot calculate silhouette index (only one cluster)\n")
}

# GMM without PCA
gmm_model <- mclust::Mclust(scaled_data)
gmm_clusters <- gmm_model$classification
gmm_nclusters <- gmm_model$G

# GMM with PCA
gmm_model_pca <- mclust::Mclust(pca_data)
gmm_clusters_pca <- gmm_model_pca$classification
gmm_nclusters_pca <- gmm_model_pca$G

# Helper function to calculate silhouette index

calculate_silhouette_single <- function(scaled_data, clusters) {
  if (length(unique(clusters)) > 1) {
    sil <- silhouette(clusters, dist(scaled_data))
    return(mean(sil[, "sil_width"]))
  }
  return(NA)
}

# Calculate silhouette indices
gmm_silhouette_avg <- calculate_silhouette_single(scaled_data, gmm_clusters)
gmm_silhouette_avg_pca <- calculate_silhouette_single(pca_data, gmm_clusters_pca)

# Visualize GMM without PCA
GMM_without_PCA <- ggplot(as.data.frame(scaled_data), 
       aes(x = minimum_frequent_depth, y = maximum_frequent_depth, color = factor(gmm_clusters))) +
  geom_point() +
  labs(title = "GMM Clustering without PCA", 
       x = "Minimum frequent depth (scaled)", 
       y = "Maximum frequent depth (scaled)") +
  theme_bw()

ggsave("GMM_without_PCA.pdf", plot = GMM_without_PCA, width = 6, height = 6, device = "pdf")

# Visualize GMM with PCA

GMM_with_PCA <- ggplot(as.data.frame(pca_data), 
       aes(x = PC1, y = PC2, color = factor(gmm_clusters_pca))) +
  geom_point() +
  labs(title = "GMM Clustering with PCA", x = "PC1", y = "PC2") +
  theme_bw()

ggsave("GMM_with_PCA.pdf", plot = GMM_with_PCA, width = 6, height = 6, device = "pdf")

# Print comparison results

cat("Number of clusters:\n")
cat("Without PCA:", gmm_nclusters, "\n")
cat("With PCA:", gmm_nclusters_pca, "\n\n")

# Print silhouette results

if (!is.na(gmm_silhouette_avg)) {
  cat("Average silhouette index for GMM without PCA:", gmm_silhouette_avg, "\n")
} else {
  cat("Cannot calculate silhouette index for GMM without PCA (only one cluster)\n")
}

if (!is.na(gmm_silhouette_avg_pca)) {
  cat("Average silhouette index for GMM with PCA:", gmm_silhouette_avg_pca, "\n")
} else {
  cat("Cannot calculate silhouette index for GMM with PCA (only one cluster)\n")
}

# k-means Clustering

# Function for Calculation of WSS for Varying Numbers of Clusters
calculate_wss <- function(scaled_data, max_k = 2) {
  wss <- sapply(2:max_k, function(k) {
    tryCatch({
      km <- kmeans(scaled_data, centers = k)
      km$tot.withinss
    }, error = function(e) {
      print(paste0("Error for k=", k, ": ", e))
      NA
    })
  })
  return(wss)
}

# Determining the Optimal Cluster Number

# Calculation of WSS for Varying Cluster Numbers
wss_values <- calculate_wss(scaled_data, max_k = 10)

wss_values

# Elbow Method Visualization

# Creating a Data Frame for ggplot2
elbow_data <- data.frame(
  clusters = 2:10,
  wss = wss_values
)

# Visualization Using ggplot2
Elbow_Method <- ggplot(elbow_data, aes(x = clusters, y = wss)) +
  geom_line() +
  geom_point(size = 3) +
  labs(title = "Elbow Method",
       x = "Number of clusters",
       y = "Within-Cluster Sum of Squares (WSS)") +
  ylim(0, max(wss_values) * 1.1) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Saving Plots
ggsave("Elbow_Method.pdf", plot = Elbow_Method, width = 6, height = 6, device = "pdf")

# Silhouette Score

# Function for Silhouette Score Calculation for Different Numbers of Clusters
calculate_silhouette_range <- function(scaled_data, max_k = 10) {
  silhouette_scores <- sapply(2:max_k, function(k) {
    km <- kmeans(scaled_data, centers = k)
    sil <- silhouette(km$cluster, dist(scaled_data))
    mean(sil[, "sil_width"])
  })
  return(silhouette_scores)
}

# Calculating Silhouette Scores for Different Numbers of Clusters
silhouette_scores <- calculate_silhouette_range(scaled_data, max_k = 10)

# Creating a Data Frame for ggplot2
df_silhouette <- data.frame(
  clusters = 2:10,
  score = silhouette_scores
)

# Visualization Using ggplot2
Plot_silhouette <- ggplot(df_silhouette, aes(x = clusters, y = score)) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 2:10) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Silhouette Score",
    x = "Number of clusters",
    y = "Average silhouette width"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "grey70"),
    plot.title = element_text(hjust = 0.5)
  )

# Saving Plots
ggsave("Plot_silhouette.pdf", plot = Plot_silhouette, width = 6, height = 6)

# Calinski–Harabasz Index

# Function for Calinski–Harabasz Index Calculation for Different Numbers of Clusters
calculate_ch_index <- function(scaled_data, max_k = 10) {
  ch_indices <- sapply(2:max_k, function(k) {
    km <- kmeans(scaled_data, centers = k)
    cluster.stats(dist(scaled_data), km$cluster)$ch
  })
  return(ch_indices)
}

# Calinski–Harabasz Index Calculation for Different Numbers of Clusters
ch_indices <- calculate_ch_index(scaled_data, max_k = 10)

# Creating a Data Frame for ggplot2
ch_data <- data.frame(
  clusters = 2:10,
  ch_index = ch_indices
)

# Visualization Using ggplot2
Calinski_Harabasz_plot <- ggplot(ch_data, aes(x = clusters, y = ch_index)) +
  geom_line() +
  geom_point(size = 3, pch = 19) +
  labs(title = "Calinski-Harabasz Index",
       x = "Number of clusters",
       y = "Index Value") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Saving Plots
ggsave("Calinski_Harabasz.pdf", plot = Calinski_Harabasz_plot, width = 6, height = 6, device = "pdf")

# The Selected Optimal Number of Clusters
optimal_k_kmeans <- as.integer(readline(prompt="Optimal Number of Clusters for k-means: "))

# k-means Without PCA
kmeans_without_pca <- kmeans(scaled_data, centers = optimal_k_kmeans, nstart = 25)

# Assessment of k-means Clustering Quality Without PCA
wss_before_pca <- sum(kmeans_without_pca$withinss)

silhouette_output_before_pca <- silhouette(kmeans_without_pca$cluster, dist(scaled_data))
silhouette_score_before_pca <- mean(silhouette_output_before_pca[, "sil_width"])

ch_index_before_pca <- cluster.stats(dist(scaled_data), kmeans_without_pca$cluster)$ch

# k-means with PCA 
kmeans_with_pca <- kmeans(pca_data, centers = optimal_k_kmeans, nstart = 20)

# Assessment of k-means Clustering Quality with PCA
wss_after_pca <- sum(kmeans_with_pca$withinss)

silhouette_output_after_pca <- silhouette(kmeans_with_pca$cluster, dist(pca_data))
silhouette_score_after_pca <- mean(silhouette_output_after_pca[, "sil_width"])

ch_index_after_pca <- cluster.stats(dist(pca_data), kmeans_with_pca$cluster)$ch

# Comparison of Results
cat("Within-Cluster Sum of Squares (WSS):\n")
cat("Before PCA:", wss_before_pca, "\n")
cat("After PCA:", wss_after_pca, "\n\n")

cat("Silhouette Score:\n")
cat("Before PCA:", silhouette_score_before_pca, "\n")
cat("After PCA:", silhouette_score_after_pca, "\n\n")

cat("Calinski Harabasz Index:\n")
cat("Before PCA:", ch_index_before_pca, "\n")
cat("After PCA:", ch_index_after_pca, "\n\n")

kmeans_Plot <- fviz_cluster(object = kmeans_without_pca, data = scaled_data, 
                            geom = "point", 
                            show.clust.cent = FALSE, 
                            repel = TRUE, 
                            labelsize = 0, 
                            pointsize = 1,
                            alpha = 0.7) +
  scale_fill_bmj() +
  scale_color_bmj() +
  labs(title = "k-means Clusters ",
       subtitle = paste0("Number of Clusters: ", length(unique(kmeans_without_pca$cluster)))) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        legend.position = "bottom")

ggsave("kmeans_Plot.pdf", plot = kmeans_Plot, width = 8, height = 8, dpi = 300)

# Creating a Data Frame with Results
results_df <- data.frame(
  species = mydata$species,
  minimum_frequent_depth = mydata$minimum_frequent_depth,
  maximum_frequent_depth = mydata$maximum_frequent_depth,
  minimum_depth = mydata$minimum_depth,
  maximum_depth = mydata$maximum_depth,
  cluster = kmeans_without_pca$cluster
)

# Exporting Results to a CSV File
write.csv(results_df, "clustering_results.csv", row.names = FALSE)

summary_table <- do.call(rbind, lapply(unique(results_df$cluster), function(cluster_num) {
  cluster_data <- subset(results_df, cluster == cluster_num)
  min_vals <- apply(cluster_data[, 2:5], 2, min)
  max_vals <- apply(cluster_data[, 2:5], 2, max)
  data.frame(
    cluster = cluster_num,
    depth_col = names(cluster_data[, 2:5]),
    min_value = round(min_vals, 2),
    max_value = round(max_vals, 2)
  )
}))

summary_table
