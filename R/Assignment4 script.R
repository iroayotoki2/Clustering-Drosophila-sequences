# BINF6210 Assignment 4: ----
# Author: Iroayo Toki
# Published: December 5, 2025
# Last updated: December 5, 2025
# Model Taxa:Drosophilidae
#Model Genes: ND2, Cytb
#Difference in Clustering patterns of ND2 and Cytb genes for the Drosophilidae Taxonomic family 
library(tidyverse)
library(vegan)
library(viridis)
library(ggplot2)
library(rentrez)
library(Biostrings)
library(cluster)

##_. Entrez search----
#Carrying out entrez search to find data sets within the needed paramaters, 2 genes of the same taxonomic group with over 1000 sequences and in similar range, I eventually decided to use ND2 and Cytb genes for the Drosophilidae family.
gene1_search <- entrez_search(db = "nucleotide", 
                              term = "drosophilidae [ORGN] AND ND2 [gene]",
                              use_history = T)
gene1_search

gene2_search <- entrez_search(db = "nucleotide", 
                              term = "drosophilidae [ORGN] AND Cytb [gene]",
                              use_history = T)
gene2_search

#1. Read in data for both genes----
st_ND2 <- readDNAStringSet("../data/sequenceND2.fasta")
st_Cytb <- readDNAStringSet("../data/sequenceCytb.fasta")

rm(gene1_search, gene2_search)
names(st_ND2)
#Since we will be working with 2 different datasets it would be more efficient to carry out some reproducible filtering and analysis steps by creating functions to carry them out.
#Function for creating the dataframes and adding metadata
fn_dataframe <- function(st) {
  df_new <-  data.frame(title = names(st),sequence = paste(st))
  df_new = df_new %>% 
  mutate(species_name = word(title, 2L, 3L)) %>% 
  mutate(unique_id = word(title, 1L))
  return(df_new)
}
#Create dataframes
dfND2 <- fn_dataframe(st_ND2)
dfCytb <- fn_dataframe(st_Cytb)

#2. Nucleotide Filtering steps----
summary(str_count(dfND2$sequence))  #median= 893
summary(str_count(dfCytb$sequence)) #median= 1026

#Setting missing data to 1% internal N's and length variability to 30% of our median value to keep good amount of our data and remove outliers
missing.data <- 0.01
length.var <- 300

#Filtering function removes N's at start and end, filters for over 1% internal N's and keeps length within a range of 300 from the median
fn_nuc_filter <- function(df) {
  df <- df %>%
    filter(!is.na(sequence)) %>%
    mutate(nucleotides2 = str_remove_all(sequence, "^N+|N+$|-")) %>%
    filter(str_count(nucleotides2, "N") <= (missing.data * str_count(nucleotides2))) %>%
    filter(str_count(nucleotides2) >= median(str_count(nucleotides2)) - length.var &
    str_count(nucleotides2) <= median(str_count(nucleotides2)) + length.var)
  return(df)
}

dfND2 <- fn_nuc_filter(dfND2)
dfCytb <- fn_nuc_filter(dfCytb)

#Graphical exploration of filtered data
hist(str_count(dfND2$sequence))
hist(str_count(dfCytb$sequence)) #Aggregate around the median


#3. Adding Kmer Frequencies (Dinucleotide and Trinucleotide)----

#Creating function
Fn_Kmers <- function(df) {
  df <- as.data.frame(df)
  df$nucleotides2 <- DNAStringSet(df$nucleotides2)
  df <- cbind(df, as.data.frame(dinucleotideFrequency(df$nucleotides2, as.prob = TRUE)))
  df <- cbind(df, as.data.frame(trinucleotideFrequency(df$nucleotides2, as.prob = TRUE)))
  return(df)
}
#Apply
dfND2 <- Fn_Kmers(dfND2)
dfCytb <- Fn_Kmers(dfCytb)


#4. Clustering----
# Creating function to create  distance matrix and clusters with kmer frequencies
Fn_cluster <- function(df) {
  dist_mat= dist(df[,-(1:5)], method = "euclidean")
  Hier_cluster = hclust(dist_mat, method = "single")
  return(list(
    Distance_matrix = dist_mat,
    Cluster= Hier_cluster))
}

#Apply and create list containing distance matrix and cluster
ClustND2 <- Fn_cluster(dfND2)
ClustCytb <- Fn_cluster(dfCytb)


#Visualize clusters labelled with species names 
#Both dendograms show good grouping according to species 
plot(ClustND2$Cluster, labels = dfND2$species_name, cex = 0.3)

plot(ClustCytb$Cluster, labels = dfCytb$species_name, cex = 0.3)

#5. Comparison using internal measures of cluster strength----

#Seperating Cluster object and distance matrices

hc_ND2 <- ClustND2$Cluster
dist_ND2 <- ClustND2$Distance_matrix

hc_Cytb <- ClustCytb$Cluster
dist_Cytb <- ClustCytb$Distance_matrix

#removing Cluster list
rm(ClustCytb, ClustND2)

#Identifying optimal number of clusters(k) for Silhouette index
#This will be done by iterating through values for k to find the highest average silhouette index
#Iteration stops at 19 as √(n/2) can be used as good assumption to calculate max clusters
#Writing function to determine optimal k

Fn_Silhouette_K <- function(hc, dist) {
  sil_values <- numeric(19)
  for (k in 2:19) {
    clusters_k <- cutree(hc, k = k)
    sil_k <- silhouette(clusters_k, dist)
    sil_values[k] <- mean(sil_k[, 3])   # column 3 = silhouette width
    }
  # Show silhouette widths
  print(sil_values)
  
  # Identify best k
  best_k <- which.max(sil_values)
  cat("Optimal k based on silhouette width:", best_k, "\n")

  #Plot best k Silhouette
  clusters_best <- cutree(hc, k = best_k)
  sil_best <- silhouette(clusters_best, dist)
  plot(sil_best, main = paste("Silhouette Plot for Optimal k =", best_k))
  return(list(sil_best= sil_best, Cluster = clusters_best))
}

#Apply 
#Cytb shows a much stronger clustering signal than ND2 due to its  high average silhouette width(0.81) . The samples form tighter, more internally consistent groups with clearer separation between them.

#ND2, with its much lower silhouette values(0.15), shows weaker structure. This infers that Cytb might be better for cluster based classification. 
Sil_ND2 <- Fn_Silhouette_K(hc_ND2, dist_ND2)
Sil_Cytb <- Fn_Silhouette_K(hc_Cytb, dist_Cytb)


#PCA on K-mer Frequencies to visualize 
pca_ND2 <- prcomp(dfND2[ , -(1:5)], scale. = TRUE)
pca_Cytb <- prcomp(dfCytb[ , -(1:5)], scale. = TRUE)
pca_Cytb
#Function for PCA plot of first 2 principal components

plot_pca <- function(pca, clusters, title) {
  pca_df <- data.frame(
    PC1 = pca$x[,1],
    PC2 = pca$x[,2],
    cluster = as.factor(clusters)
  )
  
  ggplot(pca_df, aes(PC1, PC2, color = cluster)) +
    geom_point(size = 2, alpha = 0.8) +
    theme_minimal() +
    ggtitle(title)
}
#Plots
plot_pca(pca_ND2, Sil_ND2$Cluster, "ND2 PCA")
plot_pca(pca_Cytb, Sil_Cytb$Cluster, "CytB PCA")

