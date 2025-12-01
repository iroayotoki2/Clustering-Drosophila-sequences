# BINF6210 Assignment 4: ----
# Author: Iroayo Toki
# Published: December 5, 2025
# Last updated: November 29, 2025
# Model Organism:

library(tidyverse)
library(vegan)
library(viridis)
library(ggplot2)
library(rentrez)
library(Biostrings)

#Carrying out entrez search to find data sets within the needed paramaters, 2 genes of the same taxonomic group with over 1000 sequences and in similar range, I eventually decided to use ND2 and cytb genes for the Drosophilidae family.
gene1_search <- entrez_search(db = "nucleotide", 
                              term = "drosophilidae [ORGN] AND ND2 [gene]",
                              use_history = T)
gene1_search

gene2_search <- entrez_search(db = "nucleotide", 
                              term = "drosophilidae [ORGN] AND cytb [gene]",
                              use_history = T)
gene2_search

#1. Read in data for both genes----
st_ND2 <- readDNAStringSet("../data/sequenceND2.fasta")
st_Cytb <- readDNAStringSet("../data/sequenceCYTB.fasta")

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

#2. Filtering steps----
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
hist(str_count(dfCytb$sequence))  #Aggregate around the median
