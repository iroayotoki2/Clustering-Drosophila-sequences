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

#Carry out entrez search to find datasets within the needed paramaters, 2 genes of the same taxonomic group with over 1000 sequences and in similar range, I eventually decided to use ND2 and Adh genes for the Drosophilidae family.
gene1_search <- entrez_search(db = "nucleotide", 
                              term = "drosophilidae [ORGN] AND ND2 [gene]",
                              use_history = T)
gene1_search

gene2_search <- entrez_search(db = "nucleotide", 
                              term = "drosophilidae [ORGN] AND Adh [gene]",
                              use_history = T)
gene2_search

#Read in data for both genes
st_ND2 <- readDNAStringSet("../data/sequenceND2.fasta")
st_Adh <- readDNAStringSet("../data/sequenceAdh.fasta")
names(st_Adh)
