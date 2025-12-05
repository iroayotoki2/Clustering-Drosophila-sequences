#Specie and sample size filtering----
dfND2 %>% select(species_name) %>% count(species_name, sort = T)
dfCytb %>% select(species_name) %>% count(species_name, sort = T)

#ND2 has 726 samples with 345 species while Cytb has 541 samples with 115 samples 
#Removing 180 lowest sampled species for ND2 i.e 1 sample per specie would make our sample sizes and number of species more comparable 
species_counts <- dfND2 %>%
  group_by(species_name) %>%
  summarise(n_samples = n()) %>%
  arrange(n_samples)

# Get species to remove
species_to_remove <- species_counts$species_name[1:180]

# Filter them out
dfND2<- dfND2 %>%
  filter(!species_name %in% species_to_remove)

rm(species_counts)

#165 Species and 546 samples for ND2
length(unique(dfND2$species_name))


#116 Species and 541 samples for Cytb
length(unique(dfCytb$species_name))

