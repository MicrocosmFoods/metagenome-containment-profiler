library(tidyverse)
library(ggalluvial)

# sylph results
sylph_results <- read_tsv("results/2025-04-03-results/combined_sylph_profiles.tsv") %>% 
  mutate(accession_name = sample_name) %>% 
  select(-sample_name)

# sample metadata
sample_metadata <- read_tsv("inputs/2025-04-03-selected-carlino-samples-input-table.tsv") %>% 
  unique()
colnames(sample_metadata) <- c("sample_name", "accession_name")

# reference isolate metadata
isolate_metadata <- read_csv("metadata/Probiotic_Species_Strains.csv") %>% 
  mutate(isolate = `Species/Strain`) %>% 
  mutate(isolate_accession = Accession) %>% 
  select(isolate, isolate_accession)

# joined dfs
sylph_results_info <- left_join(sylph_results, sample_metadata, by="accession_name") %>% 
  select(accession_name, sample_name, Genome_file, Contig_name, Adjusted_ANI, Sequence_abundance) %>% 
  mutate(isolate_accession = gsub(".fna", "", Genome_file)) %>% 
  left_join(isolate_metadata) %>% 
  select(accession_name, sample_name, isolate_accession, isolate, Adjusted_ANI, Sequence_abundance) %>%   mutate(sample_code = paste0(accession_name, "_", sample_name)) %>% 
  filter(Sequence_abundance > 25)

sylph_results_info %>% 
  ggplot(aes(x=sample_code, y=isolate, fill=Sequence_abundance)) +
  geom_tile()

sylph_results_info %>% 
  ggplot(aes(axis1 = isolate, axis2 = sample_code, y = Sequence_abundance)) +
  geom_alluvium(aes(fill = isolate), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey80", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Strain", "Food Sample"), expand = c(.05, .05)) +
  labs(title = "Industrial Strains Detected in Fermented Foods",
       y = "Sequence Abundance (%)") +
  theme_minimal()
