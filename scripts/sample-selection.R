library(tidyverse)

# all curated MAG metadata
sample_sraruninfo_metadata <- read_tsv("metadata/2025-05-22-sample-sraruninfo-metadata.tsv")

filtered_samples <- sample_sraruninfo_metadata %>% 
  filter(library_layout == "PAIRED") %>% 
  filter(run_total_bases > 200000000)


# random sampling of maximum of 10 samples from each individual food name, prioritizing by run_total_bases for those with higher sequencing depth

top_samples_metadata <- filtered_samples %>%
  group_by(food_name) %>%
  arrange(desc(run_total_bases), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

top_samples_table <- top_samples_metadata %>% 
  select(run_accession, food_name) %>% 
  distinct(run_accession, .keep_all = TRUE) %>%  # make sure each run_accession is actually unique and not repeated
  mutate(accession = run_accession) %>% 
  mutate(sample_name = food_name) %>% 
  select(sample_name, accession)


write_tsv(top_samples_table, "inputs/2025-05-22-selected-samples.tsv")
