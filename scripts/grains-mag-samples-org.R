library(tidyverse)

# Find fermented grains samples to profile metagenomes/genome abundance

# list of fermented grains foods
grain_foods <- read_tsv("metadata/african-foods-list.txt", col_names = c("food")) %>% 
  mutate(fermented_food = tolower(food)) %>% 
  select(-food)

food_list <- grain_foods %>% 
  pull(fermented_food)

# mag/sample metadata
# original metadata
mag_sample_metadata <- read_csv("metadata/Food_MAGs_curated_metadata_250421_corrected_merged_final_v2_corrected.csv")

# metadata filtered with african foods hits
grains_mag_sample_metadata <- mag_sample_metadata %>% 
  filter(fermented_food %in% food_list)

grains_mag_sample_metadata %>% 
  group_by(fermented_food, run_accession) %>% 
  count() %>% 
  print(n=100) %>% 
  arrange(desc(n))

grains_mag_sample_metadata %>% 
  select(fermented_food, run_accession) %>% 
  mutate(run_accession = str_remove_all(run_accession,
                                        "\\[|\\]|'|\\s")) %>% 
  separate_rows(run_accession, sep = ",") %>%  
  filter(str_detect(run_accession, "^(ERR|SRR)")) %>% 
  distinct(fermented_food, run_accession) %>% 
  count(fermented_food, name = "n_samples") %>% 
  arrange(desc(n_samples))

grains_sample_runs <- grains_mag_sample_metadata %>% 
  select(fermented_food, run_accession) %>% 
  mutate(run_accession = str_remove_all(run_accession,
                                        "\\[|\\]|'|\\s")) %>% 
  separate_rows(run_accession, sep = ",") %>%  
  filter(str_detect(run_accession, "^(ERR|SRR)")) %>% 
  distinct(fermented_food, run_accession)

write_tsv(grains_sample_runs, "inputs/grains_sample_runs.tsv")
