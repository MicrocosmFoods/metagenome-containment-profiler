library(tidyverse)

# Find fermented grains samples to profile metagenomes/genome abundance

# list of fermented grains foods
grain_foods <- read_tsv("metadata/grains-foods-list.txt", col_names = c("food")) %>% 
  mutate(fermented_food = tolower(food)) %>% 
  select(-food)

food_list <- grain_foods %>% 
  pull(fermented_food)

# mag/sample metadata
# original metadata
mag_sample_metadata <- read_csv("metadata/Food_MAGs_curated_metadata_250421_corrected_merged_final_v2_corrected.csv")

# metadata filtered with grains foods hits
grains_mag_sample_metadata <- mag_sample_metadata %>% 
  filter(fermented_food %in% food_list)

grains_mag_sample_metadata %>% 
  group_by(fermented_food, run_accession) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  print(n=150)

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

# profiling results
grains_profiling_results_batch1 <- read_tsv("results/2025-05-14-updated-grains-profiles/combined_large_batch_profiles.tsv") %>% 
  mutate(run_accession = gsub("_trimmed_1.fastq.gz", "", Sample_file)) %>% 
  mutate(genome_accession = gsub(".fna", "", Genome_file))

grains_profiling_results_batch2 <- read_tsv("results/2025-05-14-updated-grains-profiles/combined_sylph_profiles.tsv") %>% 
  mutate(run_accession = sample_name) %>% 
  select(-sample_name) %>% 
  mutate(genome_accession = gsub(".fna", "", Genome_file))

combined_grains_sylph_results <- rbind(grains_profiling_results_batch1, grains_profiling_results_batch2)

strain_metadata <- read_tsv("inputs/industrial-strain-accessions.tsv") %>% 
  mutate(genome_accession = accession) %>% 
  select(genome_name, genome_accession)

grains_profiling_metadata <- left_join(combined_grains_sylph_results, grains_sample_runs) %>% 
  left_join(strain_metadata) %>% 
  mutate(food_sample = paste0(fermented_food, " (", run_accession, ")")) %>% 
  filter(Adjusted_ANI > 98)

grains_profiling_metadata %>% 
  ggplot(aes(x=food_sample, y=genome_name)) + 
  geom_tile(aes(fill=Sequence_abundance)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

grains_profiling_metadata %>% 
  filter(fermented_food != "sourdough") %>% 
  ggplot(aes(x=food_sample, y=genome_name)) + 
  geom_tile(aes(fill=Sequence_abundance)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
