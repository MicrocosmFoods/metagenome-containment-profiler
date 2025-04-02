library(tidyverse)

# all curated MAG metadata
curated_mag_metadata <- read_tsv("metadata/2025-03-24-all-ff-mag-metadata-cleaned-curated.tsv")

# all Carlino 2024 cFMD v1.0.0 sample metadata
carlino_sample_metadata <- read_tsv("metadata/2024-11-04-Carlino-sample-metadata.tsv") %>% 
  mutate(source = paste0(dataset_name, "__", sample_id))

# curated Carlino 2024 MAGs joined with sample metadata
carlino_mag_sample_info <- curated_mag_metadata %>% 
  filter(study_catalog == "Carlino2024") %>% 
  left_join(carlino_sample_metadata)

filtered_carlino_samples <- carlino_mag_sample_info %>% 
  filter(library_layout == "PAIRED") %>% 
  filter(sequencing_platform == "Illumina_NovaSeq_6000" | sequencing_platform == "Illumina_HiSeq_4000" | sequencing_platform == "Illumina_HiSeq_2500" | sequencing_platform == "NextSeq_500" | sequencing_platform == "Illumina_HiSeq_2000" | sequencing_platform == "Illumina_HiSeq_1500") %>% 
  filter(database_origin == "NCBI" | database_origin == "ENA") %>% 
  filter(n_of_reads > 10000000)

filtered_fermented_foods <- filtered_carlino_samples %>% 
  group_by(fermented_food) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  pull(fermented_food)

filtered_food_samples <- filtered_carlino_samples %>% 
  filter(fermented_food %in% filtered_fermented_foods) %>% 
  filter(completeness > 90 & contamination < 10)

# random sampling of at least 3 samples per fermented food
selected_samples <- filtered_food_samples %>%
  group_by(fermented_food) %>%
  group_modify(~ {
    n_samples <- min(nrow(.x), 3)  # take 3 or all samples if less than 3
    slice_sample(.x, n = n_samples)
  }) %>%
  ungroup()

selected_samples %>%
  group_by(fermented_food) %>%
  summarise(n_samples = n()) %>%
  arrange(desc(n_samples)) %>%
  print(n = Inf)

selected_samples_df <- selected_samples  %>% 
  select(run_accession, fermented_food, substrate_category, general_category, source)

selected_samples_input_table <- selected_samples_df %>% 
select(fermented_food, run_accession)  %>% 
mutate(sample_name = fermented_food)  %>% 
mutate(accession = run_accession)  %>% 
select(sample_name, accession)

write_tsv(selected_samples_input_table, "inputs/2025-04-03-selected-carlino-samples-input-table.tsv")
