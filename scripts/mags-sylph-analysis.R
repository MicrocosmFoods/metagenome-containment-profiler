library(tidyverse)

# MAG metadata from curation repo
mag_metadata_url <- "https://raw.githubusercontent.com/MicrocosmFoods/fermentedfood_metadata_curation/refs/heads/main/data/2025-05-21-genome-metadata-food-taxonomy.tsv"

mag_metadata <- read_tsv(mag_metadata_url) %>% 
  mutate(genome_accession = mag_id) %>% 
  select(genome_accession, completeness, contamination, contigs, taxonomy, species, rep_95id, food_name, main_ingredient, ingredient_group, origin, food_type)

rep_mags_metadata <- mag_metadata %>% 
  filter(genome_accession == rep_95id) %>% 
  select(-rep_95id) %>% 
  mutate(species = case_when(
    is.na(species) | str_to_lower(species) == "unknown" ~ str_c(
      str_extract(taxonomy, "[^;]+$"),
      " spp."
    ),
    TRUE ~ species
  )) %>% 
  select(genome_accession, completeness, contamination, contigs, taxonomy, species)

# food taxonomy
food_taxonomy <- mag_metadata %>% 
  select(food_name, main_ingredient, ingredient_group, food_type) %>% 
  distinct()

# food metadata
food_metadata <- read_tsv("inputs/2025-05-22-selected-samples.tsv") %>% 
  mutate(food_name = sample_name) %>% 
  mutate(accession_name = accession) %>% 
  select(food_name, accession_name) %>% 
  distinct(accession_name, .keep_all = TRUE) %>% 
  left_join(food_taxonomy) %>% 
  distinct(food_name, accession_name, .keep_all = TRUE)

# sylph profiling results
sylph_profiles <- read_tsv("results/2025-05-28-mags-profiling/2025-05-28-mags-rep-samples-profiles.tsv") %>% 
  mutate(accession_name = gsub("_trimmed_1.fastq.gz", "", Sample_file)) %>% 
  mutate(genome_accession = gsub(".fa", "", Genome_file)) %>% 
  select(accession_name, genome_accession, Sequence_abundance, Adjusted_ANI, Eff_cov, Contig_name)

# join profiling results with metadata 
# these are results by default at 95% ANI so don't need to filter further since doing "species" profiling

sylph_profiles_metadata <- left_join(sylph_profiles, food_metadata, by="accession_name")

# summary stats per sample
sylph_profiles_stats <- sylph_profiles_metadata %>% 
  group_by(accession_name) %>%
  summarise(
    n_genomes = n_distinct(genome_accession),
    percent_mapped = sum(Sequence_abundance, na.rm = TRUE),
    percent_unmapped = round(100 - sum(Sequence_abundance, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  left_join(sylph_profiles_metadata %>% select(accession_name, food_name, main_ingredient, ingredient_group) %>% distinct(), by = "accession_name")

missing_samples <- food_metadata %>%
  anti_join(sylph_profiles_metadata, by = "accession_name")

# top genomes across samples and within ingredient groups 
top_genomes_overall <- sylph_profiles_metadata %>%
  filter(Sequence_abundance > 1) %>%
  distinct(accession_name, genome_accession) %>%
  count(genome_accession, sort = TRUE) %>%
  rename(n_samples = n) %>%
  left_join(
    rep_mags_metadata %>%
      select(genome_accession, completeness, contamination, contigs, taxonomy, species),
    by = "genome_accession"
  ) %>%
  left_join(
    sylph_profiles_metadata %>%
      filter(Sequence_abundance > 1) %>%
      group_by(genome_accession) %>%
      summarise(
        min_abundance = min(Sequence_abundance, na.rm = TRUE),
        median_abundance = median(Sequence_abundance, na.rm = TRUE),
        max_abundance = max(Sequence_abundance, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "genome_accession"
  )

top_genomes_by_group <- sylph_profiles_metadata %>%
  filter(Sequence_abundance > 1) %>%
  distinct(accession_name, genome_accession, ingredient_group) %>%
  count(ingredient_group, genome_accession, sort = TRUE) %>%
  group_by(ingredient_group) %>%
  arrange(desc(n), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  left_join(rep_mags_metadata %>%
              select(genome_accession, completeness, contamination, contigs, taxonomy, species),
            by = "genome_accession") %>%
  left_join(
    sylph_profiles_metadata %>%
      filter(Sequence_abundance > 1) %>%
      group_by(genome_accession) %>%
      summarise(
        min_abundance = min(Sequence_abundance, na.rm = TRUE),
        median_abundance = median(Sequence_abundance, na.rm = TRUE),
        max_abundance = max(Sequence_abundance, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "genome_accession"
  )

