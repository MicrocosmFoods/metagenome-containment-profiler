library(tidyverse)

#################################
# Prep metadata and sylph profiles files
#################################

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

# sylph profiling results
sylph_profiles <- read_tsv("results/2025-12-01-zymo-grains-profiling/combined_sylph_profiles.tsv") %>%
  mutate(accession_name = gsub("_trimmed_1.fastq.gz", "", Sample_file)) %>% 
  mutate(genome_accession = gsub(".fa", "", Genome_file)) %>% 
  select(accession_name, genome_accession, Sequence_abundance, Adjusted_ANI, Eff_cov, Contig_name)

# sample metadata
sample_metadata <- read.csv("results/2025-12-01-zymo-grains-profiling/2025-12-01-zymo-oat-sequencing-metadata.csv") %>% 
  mutate(accession_name = gsub("_R1.fastq.gz", "", fastq_1)) %>% 
  select(sample, accession_name)

# merge with genome and sample metadata
sylph_profiles_metadata <- left_join(sylph_profiles, rep_mags_metadata) %>% 
  left_join(sample_metadata)

#################################
# Basic summary stats
#################################

# summary stats per sample
sylph_profiles_stats <- sylph_profiles_metadata %>% 
  group_by(sample) %>%
  summarise(
    n_genomes = n_distinct(genome_accession),
    percent_mapped = round(sum(Sequence_abundance, na.rm = TRUE), 3),
    percent_unmapped = round(100 - sum(Sequence_abundance, na.rm = TRUE), 3),
    .groups = "drop"
  )
