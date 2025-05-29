library(tidyverse)
library(pheatmap)

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

# top genomes across samples
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

# top genomes within ingredient groups
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

# heatmap for top genomes within ingredient groups 
ingredient_group_counts <- sylph_profiles_metadata %>%
  distinct(accession_name, ingredient_group) %>%
  count(ingredient_group, name = "total_samples") %>%
  mutate(
    group_label = paste0(ingredient_group, " (n=", total_samples, ")")
  )

top_genomes_by_group_labeled <- top_genomes_by_group %>%
  left_join(ingredient_group_counts, by = "ingredient_group")

select_ingredient_groups <- c("Dairy", "Grain", "Vegetables_Aromatics", "Legumes", "Sugar", "Botanicals")

top_genomes_samples_plot <- top_genomes_by_group_labeled %>%
  filter(ingredient_group %in% select_ingredient_groups) %>%
  mutate(prop_detected = n / total_samples) %>% 
  group_by(species) %>%
  filter(sum(n, na.rm = TRUE) >= 10) %>%
  ungroup() %>%
  ggplot(aes(x = group_label, y = species, fill = prop_detected)) +
  geom_tile(color = "black") +
  scale_fill_viridis_c(
    name = "% of Samples Detected",
    option = "C",
    labels = scales::percent_format(accuracy = 1),
    limits = c(0.01, 1)
  ) +
  labs(
    x = "Ingredient Group (with Total Samples)",
    y = "Species",
    title = "Species Prevalence Across Ingredient Groups",
    subtitle = "Only species detected in ≥10 samples overall"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    axis.line = element_blank(), 
    panel.border = element_blank(),   
    axis.ticks = element_blank() 
  )

ggsave("figures/top-genomes-samples.png", top_genomes_samples_plot, width=7, height=11, units=c("in"))

# comparing abundance of foods in dairy vs grains
grain_dairy_abundance <- sylph_profiles_metadata %>%
  filter(ingredient_group %in% c("Grain", "Dairy")) %>%
  left_join(
    rep_mags_metadata %>%
      select(genome_accession, completeness, contamination, contigs, taxonomy, species),
    by = "genome_accession"
  ) %>%
  mutate(sample_id = paste0(food_name, "_", accession_name)) %>%
  filter(!is.na(species))

species_to_keep <- grain_dairy_abundance %>%
  filter(Sequence_abundance > 1) %>% 
  distinct(sample_id, ingredient_group, species) %>%
  count(ingredient_group, species, name = "n_samples") %>%
  filter(n_samples >= 7) %>%
  pull(species) %>%
  unique()

filtered_abundance <- grain_dairy_abundance %>%
  filter(species %in% species_to_keep)

abundance_matrix <- filtered_abundance %>%
  group_by(species, sample_id) %>%
  summarise(Sequence_abundance = sum(Sequence_abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = sample_id,
    values_from = Sequence_abundance,
    values_fill = 0
  ) %>%
  column_to_rownames("species")

sample_annotations <- filtered_abundance %>%
  distinct(sample_id, ingredient_group) %>%
  column_to_rownames("sample_id")

ordered_samples <- sample_annotations %>%
  arrange(ingredient_group) %>% 
  rownames_to_column("sample_id") %>%
  pull(sample_id)

abundance_matrix <- abundance_matrix[, ordered_samples]

grains_dairy_comps_plot <- pheatmap(
  mat = abundance_matrix,
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = sample_annotations,
  color = viridis::mako(100),
  show_colnames = TRUE,
  fontsize_col = 6,
  main = "Species Abundance (≥7 Samples per Ingredient Group)",
  na_col = "grey90"
)

ggsave("figures/grains-dairy-species-abundance-comps.png", grains_dairy_comps_plot, width=15, height=8, units=c("in"))
