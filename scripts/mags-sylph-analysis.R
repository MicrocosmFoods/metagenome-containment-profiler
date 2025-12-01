library(tidyverse)
library(pheatmap)
library(UpSetR)
library(ggridges)
library(ape)
library(ggtree)
library(patchwork)
library(cowplot)
library(circlize)

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
sylph_profiles_mag_metadata <- left_join(sylph_profiles, rep_mags_metadata) %>% 
  mutate(domain = gsub(";.*", "", taxonomy))

sylph_profiles_full_metadata_df <- sylph_profiles_mag_metadata %>% 
  left_join(food_metadata, by="accession_name") %>% 
  select(accession_name, food_name, main_ingredient, ingredient_group, genome_accession, species, Sequence_abundance, Adjusted_ANI, Eff_cov, completeness, contamination, contigs, taxonomy)

write_tsv(sylph_profiles_full_metadata_df, "results/all-sylph-profiles-results.tsv")

# filter to just eukaryotic hits
euk_sylph_profiles_metadata <- sylph_profiles_mag_metadata %>% 
  filter(domain %in% c("Ascomycota", "Basidiomycota", "Mucoromycota")) %>% 
  left_join(food_metadata, by="accession_name") %>% 
  select(accession_name, food_name, genome_accession, species, Sequence_abundance, Adjusted_ANI, Eff_cov, completeness, contamination, contigs, taxonomy)

write_tsv(euk_sylph_profiles_metadata, "results/euk-sylph-profiles-results.tsv")

#################################
# Basic summary stats
#################################

# summary stats per sample
sylph_profiles_stats <- sylph_profiles_metadata %>% 
  group_by(accession_name) %>%
  summarise(
    n_genomes = n_distinct(genome_accession),
    percent_mapped = round(sum(Sequence_abundance, na.rm = TRUE), 3),
    percent_unmapped = round(100 - sum(Sequence_abundance, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  left_join(sylph_profiles_metadata %>% select(accession_name, food_name, main_ingredient, ingredient_group) %>% distinct(), by = "accession_name")

missing_samples <- food_metadata %>%
  anti_join(sylph_profiles_metadata, by = "accession_name")

# genome stats
genome_stats <- sylph_profiles_metadata %>%
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

#################################
# Plot viz
#################################

# heatmap for top genomes within ingredient groups 
ingredient_group_counts <- sylph_profiles_metadata %>%
  distinct(accession_name, ingredient_group) %>%
  count(ingredient_group, name = "total_samples") %>%
  mutate(
    group_label = paste0(ingredient_group, " (n=", total_samples, ")")
  )

top_genomes_by_group_labeled <- top_genomes_by_group %>%
  left_join(ingredient_group_counts, by = "ingredient_group")

select_ingredient_groups <- c("Dairy", "Grain", "Vegetables_Aromatics", "Legumes", "Sugar", "Botanical", "Fruit", "Roots_Tubers")

top_genomes_list_filtered <- top_genomes_by_group_labeled %>% 
  filter(!is.na(species)) %>% 
  filter(n >= 10) %>%
  pull(genome_accession)

# top samples heatmap plot
top_genomes_samples_plot <- top_genomes_by_group_labeled %>%
  filter(ingredient_group %in% select_ingredient_groups) %>%
  mutate(prop_detected = n / total_samples) %>% 
  group_by(species) %>%
  filter(sum(n, na.rm = TRUE) >= 10) %>%
  ungroup() %>%
  ggplot(aes(x = group_label, y = fct_rev(species), fill = prop_detected)) +
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
    subtitle = "Only species detected in ≥10 samples"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.text.y = element_text(face="italic"),
    title = element_text(face="bold", size=14),
    legend.position = "right",
    axis.line = element_blank(), 
    panel.border = element_blank(),   
    axis.ticks = element_blank() 
  ) +
  scale_x_discrete(position = "bottom") + 
  theme(rect = element_rect(fill = "transparent"))

# ridges plot of top genomes and abundance
top_genomes_abundance_data <- sylph_profiles_metadata %>% 
  filter(genome_accession %in% top_genomes_list_filtered,
         ingredient_group %in% select_ingredient_groups) %>% 
  left_join(rep_mags_metadata %>% select(genome_accession, species), by="genome_accession") %>% 
  filter(!is.na(species)) %>% 
  distinct(accession_name, genome_accession, species, ingredient_group, Sequence_abundance)

top_genomes_ridges_plot <- ggplot(top_genomes_abundance_data, aes(x = Sequence_abundance, y = fct_rev(species), fill = ingredient_group)) +
  geom_density_ridges(
    alpha = 0.7,
    scale = 1,
    rel_min_height = 0.05
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Sequence Abundance",
    y = "Species",
    fill = "Ingredient Group",
    title = "Abundance Distribution of Top Genomes Across Ingredient Groups", 
    subtitle = "Only species detected in ≥10 samples within individual ingredient groups"
  ) +
  theme_ridges() +
  theme(legend.position = "right") +
  theme(axis.text.y=element_text(face="italic"))

# save heatmap and ridges plots
ggsave("figures/top-genomes-samples.png", top_genomes_samples_plot, width=8, height=11, units=c("in"))
ggsave("figures/top-genomes-ridges.png", top_genomes_ridges_plot, width=11, height=8, units=c("in"))

# split main abundance ridges plot by bacteria vs eukaryotes with corresponding trees plotted
# first split out bacteria and eukaryotes
top_genomes_tax_info <- top_genomes_abundance_data %>% 
  select(genome_accession) %>% 
  unique() %>% 
  left_join(rep_mags_metadata) %>% 
  filter(genome_accession != "LeechJ_xxxx__AF51__bin.9") %>% 
  mutate(phylum = sapply(strsplit(as.character(taxonomy), ";"), `[`, 1))

euk_df <- top_genomes_tax_info %>% 
  filter(phylum == "Ascomycota")

bac_df <- top_genomes_tax_info %>% 
  filter(phylum != "Ascomycota")

# write out list of bacterial genome accessions to make tree
write.table(bac_df$genome_accession, "metadata/bacterial_genome_accessions.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# species names with updated genome accession formatting
species_names <- top_genomes_abundance_data %>% 
  mutate(genome_accession = gsub("-", "_", genome_accession)) %>% 
  filter(genome_accession != "LeechJ_xxxx__AF51__bin.9") %>% 
  select(genome_accession, species) %>% 
  unique()

# read in the bacterial concatenated tree and plot
bac_tree <- read.tree("results/2025-05-28-mags-profiling/bacteria-fastTree-ribosomal-tree.tre")

label_map <- setNames(species_names$species, species_names$genome_accession)

bac_tree$tip.label <- label_map[bac_tree$tip.label]

tree_plot_obj <- ggtree(bac_tree, branch.length = "none") +
  geom_tiplab(
    align = TRUE,
    linetype = "dotted",
    size = 4,
    linesize = 0.3
  ) +
  xlim_tree(30)


tree_plot_obj

tip_order <- tree_plot_obj$data %>% 
  filter(isTip) %>% 
  arrange(y) %>% 
  pull(label)

bacterial_plot_data <- top_genomes_abundance_data %>%
  filter(species %in% bac_tree$tip.label) %>%
  mutate(species = factor(species, levels = tip_order))

ordered_bacterial_ridges_plot <- ggplot(bacterial_plot_data, aes(x = Sequence_abundance, y = species, fill = ingredient_group)) +
  geom_density_ridges(alpha = 0.7, scale = 1, rel_min_height = 0.05) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    expand = c(0,0)
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Sequence Abundance",
    y = NULL,
    fill = "Ingredient Group"
  ) +
  theme_ridges() +
  theme(
    axis.text.y = element_text(face="italic", size=14),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(t = 5, r = 10, b = 5, l = 80)
  )

# eukaryotic ridges plot
euks_list <- euk_df %>% 
  mutate(genome_accession = gsub("-", "_", genome_accession)) %>% 
  pull(genome_accession)
  

euk_abundance_data <- top_genomes_abundance_data %>% 
  filter(genome_accession %in% euks_list)

euk_ridges_plot <- ggplot(euk_abundance_data, aes(x = Sequence_abundance, y = species, fill = ingredient_group)) +
  geom_density_ridges(alpha = 0.7, scale = 1, rel_min_height = 0.05) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20)
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Sequence Abundance",
    y = NULL,
    fill = "Ingredient Group"
  ) +
  theme_ridges() +
  theme(
    axis.text.y = element_text(face="italic"),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(l = -20)
  )

# save individual plots
# tree
ggsave("figures/core-bac-genomes-tree.png", tree_plot_obj)

# bacterial ridges plot
ggsave("figures/ordered-bacterial-ridges-plot.png", ordered_bacterial_ridges_plot, width=15, height=10, units=c("in"))

# euk ridges plot
ggsave("figures/euk-ridges-plot.png", euk_ridges_plot, width=15, height=8, units=c("in"))

# Upset plot of overlap of number of genomes across ingredient groups

sylph_binary <- sylph_profiles_metadata %>%
  select(genome_accession, ingredient_group) %>%
  filter(ingredient_group %in% select_ingredient_groups) %>% 
  distinct() %>%  # Remove any duplicates
  mutate(present = 1) %>%
  pivot_wider(
    names_from = ingredient_group,
    values_from = present,
    values_fill = list(present = 0)
  )
  
upset_input <- as.data.frame(sylph_binary %>% select(-genome_accession))

upset_plot <- upset(upset_input, 
      sets = colnames(upset_input),
      order.by = "freq")

png("figures/upset-plot-genomes.png", width = 2500, height = 1600, res = 300)
upset(upset_input, 
      sets = colnames(upset_input),
      order.by = "freq")
dev.off()


# heatmap of abundance of genomes that are in all select ingredient groups
genomes_in_all <- sylph_binary %>%
  filter(if_all(all_of(select_ingredient_groups), ~ . == 1)) %>%
  pull(genome_accession)

abundance_subset <- sylph_profiles_metadata %>%
  filter(genome_accession %in% genomes_in_all,
         ingredient_group %in% select_ingredient_groups) %>%
  left_join(rep_mags_metadata %>% select(genome_accession, species), by = "genome_accession") %>%
  mutate(sample_id = paste0(food_name, "_", accession_name)) %>%
  filter(!is.na(species))

species_labels <- abundance_subset %>%
  distinct(species, genome_accession) %>%
  group_by(species) %>%
  mutate(
    species_unique = if (n() > 1) {
      paste0(species, " (", LETTERS[seq_len(n())], ")")
    } else {
      species
    }
  ) %>%
  ungroup()

abundance_subset <- abundance_subset %>%
  left_join(species_labels, by = c("species", "genome_accession")) %>%
  mutate(species = species_unique) %>%
  select(species, sample_id, ingredient_group, Sequence_abundance)

abundance_matrix <- abundance_subset %>%
  group_by(species, sample_id) %>%
  summarise(Sequence_abundance = sum(Sequence_abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = sample_id,
    values_from = Sequence_abundance,
    values_fill = list(Sequence_abundance = 0)
  ) %>%
  column_to_rownames("species")

nonzero_samples <- colnames(abundance_matrix)[colSums(abundance_matrix) > 0]
abundance_matrix <- abundance_matrix[, nonzero_samples]

sample_annotations <- abundance_subset %>%
  distinct(sample_id, ingredient_group) %>%
  filter(sample_id %in% nonzero_samples) %>%
  column_to_rownames("sample_id")

ordered_samples <- sample_annotations %>%
  arrange(ingredient_group) %>%
  rownames_to_column("sample_id") %>%
  pull(sample_id)

abundance_matrix <- abundance_matrix[, ordered_samples]
sample_annotations <- sample_annotations[ordered_samples, , drop = FALSE]


# color palettes
heatmap_colors <- viridis::mako(100)

group_colors <- RColorBrewer::brewer.pal(length(select_ingredient_groups), "Set2")
names(group_colors) <- select_ingredient_groups
annotation_colors <- list(ingredient_group = group_colors)

top_abundance_heatmap <- pheatmap(
  mat = as.matrix(abundance_matrix),
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = sample_annotations,
  annotation_colors = annotation_colors,
  show_colnames = TRUE,
  fontsize_col = 6,
  main = "Abundance of Shared Species Across Samples",
  color = heatmap_colors,
  na_col = "grey90")

ggsave("figures/top-genomes-abundance-heatmap.png", top_abundance_heatmap, width=25, height=8, units=c("in"))


# comparing abundance of genomes in veg vs grains
veg_dairy_abundance <- sylph_profiles_metadata %>%
  filter(ingredient_group %in% c("Dairy", "Vegetables_Aromatics")) %>%
  left_join(
    rep_mags_metadata %>%
      select(genome_accession, completeness, contamination, contigs, taxonomy, species),
    by = "genome_accession"
  ) %>%
  mutate(sample_id = paste0(food_name, "_", accession_name)) %>%
  filter(!is.na(species))

species_to_keep <- veg_dairy_abundance %>%
  filter(Sequence_abundance > 1) %>% 
  distinct(sample_id, ingredient_group, species) %>%
  count(ingredient_group, species, name = "n_samples") %>%
  filter(n_samples >= 7) %>%
  pull(species) %>%
  unique()

filtered_abundance <- veg_dairy_abundance %>%
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

veg_dairy_comps_plot <- pheatmap(
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

ggsave("figures/veg-dairy-species-abundance-comps.png", veg_dairy_comps_plot, width=15, height=8, units=c("in"))

# ridge plot of genomes found in all ingredient groups
genomes_in_all <- sylph_binary %>%
  filter(if_all(all_of(select_ingredient_groups), ~ . == 1)) %>%
  pull(genome_accession)

abundance_data <- sylph_profiles_metadata %>%
  filter(genome_accession %in% genomes_in_all,
         ingredient_group %in% select_ingredient_groups) %>%
  left_join(rep_mags_metadata %>% select(genome_accession, species), by = "genome_accession") %>%
  filter(!is.na(species)) %>%
  distinct(accession_name, genome_accession, species, ingredient_group, Sequence_abundance) %>% 
  mutate(species = case_when(
    genome_accession == "MAG_116" ~ paste0(species, " (A)"),
    genome_accession == "MortensenS_xxxx__WK023-1-L-48h-02G3_S117__bin.1" ~ paste0(species, " (B)"),
    TRUE ~ species
  ))

abundance_data_indv_food <- sylph_profiles_metadata %>%
  left_join(rep_mags_metadata %>% select(genome_accession, species), by = "genome_accession") %>%
  filter(!is.na(species)) %>%
  distinct(accession_name, genome_accession, species, food_name, Sequence_abundance) %>% 
  mutate(species = case_when(
    genome_accession == "MAG_116" ~ paste0(species, " (A)"),
    genome_accession == "MortensenS_xxxx__WK023-1-L-48h-02G3_S117__bin.1" ~ paste0(species, " (B)"),
    TRUE ~ species
  ))

core_genomes_ridges <- ggplot(abundance_data, aes(x = Sequence_abundance, y = fct_rev(species), fill = ingredient_group)) +
  geom_density_ridges(
    alpha = 0.7,
    scale = 1,
    rel_min_height = 0.05
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Sequence Abundance",
    y = "Species",
    fill = "Ingredient Group",
    title = "Abundance Distribution of Shared Genomes Across Ingredient Groups"
  ) +
  theme_ridges() +
  theme(legend.position = "bottom",
        axis.text.y=element_text(face="italic", size=14))

core_genomes_ridges

ggsave("figures/core-genomes-ridges-plot.png", core_genomes_ridges, width=11, height=8, units=c("in"))

