library(tidyverse)
library(BiocManager)
library(ggtree)
library(treeio)
library(ggnewscale)
library(grid)

############################
# industrial strain profiling
############################

# sylph results
industrial_sylph_results <- read_tsv("results/2025-05-27-representative-samples-industrial-strain-profiling/2025-05-27-all-profiles.tsv") %>% 
  mutate(accession_name = gsub("_trimmed_1.fastq.gz", "", Sample_file)) %>% 
  mutate(genome_accession = gsub(".fna", "", Genome_file)) %>% 
  select(accession_name, genome_accession, Sequence_abundance, Adjusted_ANI, Eff_cov, Contig_name)

# sample metadata
sample_metadata <- read_tsv("inputs/2025-05-22-selected-samples.tsv") %>% 
  unique() %>% 
  mutate(accession_name = accession) %>% 
  select(sample_name, accession_name, food_group)

# reference isolate metadata
isolate_metadata <- read_csv("metadata/Probiotic_Species_Strains_Manually_Curated.csv")

colnames(isolate_metadata) <- c("species", "updated_name", "typical_use", "genome_accession", "contigs", "notes", "cpg_products", "kids_products", "GRAS_approved", "considered_probiotic", "GRAS_link", "commercial_product")

isolate_metadata_modf <- isolate_metadata %>% 
  select(updated_name, genome_accession, contigs, typical_use, cpg_products, kids_products, GRAS_approved, considered_probiotic)

# joined dfs
industrial_sylph_results_info <- left_join(industrial_sylph_results, sample_metadata, by="accession_name") %>% 
  left_join(isolate_metadata_modf)

fungi_list <- c("Saccharomyces boulardii strain KCTC", "Saccharomyces cerevisiae S288C", "Aspergillus oryzae", "Aspergillus niger", "Rhizopus arrhizus strain Z10C7", "Monascus purpureus", "Sanghuangporus sanghuang", "Pleurotus ostreatus", "Penicillium roqueforti strain LCP96")

industrial_sylph_results_info_filtered <- industrial_sylph_results_info %>% 
  filter(Adjusted_ANI > 99) %>% 
  mutate(sample_code = paste0(accession_name, "_", sample_name)) %>% 
  filter(!updated_name %in% fungi_list) %>% 
  filter(!is.na(updated_name))

industrial_sylph_results_info_filtered %>% 
  ggplot(aes(x=sample_code, y=updated_name, fill=Sequence_abundance)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# aggregate by food type
# keep foods with at least 3 samples
top_foods <- industrial_sylph_results_info_filtered %>% 
  select(accession_name, food_group) %>% 
  distinct() %>% 
  group_by(food_group) %>% 
  count() %>% 
  filter(n >= 3) %>% 
  pull(food_group)

aggregated_results <- industrial_sylph_results_info_filtered %>% 
  filter(food_group %in% top_foods) %>% 
  group_by(food_group, updated_name) %>% 
  summarize(median_abundance = median(Sequence_abundance, na.rm = TRUE)) %>% 
  pivot_wider(names_from = food_group, values_from = median_abundance, values_fill = 0)

aggregated_ordered <- aggregated_results %>% 
  filter(updated_name %in% tree$tip.label) %>% 
  slice(match(tree$tip.label, updated_name)) %>% 
  column_to_rownames(var="updated_name")


# accessions to make corresponding phylogenetic tree from
accessions <- strain_names %>% 
  pull(genome_accession)

write_lines(accessions, "results/2025-05-27-representative-samples-industrial-strain-profiling/accessions-list.txt")

# bacterial isolate phylogenetic tree
tree <- read.tree("results/2025-05-27-representative-samples-industrial-strain-profiling/bacteria-fastTree-ribosomal-tree.tre")

strain_names <- industrial_sylph_results_info_filtered %>% 
  select(genome_accession, updated_name) %>% 
  unique()

label_map <- setNames(strain_names$updated_name, strain_names$genome_accession)

tree$tip.label <- label_map[tree$tip.label]

industrial_strain_info <- industrial_sylph_results_info_filtered %>% 
  select(updated_name, cpg_products, kids_products, GRAS_approved, considered_probiotic) %>% 
  unique()

industrial_strain_info_ordered <- industrial_strain_info %>% 
  filter(updated_name %in% tree$tip.label) %>% 
  slice(match(tree$tip.label, updated_name)) %>% 
  column_to_rownames(var = "updated_name")

# tree and metadata figure

p1 <- ggtree(tree, branch.length="none") +
  geom_tiplab(align=TRUE,
              linetype="dotted",
              size = 4,
              linesize=0.3) +
  xlim_tree(30)

p1

p1 <- gheatmap(p1, industrial_strain_info_ordered, offset=7, width=0.15, font.size=4, colnames_angle = 45, hjust=0, colnames_position = "top") +
  scale_fill_manual(values = c("Yes" = "steelblue", "No" = "darkgray", "Not this exact strain" = "grey90"), name = "Metadata Response") +
  theme(legend.position = "bottom")

p1 <- p1 + ggnewscale::new_scale_fill()

combined_heatmap <- gheatmap(p1, aggregated_ordered, 
         offset = 10,
         width = 1, 
         font.size = 4,
         colnames_position = "top",
         colnames_angle = 45,
         hjust = 0) + 
  scale_fill_gradient(low = "ivory", high = "navyblue", na.value = "ivory", name = "Median Abundance") +
  theme(legend.position = "bottom")

all_heatmap <- combined_heatmap +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(t = 60, r = 10, b = 10, l = 10),
    axis.text.x.top = element_text(margin = margin(b = 8))
  )
all_heatmap

ggsave("figures/industrial-strain-foods-profiles.png", all_heatmap, width=25, height=15, units=c("in"), limitsize = FALSE)
