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

colnames(isolate_metadata) <- c("species", "updated_name", "typical_use", "ncbi_accession", "ncbi_link", "genome_accession", "contigs", "chromosome", "notes", "kids_products", "GRAS_approved", "GRAS_link", "commercial_product")

isolate_metadata_modf <- isolate_metadata %>% 
  select(updated_name, genome_accession, contigs, typical_use, kids_products, GRAS_approved)

# joined dfs
industrial_sylph_results_info <- left_join(industrial_sylph_results, sample_metadata, by="accession_name") %>% 
  left_join(isolate_metadata_modf)

fungi_list <- c("Saccharomyces boulardii strain KCTC", "Saccharomyces cerevisiae S288C", "Aspergillus oryzae", "Aspergillus niger", "Rhizopus arrhizus strain Z10C7", "Monascus purpureus", "Sanghuangporus sanghuang", "Pleurotus ostreatus", "Penicillium roqueforti strain LCP96")

industrial_sylph_results_info_filtered <- industrial_sylph_results_info %>% 
  filter(Adjusted_ANI > 99) %>% 
  mutate(sample_code = paste0(accession_name, "_", sample_name)) %>% 
  filter(!updated_name %in% fungi_list)

industrial_sylph_results_info_filtered %>% 
  ggplot(aes(x=sample_code, y=updated_name, fill=Sequence_abundance)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# aggregate by food type

grains_list <- c("ogi", "koko", "amazake", "ugi", "fura", "brukina", "zonkom", "kunu", "massa", "boza", "burukutu", "pito", "pozol")

aggregated_results <- industrial_sylph_results_info_filtered %>% 
  group_by(food_group, updated_name) %>% 
  summarize(median_abundance = median(Sequence_abundance, na.rm = TRUE)) %>% 
  pivot_wider(names_from = food_group, values_from = median_abundance, values_fill = 0)

aggregated_ordered <- aggregated_results %>% 
  filter(updated_name %in% tree$tip.label) %>% 
  slice(match(tree$tip.label, updated_name)) %>% 
  column_to_rownames(var="updated_name")

aggregated_grains_ordered <- aggregated_ordered[, names(aggregated_ordered) %in% grains_list]
  

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
  select(updated_name, kids_products, GRAS_approved) %>% 
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


p1 <- gheatmap(p1, industrial_strain_info_ordered, offset=12, width=0.10, font.size=3, colnames_angle = 45, hjust=0, colnames_position = "top") +
  scale_fill_manual(values = c("Yes" = "steelblue", "No" = "darkgray", "Not this exact strain" = "grey90"), name = "Metadata Response") +
  theme(legend.position = "bottom")

p1 <- p1 + ggnewscale::new_scale_fill()

combined_heatmap <- gheatmap(p1, aggregated_grains_ordered, 
         offset = 14,
         width = 1, 
         font.size = 3,
         colnames_position = "top",
         colnames_angle = 45,
         hjust = 0) + 
  scale_fill_gradient(low = "ivory", high = "navyblue", name = "Median Abundance")

all_heatmap <- gheatmap(p1, aggregated_ordered,
                        offset = 14,
                        width = 1, 
                        font.size = 3,
                        colnames_position = "top",
                        colnames_angle = 45,
                        hjust = 0) + 
  scale_fill_gradient(low = "ivory", high = "navyblue", name = "Median Abundance")

all_heatmap
ggsave("figures/industrial-strain-grains-profiles.png", combined_heatmap, width=20, height=10, units=c("in"), limitsize = FALSE)
