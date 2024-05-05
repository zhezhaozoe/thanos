#!/usr/bin/env Rscript --vanilla

library(thanos)
library(patchwork)

# Import depths
mags_depths_files <- "inst/extdata/mags_example/depths/mag_depths_summary.tsv"
mags_otus <- import_mag_depths(mags_depths_files)

# Read sample metadata
samplesheet <- read.table(
  "inst/extdata/samplesheets/samplesheet.tsv",
  header = TRUE,
  row.names = 1
)
samplesheet$Depth <- factor(samplesheet$Depth, ordered = T)

# Read taxonomy
gtdb_taxonomy <- read_gtdbtk("inst/extdata/mags_example/taxonomy/gtdbtk_summary.tsv")

# Build a phyloseq object from the OTUs with sample data and tax table
mags_ps <- phyloseq(
  mags_otus,
  sample_data(samplesheet),
  tax_table(gtdb_taxonomy)
)

# Path to the control gene profile
control_hmm <- "inst/extdata/controls/bac120_r214_reps_PF01025.20.hmm"

# Generate query profiles
kegg_module <- "M00176"
interesting_KOs <- get_kegg_kos_from_module(kegg_module)
queries_hmm <- build_hmm_from_ko(interesting_KOs, nmax = 15)

# List of protein sequences files
mags_sequences_files <- list.files("inst/extdata/mags_example/protein_sequences", full.names = TRUE)
names(mags_sequences_files) <- sub("\\.faa$", "", basename(mags_sequences_files))

# Get hits depths
mags_hits <- get_hits_depths_from_hmm(
  queries_hmm,
  control_hmm,
  mags_ps,
  mags_sequences_files,
  linker = mags_linker,
  taxrank = "Phylum",
  parallel_processes = 48
)
# saveRDS(mags_hits, "mags_hits.Rds")

# Make plots
base_size <- 7
update_geom_defaults("text", list(size = base_size / .pt))

p1 <- barplot_depths(mags_hits, group = "Station", fill = "Phylum", wrap = "Gene") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_bw(base_size = base_size) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("../figures/mags_barplot.png", width = 7, height = 4)

p2 <- barplot_depths(mags_hits$cysNC, group = "Station", fill = "Phylum", wrap = "Province", position = "fill") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_bw(base_size = base_size) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("../figures/mags_barplot_cysnc.png")

p3 <- boxplot_depths(mags_hits$cysNC, fill = "Province", signif = TRUE, show.legend = FALSE) +
  expand_limits(y = 6.1) +
  theme_bw(base_size = base_size) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("../figures/mags_boxplot.png")

p <- ((p1 / (p3 | p2)) & theme(legend.key.size = unit(0.4, "cm"))) +
  plot_layout(guides = "collect", heights = c(4, 2.5)) +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")
ggsave("../figures/mags_patchwork.png", width = 8, height = 6)


p_facetgrid <- barplot_depths(mags_hits, group = "Station", fill = "Phylum", wrap = c("Gene", "Province")) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_bw(base_size = base_size) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("../figures/mags_barplot_facetgrid.png", width = 3, height = 7)


p_station_depth <- barplot_depths(mags_hits, group = "Depth", fill = "Phylum", wrap = c("Gene", "Province")) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_bw(base_size = base_size) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("../figures/mags_barplot_station_depth.png", width = 3, height = 7)

p_module <- keggmodule_plot(kegg_module, setNames(mags_hits, interesting_KOs))
ggsave("../figures/mags_keggmodule.png", width = 8, height = 6)
