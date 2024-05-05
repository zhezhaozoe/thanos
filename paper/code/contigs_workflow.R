#!/usr/bin/env Rscript --vanilla

library(thanos)
library(patchwork)

# Import contig depths
contigs_depths_files <- list.files("inst/extdata/contigs_example/depths/", full.names = T)
names(contigs_depths_files) <- sub("MEGAHIT-(group-\\d+)-depth.txt.gz", "\\1", basename(contigs_depths_files))
contigs_otus <- import_contig_depths(contigs_depths_files, sub_pattern = "MEGAHIT-group-\\d+-([^.]*).*", sub_replacement = "\\1", mc.cores = min(48, length(contigs_depths_files)))

# Read sample metadata
samplesheet <- read.table(
  "inst/extdata/samplesheets/samplesheet.tsv",
  header = TRUE,
  row.names = 1
)

# Build a phyloseq object from the OTUs with sample data
contigs_ps <- phyloseq(
  contigs_otus,
  sample_data(samplesheet)
)

# Read the protein sequences files
contigs_sequences_files <- list.files("inst/extdata/contigs_example/protein_sequences", full.names = T)
names(contigs_sequences_files) <- sub("\\.faa.gz", "", basename(contigs_sequences_files))

# Control hmm
control_hmm <- "inst/extdata/controls/bac120_r214_reps_PF01025.20.hmm"

# Queries hmm
glycolysis <- "M00001"
kos <- get_kegg_kos_from_module(glycolysis)
queries_hmm <- build_hmm_from_ko(kos, nmax = 10)

# Get hits
contigs_hits <- get_hits_depths_from_hmm(
  queries_hmm,
  control_hmm,
  contigs_ps,
  contigs_sequences_files,
  linker = contigs_linker,
  parallel_processes = 31,
  cpus_per_process = 2,
)
# saveRDS(contigs_hits, "contigs_hits.Rds")

# Make plots
base_size <- 7
update_geom_defaults("text", list(size = base_size / .pt))
update_geom_defaults("label", list(size = (base_size - 1) / .pt))

p1 <- keggmodule_plot(glycolysis, contigs_hits) +
  expand_limits(x = c(-0.7, 1.6)) +
  theme_void(base_size = base_size) +
  theme(legend.position = "bottom")
ggsave("../figures/contigs_graph.png")

selected_kos <- c("K00134", "K00150", "K11389", "K00927")
p2 <- barplot_depths(contigs_hits[selected_kos], group = "Station", wrap = "Gene") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_bw(base_size = base_size) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("../figures/contigs_barplot.png", width = 7, height = 4)

p <- (p1 | p2) +
  plot_layout() +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")
ggsave("../figures/contigs_patchwork.png", width = 8, height = 6)


p_facetgrid <- barplot_depths(contigs_hits, group = "Station", wrap = c("Gene", "Province")) +
#   theme_bw(base_size = base_size) +contigs_barplot_facetgrid.png", width = 6, height = 12)
