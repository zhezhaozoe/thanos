#!/usr/bin/env Rscript --vanilla

# Import depths
mags_depths_files <- "inst/extdata/mags_example/depths/mag_depths_summary.tsv"
mags_otus <- import_mag_depths(mags_depths_files)

# Read sample metadata
samplesheet <- read.table(
  "inst/extdata/samplesheets/samplesheet.tsv",
  header = TRUE,
  row.names = 1
)

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
interesting_KOs <- get_kegg_kos_from_module("M00176")
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
saveRDS(mags_hits, "/g/bork3/home/marotta/thanos2/paper/mags_hits.Rds")

# Make plots
base_size <- 7
update_geom_defaults("text", list(size = base_size / .pt))

p1 <- barplot_depths(mags_hits, group = "Station", fill = "Phylum", wrap = "Gene") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_bw(base_size = base_size) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("/g/bork3/home/marotta/thanos2/paper/mags_barplot.png", width = 7, height = 4)

p2 <- barplot_depths(mags_hits$cysNC, group = "Station", fill = "Phylum", wrap = "Province", position = "fill") +
  ggtitle("cysNC") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_bw(base_size = base_size) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("/g/bork3/home/marotta/thanos2/paper/mags_barplot_cysnc.png")

p3 <- boxplot_depths(mags_hits$cysNC, fill = "Province", signif = TRUE, show.legend = FALSE) +
  ggtitle("cysNC") +
  expand_limits(y = 6.1) +
  theme_bw(base_size = base_size) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("/g/bork3/home/marotta/thanos2/paper/mags_boxplot.png")

p <- ((p1 / (p3 | p2)) & theme(legend.key.size = unit(0.4, "cm"))) +
  plot_layout(guides = "collect", heights = c(4, 2.5)) +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")
ggsave("/g/bork3/home/marotta/thanos2/paper/mags_patchwork.png", width = 8, height = 6)

