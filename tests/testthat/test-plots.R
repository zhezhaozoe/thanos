setwd("/g/bork3/home/marotta/zzthanos")
devtools::load_all()

mag_depths_file <- "/g/scb/bork/marotta/thanos/GenomeBinning/depths/bins/bin_depths_summary.tsv"
seq_files <- Sys.glob("/g/scb/bork/marotta/thanos/Annotation/Prokka/*.faa")
names(seq_files) <- sub("\\.faa", "", basename(seq_files))
samplesheet_file <- "/g/bork3/home/marotta/samplesheet_v2.csv"
gtdb_file <- "/g/scb/bork/marotta/thanos/Taxonomy/GTDB-Tk/gtdbtk_summary.tsv"

samplesheet <- read.csv(samplesheet_file, header = T, row.names = 1)
sunagawa <- setDT(readxl::read_xlsx("/g/bork3/home/marotta/sunagawa_tables1.xlsx"))
samplesheet$ERR <- sub(".*_(ERR[^_]*)_.*", "\\1", samplesheet$short_reads_1)
samplesheet$ERS <- sapply(samplesheet$ERR, function(err) {
  sunagawa[grep(err, `INSDC run accession number(s)`)][[2]]
})
samplesheet_raw <- merge(samplesheet, sunagawa, by.x = "ERS", by.y = names(sunagawa)[2])
setDF(samplesheet_raw, rownames = sub(".*(ERX[^_]*)_.*", "\\1", samplesheet_raw$short_reads_1))
setnames(samplesheet_raw, c("Station identifier [TARA_station#]", "Latitude [degrees North]", "Longitude [degrees East]", "Sampling depth [m]", "Marine provinces  (Longhurst 2007)"), c("Station", "Latitude", "Longitude", "Depth", "Province"))
samplesheet_raw$Province <- sub("\\((.*)\\).*", "\\1", samplesheet_raw$Province)

test_that("the plots look like Michelangelo's paintings", {
  samplesheet <- sample_data(samplesheet_raw)
  gtdb <- tax_table(read_gtdbtk(gtdb_file))
  otus <- import_mag_depths(mag_depths_file)
  ps <- phyloseq(otus, samplesheet, gtdb)
  query_tblout <- fread("query_tblout_mags.tsv")
  control_tblout <- readRDS("gtdb_markers_vs_mags.tblout.Rds")[[4]]
  nullrank <- get_hits_depths(ps, query_tblout, control_tblout, mags_linker, taxrank = NULL) |> order_samples_by("Province")
  familyrank <- get_hits_depths(ps, query_tblout, control_tblout, mags_linker, taxrank = "Phylum") |> order_samples_by("Province")

  barplot_depth_by_sample(nullrank, fill = "Province") + facet_wrap(~Province, scale = "free_x")

  barplot_depth_by_sample(familyrank, fill = "Phylum") + facet_wrap(~Phylum)

  barplot_depth(familyrank, fill = "Province", wrap = "Phylum")
  barplot_depth(familyrank, fill = "Phylum", wrap = "Province") + facet_wrap(~Province, scale = "free_x")
  barplot_depth(familyrank, fill = "Province", wrap = NULL)
  barplot_depth(familyrank, fill = "Phylum", wrap = NULL)
  barplot_depth(familyrank, fill = NULL, wrap = "Phylum")
  barplot_depth(familyrank, fill = NULL, wrap = NULL)

  boxplot_depth(nullrank, "Province")

  boxplot_depth(familyrank, fill = "Province", wrap = "Phylum")

  ggsave("tmp.png")
})
