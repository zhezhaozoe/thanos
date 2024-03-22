#' Generate a basic colData file from a mOTUs file
#'
#' @param path Path to a mOTUs file
#' @param level Taxonomic level of the mOTUs in the file
make_motus_tax_table <- function(path, level = "mOTU") {
  d <- fread(path)
  tax_table <- matrix(gsub("(^[^[]*) \\[.*\\]", "\\1", d[[1]]), ncol = 1)
  colnames(tax_table) <- level
  rownames(tax_table) <- gsub(".*\\[(.*)\\]", "\\1", d[[1]])
  tax_table
}

make_motus_tax_table_from_gtdb <- function(motus_table, gtdb_map) {
  if (!file.exists(gtdb_map)) {
    error("Please download or create a mapping between mOTUs and GTDB identifier. For mOTUs3, the file can be downloaded at https://sunagawalab.ethz.ch/share/MOTU_GTDB/mOTUs_3.0.0_GTDB_tax.tsv")
  }
  gtdb <- as.matrix(fread(gtdb_map, header = F), rownames = 1)
  gtdb_tax_table <- gtdb[rownames(motus_table), ]
  tax_table(gtdb_tax_table)
}

shorten_species_name <- function(common_name) {
  # TODO
  common_name
}

#' @import data.table
read_hmmer_tblout <- function(file) {
  all_lines <- fread(file, sep = "", showProgress = F)[[1]]
  content_lines <- grep("^#", all_lines, value = T, invert = T)
  tblout <- setDT(tstrsplit(content_lines, "\\s+")[1:18])
  names(tblout) <- c("Target", "Accession_target", "Query", "Accession_query", "Full_evalue", "Full_score", "Full_bias", "Best_Evalue", "Best_score", "Best_bias", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc")
  tblout$Description <- gsub("([^ ]+ +){18}", "", content_lines)
  type.convert(tblout, as.is = T, numerals = "no.loss")
}

#' @import data.table
read_gtdbtk <- function(file) {
  d <- fread(file, select = c("user_genome", "classification"))
  d_classification <- d[, tstrsplit(classification, ";?[dpcofgs]__")[-1]]
  setattr(d_classification, "names", c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
  as.matrix(setDF(d_classification, rownames = d$user_genome))
}
