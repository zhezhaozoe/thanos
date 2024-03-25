#' Generate a basic colData file from a mOTUs file
#'
#' @param path Path to a mOTUs file
#' @param level Taxonomic level of the mOTUs in the file
#' Generate Taxonomy Table from mOTUs Data
#'
#' @return A taxonomy table with one column named after the specified taxonomic level (`level`). This column contains the taxonomic names, and the row names of the table are the mOTU identifiers extracted from the data.
#'
#' @examples
#' tax_table <- make_motus_tax_table("path/to/motus_data.txt")
#' head(tax_table)
#'
#' @import data.table
#' @export
make_motus_tax_table <- function(path, level = "mOTU") {
  d <- fread(path)
  tax_table <- matrix(gsub("(^[^[]*) \\[.*\\]", "\\1", d[[1]]), ncol = 1)
  colnames(tax_table) <- level
  rownames(tax_table) <- gsub(".*\\[(.*)\\]", "\\1", d[[1]])
  tax_table
}

#' Create taxonomic table for mOTUs using GTDB identifiers
#'
#' This function takes a user's mOTUs table and a mapping file between mOTUs and GTDB
#' identifiers to create a new taxonomic table with GTDB taxonomy.
#'
#' @param motus_table The user's mOTUs table with rownames as mOTUs identifiers.
#' @param gtdb_map Path to the GTDB mapping file.
#'
#' @return A `tax_table` object containing the GTDB taxonomy for the input mOTUs.
#'
#' @examples
#' # Assuming you've already downloaded the GTDB mapping file:
#' gtdb_map_path <- "path/to/mOTUs_3.0.0_GTDB_tax.tsv"
#' motus_table <- read.csv("your_motus_table.csv", row.names = 1)
#' gtdb_tax_table <- make_motus_tax_table_from_gtdb(motus_table, gtdb_map_path)
#'
#' @export
#'
#' @note Please ensure the GTDB mapping file is available at the specified path. For mOTUs3,
#' the file can be downloaded from https://sunagawalab.ethz.ch/share/MOTU_GTDB/mOTUs_3.0.0_GTDB_tax.tsv.
#' If the file is not found, the function will stop and prompt to download or specify the mapping file.
make_motus_tax_table_from_gtdb <- function(motus_table, gtdb_map) {
  if (!file.exists(gtdb_map)) {
    error("Please download or create a mapping between mOTUs and GTDB identifier. For mOTUs3, the file can be downloaded at https://sunagawalab.ethz.ch/share/MOTU_GTDB/mOTUs_3.0.0_GTDB_tax.tsv")
  }
  gtdb <- as.matrix(fread(gtdb_map, header = F), rownames = 1)
  gtdb_tax_table <- gtdb[rownames(motus_table), ]
  tax_table(gtdb_tax_table)
}

#' Shorten Common Name to Species Initialism
#'
#' This function takes a common species name and shortens it into an initialism
#' format consisting of the first letter of the first word followed by a period
#' and the full second word. Leading and trailing whitespaces are trimmed.
#' Common long words like "Candidatus" and "species incertae sedis" are abbreviated.
#' Useful for reducing the legend clutter when plotting.
#'
#' @param common_name A string containing the common name of a species.
#' @return A string with the abbreviated species name.
#' @examples
#' shorten_species_name("Canis lupus")
#' # [1] "C. lupus"
#' shorten_species_name("Lynx")
#' # Lynx
#' shorten_species_name("Candidatus Mycoplasma hyopneumoniae")
#' # Ca. M. hyopneumoniae
#' shorten_species_name("Proteobacteria species incertae sedis")
#' # Proteobacteria sp. inc. s.
#' @export
shorten_species_name <- function(common_name) {
  common_name <- trimws(common_name) # Remove leading/trailing whitespace
  candidatus <- ""
  if (grepl("Candidatus", common_name)) {
    candidatus <- "Ca. "
    common_name <- sub("Candidatus ", "", common_name)
  }
  if (grepl("species incertae sedis", common_name)) {
    common_name <- sub("species incertae sedis", "sp. inc. s.", common_name)
    return(paste0(candidatus, common_name))
  }
  parts <- strsplit(common_name, " ")[[1]]
  if (length(parts) == 2) {
    if (parts[2] != "sp.") {
      common_name <- paste0(substr(parts[1], 1, 1), ". ", paste(parts[2], collapse = " "))
    }
    return(paste0(candidatus, common_name))
  }
  return(paste0(candidatus, common_name))
}

#' Read HMMER tblout files
#'
#' Parses the tblout output format of HMMER tools (e.g., hmmscan, hmmsearch) and
#' creates a data table with relevant columns for target and query sequences,
#' including e-values, scores, and additional metrics, plus the description.
#'
#' @param file A string. The path to the HMMER tblout file to be read.
#'
#' @return A data.table object containing the parsed HMMER tblout data. Columns
#' include Target, Accession_target, Query, Accession_query, Full_evalue,
#' Full_score, Full_bias, Best_Evalue, Best_score, Best_bias, exp, reg, clu,
#' ov, env, dom, rep, inc, and Description.
#'
#' @import data.table
#'
#' @examples
#' read_hmmer_tblout("path/to/your/hmmer_output.tblout")
#'
#' @export
read_hmmer_tblout <- function(file) {
  all_lines <- fread(file, sep = "", showProgress = F)[[1]]
  content_lines <- grep("^#", all_lines, value = T, invert = T)
  tblout <- setDT(tstrsplit(content_lines, "\\s+")[1:18])
  names(tblout) <- c("Target", "Accession_target", "Query", "Accession_query", "Full_evalue", "Full_score", "Full_bias", "Best_Evalue", "Best_score", "Best_bias", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc")
  tblout$Description <- gsub("([^ ]+ +){18}", "", content_lines)
  type.convert(tblout, as.is = T, numerals = "no.loss")
}

#' Read GTDB-Tk classification output
#'
#' Reads a GTDB-Tk (Genome Taxonomy Database Toolkit) classification file and extracts
#' the taxonomic classification into a phyloseq-friendly format. It parses the classifications
#' into domain, phylum, class, order, family, genus, and species, and returns a matrix
#' with these categories as column names and genomes as row names.
#'
#' @param file A string specifying the path to the GTDB-Tk classification output file.
#'
#' @return A matrix where rows are genomes (user_genome) and columns are the taxonomic
#' categorizations (Domain, Phylum, Class, Order, Family, Genus, Species).
#'
#' @import data.table
#'
#' @examples
#' # Assume "gtdbtk_output.tsv" is a file in the working directory with GTDB-Tk classifications
#' gtdbtk_results <- read_gtdbtk("gtdbtk_output.tsv")
#' head(gtdbtk_results)
#'
#' @export
read_gtdbtk <- function(file) {
  d <- fread(file, select = c("user_genome", "classification"))
  d_classification <- d[, tstrsplit(classification, ";?[dpcofgs]__")[-1]]
  setattr(d_classification, "names", c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
  as.matrix(setDF(d_classification, rownames = d$user_genome))
}
