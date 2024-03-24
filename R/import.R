#' Import mOTUs data
#'
#' This function is used to import mOTUs data from a specified path.
#' Optionally, it can also create a default tax table from the mOTUs names.
#'
#' @param path A string specifying the path to the motus data file.
#' @param make_tax_table A logical value indicating whether to create a default tax table. Default is FALSE.
#'
#' @description Keep only the motus ID, not the tax name, which should be provided separately as colData.
#'
#' @return The imported motus data. If make_tax_table is TRUE, a tax table is also returned.
#'
#' @import data.table
#' @import phyloseq
#' @examples
#' motus_data <- import_motus("path/to/your/motus_file.csv", make_tax_table = TRUE)
#' @export
import_motus <- function(path, make_tax_table = F) {
  d <- fread(path)
  names(d)[1] <- "ID"
  d$ID <- gsub(".*\\[(.*)\\]", "\\1", d$ID)
  m <- as.matrix(d[, !"ID"])
  rownames(m) <- d$ID
  ps <- phyloseq(otu_table(read_motus(path), taxa_are_rows = T))
  if (make_tax_table) {
    tax_table(ps) <- tax_table(make_motus_tax_table(path))
  }
}

#' Import Contig Depths from Depth Files
#'
#' This function reads depth files, processes and normalizes their contents, and then 
#' generates a contig-based OTU table suitable for downstream analysis.
#'
#' @param depth_files A character vector of file paths, each pointing to a different depth file.
#' @param sub_pattern A pattern (regular expression) to search for within non-"bam-var" column names.
#' @param sub_replacement Replacement text for any occurrences of \code{sub_pattern} within column names.
#'
#' @return An OTU table object with taxa as rows, where each row represents a contig and columns represent
#'         processed depth values adjusted for contig length and total average depth.
#'
#' @details The function operates in several steps:
#'          - It reads the depth files and combines them into a single data table.
#'          - "bam-var" columns are identified and removed from the dataset.
#'          - Non-"bam-var" column names are processed based on \code{sub_pattern} and \code{sub_replacement}.
#'          - The data is then structured into a matrix representation where each row corresponds to a contig,
#'            and columns represent depth values adjusted for each contig's length and relative contribution
#'            to the total average depth.
#'          - A contig's name and its assembly group are concatenated and used as row names for the matrix.
#'          - The function finally returns an OTU table object with taxa as rows.
#'
#' @examples
#' depth_files <- c("depth_file1.txt", "depth_file2.txt")
#' sub_pattern <- "contig"
#' sub_replacement <- "region"
#' otu_table <- import_contig_depths(depth_files, sub_pattern, sub_replacement)
#' 
#' @import phyloseq
#' @import data.table
#' 
#' @export
import_contig_depths <- function(depth_files, sub_pattern, sub_replacement) {
  contig_depths <- rbindlist(lapply(depth_files, function(depth_file) {
    d <- fread(depth_file)
    bam_var <- grep("bam-var", names(d), value = T)
    new_names <- sub(sub_pattern, sub_replacement, grep("bam-var", names(d), value = T, invert = T))
    setNames(d[, !bam_var, with = F], new_names)
  }), use.names = T, id = "AssemblyGroup")
  contig_attributes <- as.matrix(contig_depths[, .(contigLen, totalAvgDepth)])
  contig_matrix <- as.matrix(contig_depths[, lapply(.SD, function(avg_depth) {
    avg_depth * contig_depths$contigLen / sum(avg_depth * contig_depths$contigLen)
  }), .SDcols = !c("contigName", "AssemblyGroup", "contigLen", "totalAvgDepth")])
  # contig_depths$totalAvgDepth * contig_depths
  rownames(contig_matrix) <- rownames(contig_attributes) <- paste(contig_depths$contigName, contig_depths$AssemblyGroup, sep = "@")
  otu_table(contig_matrix, taxa_are_rows = T)
  # Taxonomy table is expected to be characters, so we can't use contig_attributes... too bad.
  # phyloseq(otu_table(contig_matrix, taxa_are_rows = T), tax_table(contig_attributes))
}

#' Import MAG Depths
#'
#' This function takes a file path for a table of microbial abundance by MAG,
#' where the first column is the sample (or taxa) identifier and the remaining columns
#' are numeric depths. It reads the file, converts it into an OTU table format 
#' suitable for downstream ecological analysis.
#'
#' @param depths_file A string path to the file containing the microbial abundance data.
#'
#' @return An object of class `otu_table` with taxa as rows, suitable for ecological analysis.
#'
#' @examples
#' import_mag_depths("path/to/your/depths_file.csv")
#'
#' @importFrom data.table fread
#' @importFrom phyloseq otu_table
#'
#' @export
import_mag_depths <- function(depths_file) {
  mag_depths <- fread(depths_file)
  mag_matrix <- as.matrix(mag_depths[, -1])
  rownames(mag_matrix) <- mag_depths[[1]]
  otu_table(mag_matrix, taxa_are_rows = T)
}
