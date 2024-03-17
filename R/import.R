#' @description Keep only the motus ID, not the tax name, which should be provided separately as colData.
#'
#' @import data.table
#' @export
import_motus <- function(path, make_tax_table = F) {
  d <- fread(path)
  names(d)[1] <- "ID"
  d$ID <- gsub(".*\\[(.*)\\]", "\\1", d$ID)
  m <- as.matrix(d[, !"ID"])
  rownames(m) <- d$ID
  ps <- phyloseq(otu_table(read_motus(path), taxa_are_rows = T))
  if (make_tax_table) {
    tax_table(ps) <- tax_table(make_tax_table_from_motus(path))
  }
}

#' @export
#' @import data.table
import_contig_depths <- function(depth_files) {
  contig_depths <- rbindlist(lapply(depth_files, function(depth_file) {
    d <- fread(depth_file)
    bam_var <- grep("bam-var", names(d), value = T)
    new_names <- sub(sub_pattern, sub_replacement, grep("bam-var", names(d), value = T, invert = T))
    setNames(d[, !bam_var, with = F], new_names)
  }), use.names = T, id = "AssemblyGroup")
  contig_attributes <- contig_depths[, .(contigName, AssemblyGroup, contigLen, totalAvgDepth)]
  contig_depths[, c("contigName", "AssemblyGroup") := .(paste(contigName, AssemblyGroup, sep = "@"), NULL)]
  contig_matrix <- as.matrix(contig_depths[, !c("contigName", "contigLen", "totalAvgDepth")])
  rownames(contig_matrix) <- contig_depths$contigName
  otu_table(contig_matrix, taxa_are_rows = T)
}

#' @export
#' @import data.table
import_mag_depths <- function(depths_file) {
  mag_depths <- fread(depths_file)
  mag_matrix <- as.matrix(mag_depths[, -1])
  rownames(mag_matrix) <- mag_depths[[1]]
  otu_table(mag_matrix, taxa_are_rows = T)
}
