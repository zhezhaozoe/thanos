#' @export
import_motus <- function(path, sample_data_d = NULL, tax_table_m = NULL) {
  otu <- otu_table(read_motus(path), taxa_are_rows = T)
  tax <- if (is.null(tax_table_m)) {
    tax_table(make_tax_table_from_motus(path))
  } else {
    tax_table_m
  }
  phyloseq(otu, sample_data_d, tax)
}

#' @export
#' @import data.table
import_contig_depths <- function(depth_files, sample_data_d = NULL, seq_files = NULL, sub_pattern = "", sub_replacement = "") {
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
  ps <- phyloseq(otu_table(contig_matrix, taxa_are_rows = T), sample_data_d)
  attr(ps, "contigs") <- contig_attributes
  if (!is.null(seq_files)) {
    if (!isTRUE(all.equal(sort(names(depth_files)), sort(names(seq_files))))) {
      stop("The names of the depths_files and the names of seq_files must be identical.")
    }
    attr(ps, "seq_files") <- seq_files[names(depth_files)]
  }
  ps
}

#' @export
import_mag_depths <- function(depths_file, sample_data_d = NULL, tax_table_m = NULL, seq_files = NULL) {
  mag_depths <- fread(depths_file)
  mag_matrix <- as.matrix(mag_depths[, -1])
  rownames(mag_matrix) <- mag_depths[[1]]
  ps <- phyloseq(otu_table(mag_matrix, taxa_are_rows = T), sample_data_d, tax_table_m)
  if (!is.null(seq_files)) {
    if (!isTRUE(all.equal(sort(rownames(mag_matrix)), sort(names(seq_files))))) {
      stop("The rownames of depths_file and the names of seq_files must be identical.")
    }
    attr(ps, "seq_files") <- seq_files[rownames(mag_matrix)]
  }
  ps
}
