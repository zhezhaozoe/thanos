#' Generate a basic colData file from a mOTUs file
#'
#' @param path Path to a mOTUs file
#' @param level Taxonomic level of the mOTUs in the file
#' Generate Taxonomy Table from mOTUs Data
#'
#' @return A taxonomy table with one column named after the specified
#' taxonomic level (`level`). This column contains the taxonomic names, and
#' the row names of the table are the mOTU identifiers extracted from the
#' data.
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
#' This function takes a user's mOTUs table and a mapping file between
#' mOTUs and GTDB identifiers to create a new taxonomic table with GTDB
#' taxonomy.
#'
#' @param motus_table The user's mOTUs table with rownames as mOTUs identifiers.
#' @param gtdb_map Path to the GTDB mapping file.
#'
#' @return A `tax_table` object containing the GTDB taxonomy for the
#' input mOTUs.
#'
#' @examples
#' # Assuming you've already downloaded the GTDB mapping file:
#' gtdb_map_path <- "path/to/mOTUs_3.0.0_GTDB_tax.tsv"
#' motus_table <- read.csv("your_motus_table.csv", row.names = 1)
#' gtdb_tax_table <- make_motus_tax_table_from_gtdb(motus_table, gtdb_map_path)
#'
#' @export
#'
#' @note Please ensure the GTDB mapping file is available at the
#' specified path. For mOTUs3, the file can be downloaded from
#' https://sunagawalab.ethz.ch/share/MOTU_GTDB/mOTUs_3.0.0_GTDB_tax.tsv.
#' If the file is not found, the function will stop and prompt to
#' download or specify the mapping file.
make_motus_tax_table_from_gtdb <- function(motus_table, gtdb_map) {
  if (!file.exists(gtdb_map)) {
    error("Please download or create a mapping between mOTUs and GTDB identifier. For mOTUs3, the file can be downloaded at https://sunagawalab.ethz.ch/share/MOTU_GTDB/mOTUs_3.0.0_GTDB_tax.tsv")
  }
  gtdb <- as.matrix(fread(gtdb_map, header = FALSE), rownames = 1)
  gtdb_tax_table <- gtdb[rownames(motus_table), ]
  tax_table(gtdb_tax_table)
}

#' Shorten Common Name to Species Initialism
#'
#' This function takes a common species name and shortens it into an
#' initialism format consisting of the first letter of the first word
#' followed by a period and the full second word. Leading and trailing
#' whitespaces are trimmed. Common long words like "Candidatus" and
#' "species incertae sedis" are abbreviated. Useful for reducing the
#' legend clutter when plotting.
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
      common_name <- paste0(
        substr(parts[1], 1, 1), ". ", paste(parts[2], collapse = " ")
      )
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
  all_lines <- fread(file, sep = "", showProgress = FALSE)[[1]]
  content_lines <- grep("^#", all_lines, value = TRUE, invert = TRUE)
  tblout <- setDT(tstrsplit(content_lines, "\\s+")[1:18])
  names(tblout) <- c(
    "Target", "Accession_target", "Query", "Accession_query",
    "Full_evalue", "Full_score", "Full_bias",
    "Best_Evalue", "Best_score", "Best_bias",
    "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc"
  )
  tblout$Description <- gsub("([^ ]+ +){18}", "", content_lines)
  type.convert(tblout, as.is = TRUE, numerals = "no.loss")
}

#' Read GTDB-Tk classification output
#'
#' Reads a GTDB-Tk (Genome Taxonomy Database Toolkit) classification
#' file and extracts the taxonomic classification into a
#' phyloseq-friendly format. It parses the classifications into domain,
#' phylum, class, order, family, genus, and species, and returns a
#' matrix with these categories as column names and genomes as row
#' names.
#'
#' @param file A string specifying the path to the GTDB-Tk
#' classification output file.
#'
#' @return A matrix where rows are genomes (user_genome) and columns are
#' the taxonomic categorizations (Domain, Phylum, Class, Order, Family,
#' Genus, Species).
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
  names(d_classification) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  as.matrix(setDF(d_classification, rownames = d$user_genome))
}

#' Invert Names in a Named Vector
#'
#' Swaps the names and values in a named vector to create a new vector
#' with values as names and names as values.
#'
#' @param v A named vector.
#'
#' @return A named vector with its names and values swapped.
#'
#' @examples
#' my_vector <- c("name1" = "Alice", "name2" = "Bob", "name3" = "Charlie")
#' inverted_names(my_vector)
#' # Returns: Alice = "name1", Bob = "name2", Charlie = "name3"
#'
#' @export
inverted_names <- function(v) {
  with(stack(v), split(as.character(ind), values))
}

#' Merged Hits to OTU Table
#'
#' This function takes a merged hits data.table (with queries and
#' controls), calculates the depths ratio, and returns an OTU
#' (Operational Taxonomic Unit) table, optionally transposing the table
#' based on taxa orientation.
#' It operates on the result of get_hits_depths() when the `phyloseq`
#' argument is FALSE.
#'
#' @param mer A `data.table` with matched pairs of columns for each
#' taxon, each pair consisting of a ".x" and ".y" column, and a column
#' named "ID" that contains the taxa names or IDs.
#' @param taxa_are_rows A logical parameter indicating whether taxa
#' should be represented by rows (TRUE) or by columns (FALSE) in the
#' final OTU table.
#'
#' @details This function cleans the resulting matrix by replacing NaNs and Infs
#' with NAs.
#'
#' @return Returns an `otu_table` object with ratios of hits, with
#' NaNs and Infs set to NA. If `taxa_are_rows` is TRUE, taxa will be
#' represented by rows; otherwise, they will be represented by columns.
#'
#' @import data.table
#' @importFrom phyloseq otu_table
#' @examples
#' # Assuming `dt_hits` is a data.table with matched pairs of hit columns (e.g., "gene1.x", "gene1.y")
#' # and an "ID" column with taxa names.
#' otu_table <- merged_hits_to_otu_table(dt_hits, TRUE)
#'
#' @export
merged_hits_to_otu_table <- function(mer, taxa_are_rows) {
  res <- mer[, .SD, .SDcols = patterns("*\\.x$")] / mer[, .SD, .SDcols = patterns("*\\.y$")]
  setnames(res, names(res), sub("\\.x$", "", names(res)))
  m <- as.matrix(res, rownames.value = mer$ID)
  # XXX: set NaNs (0/0) and Infs (x/0) to NA
  m[is.nan(m)] <- NA
  m[is.infinite(m)] <- NA
  if (isTRUE(taxa_are_rows)) {
    phyloseq::otu_table(m, taxa_are_rows = TRUE)
  } else {
    phyloseq::otu_table(t(m), taxa_are_rows = FALSE)
  }
}

#' Perform HMM Search on a KEGG Orthology (KO) Alignment
#'
#' Given a KO identifier, this function retrieves the corresponding
#' protein sequences and builds a Multiple Sequence Alignment (MSA),
#' then builds a Hidden Markov Model (HMM) based on the MSA, and then
#' searches the HMM against a set of sequence databases. The function
#' allows adjusting the search method, computational resources, and
#' inclusion threshold.
#'
#' @param ko A character string representing the KEGG Orthology (KO)
#' identifier.
#' @param method Character. The method used for performing
#' the MSA. Defaults to "Muscle".
#' @param cpu Integer. The number of CPUs to use for the search.
#' Defaults to 1.
#' @param incE Numeric. The ' inclusion threshold for the HMM search.
#' Defaults to 1e-6.
#' @param ... Additional arguments passed to kegg_kegg_msa().
#'
#' @return The result of the HMM search.
#' @examples
#' tblout_from_ko("K00001")
#' @export
tblout_from_ko <- function(
    ko, dbs, method = "Muscle",
    cpu = 1, incE = 1e-6, hmmer_path = "", ...) {
  aln <- get_kegg_msa(ko, method = method, ...)
  tblout_from_afa(aln)
}

#' Perform HMM-based search from an alignment
#'
#' This function takes aligned fasta sequences (AFA) as an input, constructs
#' a Hidden Markov Model (HMM), and uses it to search in designated databases.
#'
#' @param afa Aligned fasta sequences in string format or an object that
#'  can be interpreted as such.
#' @param cpu The number of CPU cores to use for the search process.
#'  Defaults to 1 to ensure compatibility with all systems.
#' @param incE The inclusion threshold E-value for reported hits in the
#'  search. Lower values are more stringent. Defaults to 1e-6.
#' @param ... Additional arguments to be passed to the search function.
#'
#' @return Returns the path to the tblout results.
#'
#' @examples
#' # Example assuming 'sequences.afa' is your aligned fasta file
#' afa_content <- readLines("sequences.afa")
#' search_results <- tblout_from_afa(afa_content)
#' print(search_results)
#'
#' @export
tblout_from_afa <- function(
    afa, dbs, parallel_processes = 1, cpus_per_process = 1,
    incE = 1e-6, hmmer_path = "") {
  hmm <- build_hmm(afa, hmmer_path = hmmer_path)
  search_hmm(
    hmm,
    dbs,
    parallel_processes = parallel_processes,
    cpus_per_process = cpus_per_process,
    incE = incE,
    hmmer_path = hmmer_path
  )
}

#' Perform a Multiple Sequence Alignment and Find Hits in Databases
#'
#' This function first reads a fasta file containing amino acid sequences, performs
#' a multiple sequence alignment using the specified method, and then searches
#' for hits in databases using the aligned sequences. The results are returned
#' in a specific format.
#'
#' @param faa Character. The path to the fasta file containing amino acid sequences.
#' @param dbs Character or Character vector. Specifies the databases against which the
#' sequences will be searched.
#' @param method Character. The default method for sequence alignment is "Muscle".
#' Other methods supported by `msa` function can be used.
#' @param cpu Integer. The number of CPU threads to use for the search. The default is 1.
#' @param incE Numeric. The inclusion threshold E-value for considering a database hit
#' significant. Defaults to 1e-6.
#' @param ... Additional arguments passed to the `msa` function.
#'
#' @return Returns the result of searching the aligned sequences against specified
#' databases. The format of the returned object depends on the implementation of
#' `tblout_from_afa`.
#'
#' @examples
#' \dontrun{
#' tblout_from_faa("sequences.faa", dbs=c("db1", "db2"), method="ClustalW")
#' }
#'
#' @export
#' @importFrom Biostrings readAAStringSet
#' @import msa
tblout_from_faa <- function(
    faa, dbs, method = "Muscle",
    cpu = 1, incE = 1e-6, ...) {
  my_sequences <- Biostrings::readAAStringSet(faa)
  aln <- msa::msa(my_sequences, method = method, ...)
  tblout_from_afa(aln, cpu = cpu, incE = incE)
}
