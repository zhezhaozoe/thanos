#' Order samples in phyloseq object
#'
#' Reorders the samples in a `phyloseq` object based on the specified
#' metadata column.
#'
#' This function reorders the samples in a given `phyloseq` object (or
#' list thereof) according to the values of a specified metadata column
#' in ascending order. This can be useful for plotting or analyses where
#' the order of samples based on metadata is relevant.
#'
#' @param ps A `phyloseq` object (or list thereof) containing sample data.
#' @param by The name of the column in the sample metadata to order the
#' samples by.
#'
#' @return A `phyloseq` object (or list thereof) with reordered samples.
#'
#' @examples
#' # Assuming `ps` is a phyloseq object and `pH` is a metadata column
#' ps_ordered <- order_samples_by(ps, "pH")
#'
#' @export
order_samples_by <- function(ps, by) {
  if (is.list(ps)) {
    sapply(simplify = FALSE, ps, function(psi) {
      levels <- rownames(sample_data(psi))[order(sample_data(psi)[[by]])]
      set_group_order(psi, "Sample", levels)
    })
  } else {
    levels <- rownames(sample_data(ps))[order(sample_data(ps)[[by]])]
    set_group_order(ps, "Sample", levels)
  }
}

#' Set custom order for groups in phyloseq object
#'
#' This function sets a custom ordering for a specified group within the
#' sample data of a phyloseq object (or list thereof). The new ordering
#' is stored in a new column within the sample data, distinguishing
#' between a special case for the "Sample" group and other groups.
#'
#' @param ps A phyloseq object (or list thereof) containing sample data.
#' @param group The name of the group (column in the sample data) to set
#' the order for. If the group is "Sample", the function treats it as a
#' special case.
#' @param levels A vector of values defining the new order of the group.
#'
#' @return A modified phyloseq object (or list thereof) with the new
#' order set for the specified group.
#'
#' @examples
#' # Given a phyloseq object `ps`, set a new order for the "Sample" group
#' ps_new <- set_group_order(ps, "Sample", c("Sample3", "Sample1", "Sample2"))
#'
#' # For a group "Treatment", set a new order
#' ps_new <- set_group_order(ps, "Treatment", c("Control", "Treatment1", "Treatment2"))
#'
#' @export
set_group_order <- function(ps, group, levels) {
  if (is.list(ps)) {
    sapply(simplify = F, ps, function(psi) {
      if (group == "Sample") {
        sample_data(psi)$Sample_zzorder = factor(
          rownames(sample_data(psi)), levels = levels, ordered = T
        )
      } else {
        sample_data(psi)[[paste(group, "zzorder", sep = "_")]] = factor(
          sample_data(psi)[[group]], levels = levels, ordered = T
        )
      }
      psi
    })
  } else {
    if (group == "Sample") {
      sample_data(ps)$Sample_zzorder = factor(
        rownames(sample_data(ps)), levels = levels, ordered = T
      )
    } else {
      sample_data(ps)[[paste(group, "zzorder", sep = "_")]] = factor(
        sample_data(ps)[[group]], levels = levels, ordered = T
      )
    }
    ps
  }
}
