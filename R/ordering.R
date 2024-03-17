order_samples_by <- function(ps, by) {
  levels <- rownames(sample_data(ps))[order(sample_data(ps)[[by]])]
  set_group_order(ps, "Sample", levels)
}

set_group_order <- function(ps, group, levels) {
  if (group == "Sample") {
    sample_data(ps)$Sample_zzorder = factor(rownames(sample_data(ps)), levels = levels, ordered = T)
  } else {
    sample_data(ps)[[paste(group, "zzorder", sep = "_")]] = factor(sample_data(ps)[[group]], levels = levels, ordered = T)
  }
  ps
}
