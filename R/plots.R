# fill by gene for plot 5.
# this function (as opposed to *_by_sample_group) keeps the sample_data, so it's possible to do + facet_wrap(~whatever) afterwards.
barplot_depth_by_sample <- function(ps, fill = NULL, position = "stack", ...) {
  d <- if (is.list(ps)) {
    rbindlist(lapply(ps_list, function(ps) {
      setDT(psmelt(ps))
    }), id = "Gene")
  } else {
    setDT(psmelt(ps))
  }
  d <- d[, .(Abundance = sum(Abundance)), by = c("Sample", names(sample_data(ps)), fill)]
  p <- ggplot(d,
    aes(
      x = if ("Sample_zzorder" %in% names(sample_data(ps))) Sample_zzorder else Sample,
      y = Abundance,
      fill = if (is.null(fill)) NULL else if (paste(fill, "zzorder", sep = "_") %in% names(d)) .data[[paste(fill, "zzorder", sep = "_")]] else .data[[fill]])) +
    geom_col(position = position, ...) +
    labs(fill = fill, x = "Sample") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5))
  if (all(d[, sum(Abundance), by = Sample]$V1 == 1, na.rm = T) || position == "fill") {
    p <- p + scale_y_continuous(labels = scales::percent)
  }
  p
}

# barplot_depth(ps_list, group = "group", fill = "Gene")
# barplot_depth(ps_list, group = "group", fill = "Phylum")
# barplot_depth(ps_list, group = "Gene", fill = "Phylum")
# barplot_depth(ps_list, group = c("Gene", "group"), fill = "Phylum") + facet_wrap(~group, scales = "free")
barplot_depth <- function(ps, group = "Sample", fill = NULL, position = "stack", wrap = NULL, ...) {
  d <- if (is.list(ps)) {
    rbindlist(lapply(ps_list, function(ps) {
      setDT(psmelt(ps))
    }), id = "Gene")
  } else {
    setDT(psmelt(ps))
  }
  by_cols <- c(group, fill, wrap)
  ordered_by_cols <- c()
  for (col in by_cols) {
    ordered_col <- paste(col, "zzorder", sep = "_")
    if (ordered_col %in% names(d)) {
      ordered_by_cols <- c(ordered_by_cols, ordered_col)
    }
  }
  d <- d[, .(Abundance = sum(Abundance)), by = c(by_cols, ordered_by_cols)]
  p <- ggplot(d,
    aes(
      x = if (paste(group, "zzorder", sep = "_") %in% names(d)) { .data[[paste(group, "zzorder", sep = "_")]] } else { .data[[group]] },
      fill = if (is.null(fill)) { NULL } else if (paste(fill, "zzorder", sep = "_") %in% names(d)) { .data[[paste(fill, "zzorder", sep = "_")]] } else { .data[[fill]] },
      y = Abundance
  )) +
    geom_col(position = position, ...) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) +
    labs(fill = fill, x = group) +
    if (is.null(wrap)) NULL else facet_wrap(~.data[[wrap]], scale = "free_x")
  # The c() in by is necessary, otherwise group will be interpreted as a column of d (if present)
  if (all(d[, sum(Abundance), by = c(group)]$V1 == 1, na.rm = T) || position == "fill") {
    p <- p + scale_y_continuous(labels = scales::percent)
  }
  p
}

boxplot_depth <- function(ps, fill = NULL, x = fill, wrap = NULL, ...) {
  d <- if (is.list(ps)) {
    rbindlist(lapply(ps_list, function(ps) {
      setDT(psmelt(ps))
    }), id = "Gene")
  } else {
    setDT(psmelt(ps))
  }
  ggplot(d,
    aes(
      x = if (is.null(x)) { NULL } else if (paste(x, "zzorder", sep = "_") %in% names(d)) { .data[[paste(x, "zzorder", sep = "_")]] } else { .data[[x]] },
      fill = if (is.null(fill)) { NULL } else if (paste(fill, "zzorder", sep = "_") %in% names(d)) { .data[[paste(fill, "zzorder", sep = "_")]] } else { .data[[fill]] },
      y = Abundance
  )) +
    geom_boxplot(...) +
    labs(x = x, fill = fill) +
    if (is.null(wrap)) NULL else facet_wrap(~.data[[wrap]])
}
