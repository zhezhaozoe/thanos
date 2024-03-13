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
  p <- ggplot(d, aes(x = Sample, y = Abundance, fill = if (is.null(fill)) NULL else .data[[fill]])) +
    geom_col(position = position, ...) +
    labs(fill = fill) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5))
  if (all(d[, sum(Abundance), by = Sample]$V1 == 1, na.rm = T) || position == "fill") {
    p <- p + scale_y_continuous(labels = scales::percent)
  }
  p
}

# barplot_depth_by_group(ps_list, group = "group", fill = "Gene")
# barplot_depth_by_group(ps_list, group = "group", fill = "Phylum")
# barplot_depth_by_group(ps_list, group = "Gene", fill = "Phylum")
# barplot_depth_by_group(ps_list, group = c("Gene", "group"), fill = "Phylum") + facet_wrap(~group, scales = "free")
# The last one doesn't really work, sadly... maybe adding a wrap argument is a better interface.
barplot_depth_by_group <- function(ps, group = "Sample", fill = NULL, position = "stack", ...) {
  d <- if (is.list(ps)) {
    rbindlist(lapply(ps_list, function(ps) {
      setDT(psmelt(ps))
    }), id = "Gene")
  } else {
    setDT(psmelt(ps))
  }
  d <- d[, .(Abundance = sum(Abundance)), by = c(group, fill)]
  p <- ggplot(d, aes(x = .data[[group]], y = Abundance, fill = if (is.null(fill)) NULL else .data[[fill]])) +
    geom_col(position = position, ...) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) +
    labs(fill = fill, x = group) +
  # The c() in by is necessary, otherwise group will be interpreted as a column of d (if present)
  if (all(d[, sum(Abundance), by = c(group)]$V1 == 1, na.rm = T) || position == "fill") {
    p <- p + scale_y_continuous(labels = scales::percent)
  }
  p
}

boxplot_depth_by_group <- function(ps, group, x = group, ...) {
  d <- if (is.list(ps)) {
    rbindlist(lapply(ps_list, function(ps) {
      setDT(psmelt(ps))
    }), id = "Gene")
  } else {
    setDT(psmelt(ps))
  }
  # The c() in by is necessary, otherwise group will be interpreted as a column of d (if present)
  d <- d[, .(Abundance = sum(Abundance)), by = c(unique(c(group, x)))]
  ggplot(d, aes(x = factor(.data[[x]]), y = Abundance)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) +
    xlab(group)
}
