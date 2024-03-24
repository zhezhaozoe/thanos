#' Plots aggregated abundance by sample in barplot
#'
#' This function plots the aggregated abundance of genes by sample and allows for
#' customization of the filling color of the bars. Unlike `barplot_depths()`,
#' this function retains the sample data, enabling the use of `facet_wrap(~whatever)`
#' for further customization.
#'
#' @param ps A phyloseq object or a list of such objects. If a list is provided, it
#'   merges the data from all objects.
#' @param fill The variable by which to fill the bars. If `NULL`, bars won't be
#'   filled according to a variable. The default is `NULL`.
#' @param position The position adjustment of bars. Can be "stack", "fill", or
#'   others as defined by `ggplot2`. Default is "stack".
#' @param ... Additional arguments passed to `geom_col`.
#'
#' @return A `ggplot` object representing the bar plot.
#'
#' @examples
#' ## Assuming 'ps' is a phyloseq object and 'ps_list' is a list of phyloseq objects
#' barplot_depths_by_sample(ps, fill = "Kingdom")
#' barplot_depths_by_sample(ps_list, fill = "Phylum", position = "fill")
#'
#' @export
barplot_depths_by_sample <- function(ps, fill = NULL, position = "stack", ...) {
  d <- if (is.list(ps)) {
    rbindlist(lapply(ps, function(psi) {
      setDT(psmelt(psi))
    }), id = "Gene")
  } else {
    setDT(psmelt(ps))
  }
  d <- d[, .(Abundance = sum(Abundance)), by = c("Sample", names(sample_data(ps)), fill)]
  p <- ggplot(
    d,
    aes(
      x = if ("Sample_zzorder" %in% names(sample_data(ps))) Sample_zzorder else Sample,
      y = Abundance,
      fill = if (is.null(fill)) NULL else if (paste(fill, "zzorder", sep = "_") %in% names(d)) .data[[paste(fill, "zzorder", sep = "_")]] else .data[[fill]]
    )
  ) +
    geom_col(position = position, ...) +
    labs(fill = fill, x = "Sample") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5))
  if (all(d[, sum(Abundance), by = Sample]$V1 == 1, na.rm = T) || position == "fill") {
    p <- p + scale_y_continuous(labels = scales::percent)
  }
  p
}

#' Create a Bar Plot of Abundances from Phyloseq Data with Grouping and Fill Options
#'
#' @param ps A phyloseq object or a list of phyloseq objects. If a list is provided, all objects are combined.
#' @param group The name of the column to be used for grouping on the x-axis. Defaults to "Sample".
#' @param fill (Optional) The name of the column to apply a fill color to bars. NULL by default, which results in no fill.
#' @param position Determines how to position bars. "stack" stacks bars on top of each other, "dodge" places them side by side. Defaults to "stack".
#' @param wrap (Optional) The name of a column to be used for facetting the plot into separate panels.
#' @param ... Additional arguments passed on to `geom_col` from `ggplot2`.
#'
#' @return A `ggplot` object representing the bar plot. Can be printed or modified further.
#'
#' @examples
#' data("GlobalPatterns", package = "phyloseq")
#' barplot_depths(GlobalPatterns, group = "SampleType", fill = "Phylum")
#' barplot_depth(ps_list, group = "group", fill = "Gene")
#' barplot_depth(ps_list, group = "group", fill = "Phylum")
#' barplot_depth(ps_list, group = "Gene", fill = "Phylum")
#' barplot_depth(ps_list, group = c("Gene", "group"), fill = "Phylum") + facet_wrap(~group, scales = "free")
#'
#'
#' @import ggplot2
#' @import data.table
barplot_depths <- function(ps, group = "Sample", fill = NULL, position = "stack", wrap = NULL, ...) {
  d <- if (is.list(ps)) {
    rbindlist(lapply(ps, function(psi) {
      setDT(psmelt(psi))
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
  p <- ggplot(
    d,
    aes(
      x = if (paste(group, "zzorder", sep = "_") %in% names(d)) {
        .data[[paste(group, "zzorder", sep = "_")]]
      } else {
        .data[[group]]
      },
      fill = if (is.null(fill)) {
        NULL
      } else if (paste(fill, "zzorder", sep = "_") %in% names(d)) {
        .data[[paste(fill, "zzorder", sep = "_")]]
      } else {
        .data[[fill]]
      },
      y = Abundance
    )
  ) +
    geom_col(position = position, ...) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) +
    labs(fill = fill, x = group) +
    if (is.null(wrap)) NULL else facet_wrap(~ .data[[wrap]], scale = "free_x")
  # The c() in by is necessary, otherwise group will be interpreted as a column of d (if present)
  if (all(d[, sum(Abundance), by = c(group)]$V1 == 1, na.rm = T) || position == "fill") {
    p <- p + scale_y_continuous(labels = scales::percent)
  }
  p
}

#' Generate a Boxplot with Significance Indicators
#'
#' @param ps A phyloseq object or a list of phyloseq objects. The function melts the data for plotting.
#' @param fill The name of the variable to use for filling boxes. Can be set to NULL. Defaults to NULL.
#' @param x The name of the variable to use for the x-axis. Inherits the value of `fill` by default.
#' @param wrap An optional string specifying a variable name to wrap the plot into multiple panels. Not supported alongside `signif`.
#' @param signif Logical, indicating whether to perform and display significance testing between groups. Defaults to TRUE.
#' @param test The name of the statistical test function to use as a character string. Defaults to "wilcox.test".
#' @param ... Additional arguments passed to `geom_boxplot`.
#'
#' @details This function creates a boxplot of abundance data from a given phyloseq object or a list of phyloseq objects. When `signif` is TRUE, it automatically calculates and displays significance levels between the comparisons of boxplots based on the specified statistical test. `wrap` parameter allows creating multiple panels but is not compatible with `signif` parameter. If using `signif`, consider using `patchwork` for combining multiple plots.
#'
#' @return A ggplot object representing the boxplot with or without significance indicators.
#'
#' @note Using `wrap` together with `signif` will result in an error message as they are not compatible. Use `patchwork` package for layout management in such cases.
#'
#' @examples
#' # Assuming `ps` is a phyloseq object with proper setup:
#' boxplot_depths(ps, fill = "SampleType", x = "Condition", signif = TRUE)
#' # For multiple phyloseq objects in a list and using a Wilcoxon test:
#' boxplot_depths(list(ps1, ps2), fill = "SampleType", x = "Condition", test = "wilcox.test", signif = TRUE)
#'
#' @import ggplot2
#' @import data.table
#' @export
boxplot_depths <- function(ps, fill = NULL, x = fill, wrap = NULL, signif = TRUE, test = "wilcox.test", ...) {
  if (!is.null(wrap) && isTRUE(signif)) {
    stop("`wrap` is not supported with `signif`, please use patchwork instead.")
  }
  d <- if (is.list(ps)) {
    rbindlist(lapply(ps, function(psi) {
      setDT(psmelt(psi))
    }), id = "Gene")
  } else {
    setDT(psmelt(ps))
  }
  if (!is.factor(d[[x]])) {
    d[[x]] <- factor(d[[x]])
  }
  p <- ggplot(d, aes(
    x = if (is.null(x)) {
      NULL
    } else if (paste(x, "zzorder", sep = "_") %in% names(d)) {
      .data[[paste(x, "zzorder", sep = "_")]]
    } else {
      .data[[x]]
    },
    fill = if (is.null(fill)) {
      NULL
    } else if (paste(fill, "zzorder", sep = "_") %in% names(d)) {
      .data[[paste(fill, "zzorder", sep = "_")]]
    } else {
      .data[[fill]]
    },
    y = Abundance
  )) +
    geom_boxplot(...) +
    labs(x = x, fill = fill)
  if (isTRUE(signif)) {
    x_positions <- as.integer(unique(d[[x]]))
    comps <- rbindlist(apply(CJ(V1 = x_positions, V2 = x_positions)[V1 < V2], 1, function(r) {
      test_res <- match.fun(test)(as.formula(paste("Abundance", x, sep = " ~ ")), d[as.integer(d[[x]]) %in% r])
      data.table(x1 = r[1], x2 = r[2], p.value = test_res$p.value)
    }))
    lev <- 1
    for (span in 1:(length(x_positions) - 1)) {
      for (start in 1:(min(span, length(x_positions) - span))) {
        mask <- data.table(
          x1 = seq(start, length(x_positions) - span, span),
          x2 = seq(start + span, length(x_positions), span)
        )
        comps[mask, level := lev, on = c("x1", "x2")]
        lev <- lev + 1
      }
    }
    comps$max_y <- max(d$Abundance)
    comps$shift_y <- (max(d$Abundance) - min(d$Abundance)) / 20
    p <- p +
      geom_segment(data = comps, aes(x = x1, xend = x2, y = max_y + 2 * shift_y * level, yend = max_y + 2 * shift_y * level), inherit.aes = F) +
      geom_segment(data = comps, aes(x = x1, xend = x1, y = max_y + 2 * shift_y * level, yend = max_y + 2 * shift_y * level - shift_y), inherit.aes = F) +
      geom_segment(data = comps, aes(x = x2, xend = x2, y = max_y + 2 * shift_y * level, yend = max_y + 2 * shift_y * level - shift_y), inherit.aes = F) +
      geom_text(data = comps, aes(x = (x1 + x2) / 2, y = max_y + 2 * shift_y * level, label = prettyNum(`p.value`, 2)), vjust = 0, nudge_y = comps$shift_y / 2, inherit.aes = F)
  }
  if (!is.null(wrap)) {
    p <- p + facet_wrap(~ .data[[wrap]], scale = "free_x")
  }
  p
}
