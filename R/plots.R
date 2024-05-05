#' Plots aggregated abundance by sample in barplot
#'
#' This function plots the aggregated abundance of genes by sample and
#' allows for customization of the filling color of the bars. Unlike
#' `barplot_depths()`, this function retains the sample data, enabling
#' the use of `facet_wrap(~whatever)` for further customization.
#'
#' @param ps A phyloseq object or a list of such objects. If a list is
#' provided, it merges the data from all objects.
#' @param fill The variable by which to fill the bars. If `NULL`, bars
#' won't be filled according to a variable. The default is `NULL`.
#' @param position The position adjustment of bars. Can be "stack",
#' "fill", or others as defined by `ggplot2`. Default is "stack".
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
  d <- d[,
    .(Abundance = mean(Abundance, na.rm = TRUE)),
    by = c("Sample", names(sample_data(ps)), fill)
  ]
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
  if (all(d[, sum(Abundance, na.rm = TRUE), by = Sample]$V1 == 1) || position == "fill") {
    p <- p + scale_y_continuous(labels = scales::percent)
  }
  p
}

#' Create a Bar Plot of Abundances from Phyloseq Data with Grouping and Fill Options
#'
#' @param ps A phyloseq object or a list of phyloseq objects. If a list
#' is provided, all objects are combined.
#' @param group The name of the column to be used for grouping on the
#' x-axis. Defaults to "Sample".
#' @param fill (Optional) The name of the column to apply a fill color
#' to bars. NULL by default, which results in no fill.
#' @param position Determines how to position bars. "stack" stacks bars
#' on top of each other, "dodge" places them side by side. Defaults to
#' "stack".
#' @param wrap (Optional) The name of a column to be used for facetting
#' the plot into separate panels. It can also be a vector of two names,
#' in which case it will use facet_grid() with the first variable on the
#' rows and the second one on the columns.
#' @param ... Additional arguments passed on to `geom_col` from
#'`ggplot2`.
#'
#' @return A `ggplot` object representing the bar plot. Can be printed
#' or modified further.
#
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
#' @importFrom phyloseq psmelt
#'
#' @export
barplot_depths <- function(
    ps, group = "Sample", fill = NULL,
    position = "stack", wrap = NULL, ...) {
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
  d <- d[,
    .(Abundance = mean(Abundance, na.rm = TRUE)),
    by = c(by_cols, ordered_by_cols)
  ]
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
    labs(fill = fill, x = group)
  if (length(wrap) == 1) {
    p <- p + facet_wrap(~ .data[[wrap]], scale = "free_x")
  } else if (length(wrap) == 2) {
    p <- p + facet_grid(
      formula(paste(wrap, collapse = "~")),
      scales = "free_x",
      space = "free_x"
    )
  }
  # The c() in by is necessary, otherwise group will be interpreted as a column of d (if present)
  if (all(d[, sum(Abundance, na.rm = TRUE), by = c(group)]$V1 == 1) || position == "fill") {
    p <- p + scale_y_continuous(labels = scales::percent)
  }
  p
}

#' Generate a Boxplot with Significance Indicators

#' @param ps A phyloseq object or a list of phyloseq objects. The
#' function melts the data for plotting.
#' @param fill The name of the variable to use for filling boxes. Can be
#' set to NULL. Defaults to NULL.
#' @param x The name of the variable to use for the x-axis. Inherits the
#' value of `fill` by default.
#' @param wrap An optional string specifying a variable name to wrap the
#' plot into multiple panels. Not supported alongside `signif`.
#' @param signif Logical, indicating whether to perform and display
#' significance testing between groups. Defaults to TRUE.
#' @param test The name of the statistical test function to use as a
#' character string. Defaults to "wilcox.test".
#' @param ... Additional arguments passed to `geom_boxplot`.
#'
#' @details This function creates a boxplot of abundance data from a
#' given phyloseq object or a list of phyloseq objects. When `signif` is
#' TRUE, it automatically calculates and displays significance levels
#' between the comparisons of boxplots based on the specified statistical
#' test. `wrap` parameter allows creating multiple panels but is not
#' compatible with `signif` parameter. If using `signif`, consider using
#' `patchwork` for combining multiple plots.
#'
#' @return A ggplot object representing the boxplot with or without
#' significance indicators.
#'
#' @note Using `wrap` together with `signif` will result in an error
#' message as they are not compatible. Use `patchwork` package for layout
#' management in such cases.
#'
#' @examples
#' # Assuming `ps` is a phyloseq object with proper setup:
#' boxplot_depths(ps, fill = "SampleType", x = "Condition", signif = TRUE)
#' # For multiple phyloseq objects in a list and using a Wilcoxon test:
#' boxplot_depths(list(ps1, ps2), fill = "SampleType", x = "Condition", test = "wilcox.test", signif = TRUE)
#'
#' @import ggplot2
#' @import data.table
#' @importFrom phyloseq psmelt
#' @export
boxplot_depths <- function(
    ps, fill = NULL, x = fill, wrap = NULL,
    signif = TRUE, test = "wilcox.test", ...) {
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
    pairs <- CJ(V1 = x_positions, V2 = x_positions)[V1 < V2]
    comps <- rbindlist(apply(pairs, 1, function(r) {
      test_res <- match.fun(test)(
        as.formula(paste("Abundance", x, sep = " ~ ")), d[as.integer(d[[x]]) %in% r]
      )
      data.table(x1 = r[1], x2 = r[2], p.value = test_res$p.value)
    }))
    lev <- 1
    for (span in seq_len(length(x_positions) - 1)) {
      for (start in seq_len(min(span, length(x_positions) - span))) {
        mask <- data.table(
          x1 = seq(start, length(x_positions) - span, span),
          x2 = seq(start + span, length(x_positions), span)
        )
        comps[mask, level := lev, on = c("x1", "x2")]
        lev <- lev + 1
      }
    }
    comps$max_y <- max(d$Abundance, na.rm = TRUE)
    comps$shift_y <- (max(d$Abundance, na.rm = TRUE) - min(d$Abundance, na.rm = TRUE)) / 20
    p <- p +
      geom_segment(
        data = comps,
        aes(x = x1, xend = x2, y = max_y + 2 * shift_y * level, yend = max_y + 2 * shift_y * level),
        inherit.aes = FALSE
      ) +
      geom_segment(
        data = comps,
        aes(x = x1, xend = x1, y = max_y + 2 * shift_y * level, yend = max_y + 2 * shift_y * level - shift_y),
        inherit.aes = FALSE
      ) +
      geom_segment(
        data = comps,
        aes(x = x2, xend = x2, y = max_y + 2 * shift_y * level, yend = max_y + 2 * shift_y * level - shift_y),
        inherit.aes = FALSE
      ) +
      geom_text(
        data = comps,
        aes(x = (x1 + x2) / 2, y = max_y + 2 * shift_y * level, label = prettyNum(signif(`p.value`, 2))),
        vjust = 0, nudge_y = comps$shift_y / 2,
        inherit.aes = FALSE
      )
  }
  if (!is.null(wrap)) {
    p <- p + facet_wrap(~ .data[[wrap]], scale = "free_x")
  }
  p
}

#' Plot Enzyme Depths for a Given KEGG Module and Depths List
#'
#' @param kegg_module The module name as a string (e.g. "M00001") or
#' a data.table containing the KEGG pathway graph data (derived from
#' get_kegg_reactions_from_module()).
#' @param depths_list A list of phyloseq objects representing the depths.
#' @param wrap An optional column name as string for faceting.
#' If `NULL`, no wrapping is performed.
#' @param plot If FALSE, do not plot anything, just return the raw table.
#'
#' @return A ggraph plot of the specified KEGG module, annotated with
#' enzyme depths. if `plot` is FALSE, a tidygraph object with the depths
#' information will be returned instead of the plot.
#'
#' @description This function takes a KEGG pathway graph and a list of
#' depths data, aggregates the depths information, and creates an igraph
#' object representing the KEGG pathway with aggregated node abundances.
#'
#' @import data.table
#' @importFrom igraph graph_from_data_frame
#' @import ggraph
#' @importFrom tidygraph as_tbl_graph filter
#' @importFrom phyloseq psmelt
#'
#' @export
keggmodule_plot <- function(
    kegg_module, depths_list, wrap = NULL, plot = TRUE) {
  if (is.character(kegg_module)) {
    # Assume it's the module name, else assume it's the graph already
    kegg_module <- get_kegg_reactions_from_module(kegg_module)
  }
  depths_list_aggr <- unique(rbindlist(lapply(depths_list, function(psi) {
    setDT(psmelt(psi))
  }), id = "zzthanos_KO"))
  extract_kos <- function(reaction_text) {
    unique(grep("^K", strsplit(reaction_text, "[:,|+]")[[1]], value = TRUE))
  }
  reactions_to_kos <- unique(
    kegg_module[, .(zzthanos_KO = extract_kos(.BY[[1]])), by = reaction]
  )
  kegg_dt <- merge(
    merge(kegg_module, reactions_to_kos, by = "reaction"),
    depths_list_aggr, by = "zzthanos_KO", all.x = TRUE, allow.cartesian = TRUE
  )
  kegg_dt_aggr <- kegg_dt[,
    .(Abundance = mean(Abundance, na.rm = TRUE)),
    by = c("from", "to", "from_name", "to_name", "reaction", wrap)
  ]
  if (!is.null(wrap)) {
    for (var in c("from", "reaction", "to")) {
      kegg_dt_aggr[[var]] <- paste(
        kegg_dt_aggr[[var]], kegg_dt_aggr[[wrap]], sep = "@"
      )
    }
  }
  kegg_edges <- rbind(
    data.table(
      from = kegg_dt_aggr$from,
      to = kegg_dt_aggr$reaction,
      type = "blunt"
    ),
    data.table(
      from = kegg_dt_aggr$reaction,
      to = kegg_dt_aggr$to,
      type = "arrow"
    )
  )
  kegg_nodes <- unique(rbind(
    data.table(
      name = kegg_dt_aggr$from,
      common_name = kegg_dt_aggr$from_name,
      type = "compound",
      value = 0,
      wrap = if (is.null(wrap)) 0 else kegg_dt_aggr[[wrap]]
    ),
    data.table(
      name = kegg_dt_aggr$to,
      common_name = kegg_dt_aggr$to_name,
      type = "compound",
      value = 0,
      wrap = if (is.null(wrap)) 0 else kegg_dt_aggr[[wrap]]
    ),
    data.table(
      name = kegg_dt_aggr$reaction,
      common_name = "meow",
      type = "reaction",
      value = kegg_dt_aggr$Abundance,
      wrap = if (is.null(wrap)) 0 else kegg_dt_aggr[[wrap]]
    )
  ))
  graph <- igraph::graph_from_data_frame(
    kegg_edges, directed = TRUE, vertices = kegg_nodes
  )
  if (isFALSE(plot)) {
    return(tidygraph::as_tbl_graph(graph))
  }
  ggraph(graph, "sugiyama") +
    geom_edge_link2(
      data = function(x) tidygraph::filter(get_edges("long")(x), type == "blunt"),
      arrow = NULL) +
    geom_edge_link2(
      data = function(x) tidygraph::filter(get_edges("long")(x), type == "arrow"),
      angle_calc = "along",
      arrow = arrow(length = unit(2, "mm")),
      end_cap = label_rect(
        "C00000",
        cex = GeomLabel$default_aes$size / 11 * .pt) # 11/.pt is the default size
      ) +
    geom_node_label(
      data = function(x) tidygraph::filter(get_nodes()(x), type == "compound"),
      aes(label = name)
    ) +
    geom_node_label(
      data = function(x) tidygraph::filter(get_nodes()(x), type == "reaction"),
      aes(
        label = name,
        fill = value,
        color = ifelse((value - min(value)) / (max(value) - min(value)) > .5, "black", "white")
      )
    ) +
    labs(fill = "Abundance") +
    scale_fill_viridis_c() +
    scale_color_identity() +
    NULL
}
