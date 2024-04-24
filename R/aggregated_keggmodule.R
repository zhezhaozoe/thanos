aggregate_module_depths <- function(kegg_graph, depths_list, samples = NULL, taxa = NULL) {
  depths_list_aggr <- lapply(depths_list, function(psi) {
    d <- setDT(psmelt(psi))
    if (!is.null(samples)) {
      d <- d[Sample %in% samples]
    }
    if (!is.null(taxa)) {
      # TODO: maybe use NSE to pass whatever to phyloseq::subset_taxa()
      d <- d[Phylum %in% taxa] # something like this, but why phylum?
      print("meow")
    }
    sum(d$Abundance)
  })
  kegg_graph$kos_depths <- sapply(strsplit(kegg_graph$reaction, ":|,|\\|"), function(reaction) {
    reaction_kos <- grep("^K", reaction, value = T)
    sum(unlist(depths_list_aggr[reaction_kos]))
  })
  kegg_edges <- rbind(
    data.table(
      from = kegg_graph$from,
      to = kegg_graph$reaction,
      type = "blunt"
    ),
    data.table(
      from = kegg_graph$reaction,
      to = kegg_graph$to,
      type = "arrow"
    )
  )
  kegg_nodes <- unique(rbind(
    data.table(
      name = kegg_graph$from,
      common_name = kegg_graph$from_name,
      type = "compound",
      value = 0
    ),
    data.table(
      name = kegg_graph$to,
      common_name = kegg_graph$to_name,
      type = "compound",
      value = 0
    ),
    data.table(
      name = kegg_graph$reaction,
      common_name = "meow",
      type = "reaction",
      value = kegg_graph$kos_depths
    )
  ))
  graph <- igraph::graph_from_data_frame(kegg_edges, directed = T, vertices = kegg_nodes)
}

#' @import ggraph
plot_aggregated_hits <- function(graph) {
  tbl <- as_tbl_graph(graph)
  ggraph(tbl, "sugiyama") +
    geom_edge_link2(
      data = \(x) {get_edges("long")(x) |> filter(type == "blunt")},
      arrow = NULL) +
    geom_edge_link2(
      data = \(x) {get_edges("long")(x) |> filter(type == "arrow")},
      angle_calc = "along",
      arrow = arrow(length = unit(4, "mm")), 
      end_cap = label_rect("C00000")) +
    geom_node_label(
      data = \(x) {get_nodes()(x) |> filter(type == "compound")},
      aes(label = name)
    ) +
    geom_node_label(
      data = \(x) {get_nodes()(x) |> filter(type == "reaction")},
      aes(label = name, fill = value, color = ifelse((value - min(value)) / (max(value) - min(value)) > .5, "black", "white"))
    ) +
    scale_fill_viridis_c() +
    scale_color_identity() +
    NULL
}
