aggregate_module_depths <- function(kegg_graph, depths_list, wrap = NULL) {
  depths_list_aggr <- unique(rbindlist(lapply(depths_list, function(psi) {
    setDT(psmelt(psi))
  }), id = "zzthanos_KO"))
  reactions_to_kos <- unique(kegg_graph[, .(zzthanos_KO = grep("^K", strsplit(.BY[[1]], ":|,|\\|")[[1]], value = T)), by = reaction])
  kegg_dt <- merge(
    merge(kegg_graph, reactions_to_kos, by = "reaction"),
    depths_list_aggr, by = "zzthanos_KO", all.x = T, allow.cartesian = TRUE)
  kegg_dt_aggr <- kegg_dt[, .(Abundance = sum(Abundance)), by = c("from", "to", "from_name", "to_name", "reaction", wrap)]
  if (!is.null(wrap)) {
    for (var in c("from", "reaction", "to")) {
      kegg_dt_aggr[[var]] <- paste(kegg_dt_aggr[[var]], kegg_dt_aggr[[wrap]], sep = "@")
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
  graph <- igraph::graph_from_data_frame(kegg_edges, directed = T, vertices = kegg_nodes)
}

#' @import ggraph
#' @importFrom tidygraph as_tbl_graph filter
#' @export
graphplot_depths <- function(graph) {
  tbl <- tidygraph::as_tbl_graph(graph)
  ggraph(tbl, "sugiyama") +
    geom_edge_link2(
      data = \(x) {get_edges("long")(x) |> tidygraph::filter(type == "blunt")},
      arrow = NULL) +
    geom_edge_link2(
      data = \(x) {get_edges("long")(x) |> tidygraph::filter(type == "arrow")},
      angle_calc = "along",
      arrow = arrow(length = unit(4, "mm")), 
      end_cap = label_rect("C00000")) +
    geom_node_label(
      data = \(x) {get_nodes()(x) |> tidygraph::filter(type == "compound")},
      aes(label = name)
    ) +
    geom_node_label(
      data = \(x) {get_nodes()(x) |> tidygraph::filter(type == "reaction")},
      aes(label = name, fill = value, color = ifelse((value - min(value)) / (max(value) - min(value)) > .5, "black", "white"))
    ) +
    labs(fill = "Depth") +
    scale_fill_viridis_c() +
    scale_color_identity() +
    NULL
}
