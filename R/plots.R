# ps to data table:
# d <- merge(
#   melt(
#     merge(
#       as.data.table(otu_table(ps))[, ID := rownames(otu_table(ps))],
#       as.data.table(tax_table(ps)@.Data)[, ID := rownames(tax_table(ps))],
#       by = "ID"
#     ),
#     id.vars = c("ID", colnames(tax_table(ps))),
#     variable.name = "Sample",
#     value.name = "Depth"
#   ),
#   setattr(as.data.table(sample_data(ps)@.Data), "names", colnames(sample_data(ps)))[, Sample := rownames(sample_data(ps))],
#   by = "Sample"
# )

depth_by_sample <- function(out) {
  out <- merge(out, setattr(as.data.table(sample_data(ps)@.Data), "names", colnames(sample_data(ps)))[, Sample := rownames(sample_data(ps))], by = "Sample")
  ggplot(out[, .(RelativeDepth = sum(RelativeDepth)), by = c("Sample", "Taxon")], aes(x = Sample, y = RelativeDepth, fill = factor(Sample))) +
    geom_col()
}

depth_by_sample_group <- function(out, ps, group = NA) {
  out <- merge(out, setattr(as.data.table(sample_data(ps)@.Data), "names", colnames(sample_data(ps)))[, Sample := rownames(sample_data(ps))], by = "Sample")
  ggplot(out2[, .(RelativeDepth = sum(RelativeDepth)), by = .(Group = group, Taxon = Taxon)], aes(x = factor(Group), y = RelativeDepth, fill = Taxon)) +
    geom_boxplot()
}

depth_by_gene <- function(out_list) {

}

percent_depth_by_gene <- function() {

}
