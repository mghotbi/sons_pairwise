# install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap")
# install.packages("circlize")

library(ComplexHeatmap)
library(circlize)

hc <- hclust(as.dist(D), method = "average")  # UPGMA-style

col_fun <- colorRamp2(
  c(0, median(D), max(D)),
  c("#08306B", "#FDBE85", "#67000D")
)

ht <- Heatmap(
  D,
  name = "SONS distance\n(1 - θYC)",
  col = col_fun,
  cluster_rows = hc,
  cluster_columns = hc,
  rect_gp = grid::gpar(col = "white", lwd = 0.6),
  row_names_gp = grid::gpar(fontsize = 10),
  column_names_gp = grid::gpar(fontsize = 10),
  heatmap_legend_param = list(
    title_gp = grid::gpar(fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 9)
  )
)

draw(ht)


pdf("sons_thetaYC_heatmap.pdf", width = 7.5, height = 6.5, useDingbats = FALSE)
draw(ht)
dev.off()

png("sons_thetaYC_heatmap.png", width = 2400, height = 2000, res = 300)
draw(ht)
dev.off()


sample_groups <- data.frame(
  Group = rep(c("Zoo", "Wild"), each = 4),
  row.names = colnames(otu_mat)
)

ha <- HeatmapAnnotation(
  df = sample_groups,
  col = list(Group = c(Zoo = "#1b9e77", Wild = "#d95f02"))
)

ht2 <- Heatmap(
  D,
  name = "SONS distance\n(1 - θYC)",
  col = col_fun,
  cluster_rows = hc,
  cluster_columns = hc,
  top_annotation = ha
)

draw(ht2)
