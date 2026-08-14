suppressMessages({library(pheatmap)})

fig <- "/root/autodl-tmp/R_figures/Fig3_single"
dir.create(fig, recursive=TRUE, showWarnings=FALSE)

Nred <- "#B2182B"; Nblue <- "#2166AC"

# Fig3A: PC+GES 边一致性热图（精细化）
r <- readRDS("/root/autodl-tmp/P3_supp_analysis/causal_comparison.rds")
mat_cons <- (r$pc != 0) + (r$ges != 0)  # 0/1/2

png(file.path(fig,"Fig3A_method_consensus.png"), width=9, height=8, units="in", res=300)
pheatmap(mat_cons,
  color=c("#FFFFFF","#FDDBC7","#B2182B"),
  breaks=c(-0.5,0.5,1.5,2.5),
  legend_breaks=c(0,1,2),
  legend_labels=c("None","One method","Both (PC+GES)"),
  cluster_rows=TRUE, cluster_cols=TRUE,
  clustering_method="complete",
  fontsize=9, fontsize_row=8, fontsize_col=8,
  border_color=NA,
  legend=TRUE)
dev.off()
cat("Fig3A 完成\n")

# Fig3B: Bootstrap 边稳定概率热图（精细化）
b <- readRDS("/root/autodl-tmp/P3_supp_analysis/bootstrap_edge_prob.rds")

png(file.path(fig,"Fig3B_bootstrap_edge_prob.png"), width=9, height=8, units="in", res=300)
pheatmap(b,
  color=colorRampPalette(c("#FFFFFF","#FDDBC7","#B2182B"))(100),
  cluster_rows=TRUE, cluster_cols=TRUE,
  clustering_method="complete",
  fontsize=9, fontsize_row=8, fontsize_col=8,
  border_color=NA,
  legend=TRUE)
dev.off()
cat("Fig3B 完成\n")
