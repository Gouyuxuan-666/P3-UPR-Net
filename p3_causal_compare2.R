suppressMessages({library(pcalg); library(igraph); library(ggplot2); library(Seurat); library(reshape2)})

fig <- "/root/autodl-tmp/P3_supp_analysis"
dir.create(fig, recursive=TRUE, showWarnings=FALSE)

genes <- c("SP110","EIF2AK3","ERN1","ATF6","HSPA5","ATF4","DDIT3","XBP1","TRAF2",
           "MAPK8","NFKB1","ATF6B","PPP1R15A","TRIB3","ATF3","EIF2S1","MBTPS1",
           "MBTPS2","BAX","BCL2","BCL2L1","CASP3","CASP9","BECN1","ATG5","SQSTM1",
           "MAP1LC3B","NFE2L2","HMOX1","MTOR")

s <- readRDS("/root/autodl-tmp/R_data/p3_seurat.rds")
set.seed(42)
sub <- sample(colnames(s)[s$group == "DR"], 8000)
expr <- as.matrix(GetAssayData(s, layer="data")[genes, sub])
expr <- expr[complete.cases(expr), ]
cat("表达矩阵:", nrow(expr), "x", ncol(expr), "\n")

# ===== PC =====
suffStat <- list(C = cor(t(expr)), n = ncol(expr))
pc.fit <- pc(suffStat=suffStat, labels=rownames(expr), indepTest=gaussCItest, alpha=0.01, skel.method="stable")
amat_pc <- as(pc.fit, "amat")
cat("PC 边数:", sum(amat_pc != 0), "\n")

# ===== GES =====
score <- new("GaussL0penObsScore", data = t(expr))
fit_ges <- ges(score)
nodes_ges <- fit_ges$repr$.nodes
amat_ges <- matrix(0, length(nodes_ges), length(nodes_ges), dimnames=list(nodes_ges, nodes_ges))
for (j in seq_along(nodes_ges)) {
  p <- fit_ges$repr$.in.edges[[j]]
  if (length(p) > 0) amat_ges[p, nodes_ges[j]] <- 1
}
cat("GES 边数:", sum(amat_ges != 0), "\n")

extract_edges <- function(amat, labels) {
  e <- which(amat != 0, arr.ind=TRUE)
  if (nrow(e)==0) return(data.frame(from=character(0), to=character(0)))
  data.frame(from=labels[e[,1]], to=labels[e[,2]])
}
pc_edges <- extract_edges(amat_pc, rownames(expr))
ges_edges <- extract_edges(amat_ges, rownames(expr))

# ===== NOTEARS 边恢复验证 =====
notears <- data.frame(from=c("CASP3","HSPA5","EIF2S1"), to=c("HSPA5","XBP1","HSPA5"))
cat("\n=== NOTEARS 3 边恢复验证 ===\n")
for (i in 1:nrow(notears)) {
  f <- notears$from[i]; t <- notears$to[i]
  in_pc <- any((pc_edges$from==f & pc_edges$to==t) | (pc_edges$from==t & pc_edges$to==f))
  in_ges <- any((ges_edges$from==f & ges_edges$to==t) | (ges_edges$from==t & ges_edges$to==f))
  cat(sprintf("  %s->%s : PC=%s, GES=%s\n", f, t, in_pc, in_ges))
}

cat("\n=== SP110 因果边 ===\n")
cat("PC 出边:", paste(pc_edges$to[pc_edges$from=="SP110"], collapse=","), "\n")
cat("GES 出边:", paste(ges_edges$to[ges_edges$from=="SP110"], collapse=","), "\n")

# ===== 出图 =====
g_pc <- graph_from_adjacency_matrix(amat_pc != 0, mode="directed")
V(g_pc)$color <- ifelse(names(V(g_pc))=="SP110", "#B2182B", "#4575B4")
pdf(file.path(fig, "FigS1_PC_network.pdf"), width=9.5, height=7.5)
plot(g_pc, layout=layout_with_fr(g_pc), vertex.size=9, edge.arrow.size=0.35,
     vertex.label.cex=0.8, vertex.label.color="black", vertex.frame.color="grey30",
     main="PC algorithm causal network (DR, n=8000)")
dev.off()

g_ges <- graph_from_adjacency_matrix(amat_ges != 0, mode="directed")
V(g_ges)$color <- ifelse(names(V(g_ges))=="SP110", "#B2182B", "#4575B4")
pdf(file.path(fig, "FigS2_GES_network.pdf"), width=9.5, height=7.5)
plot(g_ges, layout=layout_with_fr(g_ges), vertex.size=9, edge.arrow.size=0.35,
     vertex.label.cex=0.8, vertex.label.color="black", vertex.frame.color="grey30",
     main="GES causal network (DR, n=8000)")
dev.off()

mat_cons <- matrix(0, length(genes), length(genes), dimnames=list(genes, genes))
if (nrow(pc_edges)>0) for (i in 1:nrow(pc_edges)) mat_cons[pc_edges$from[i], pc_edges$to[i]] <- mat_cons[pc_edges$from[i], pc_edges$to[i]] + 1
if (nrow(ges_edges)>0) for (i in 1:nrow(ges_edges)) mat_cons[ges_edges$from[i], ges_edges$to[i]] <- mat_cons[ges_edges$from[i], ges_edges$to[i]] + 1

pdf(file.path(fig, "FigS3_method_consensus.pdf"), width=9, height=8)
pheatmap::pheatmap(mat_cons, cluster_rows=FALSE, cluster_cols=FALSE,
  color=c("white","#FDDBC7","#B2182B"), breaks=c(-0.5,0.5,1.5,2.5),
  legend_breaks=c(0,1,2), legend_labels=c("无","单方法","双方法一致"),
  main="Causal edge consensus (PC + GES)", fontsize=9)
dev.off()

saveRDS(list(pc=amat_pc, ges=amat_ges, pc_edges=pc_edges, ges_edges=ges_edges),
        file.path(fig, "causal_comparison.rds"))
cat("\n=== 全部完成 ===\n")
