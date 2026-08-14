suppressMessages({library(pcalg); library(Seurat); library(reshape2)})

fig <- "/root/autodl-tmp/P3_supp_analysis"
genes <- c("SP110","EIF2AK3","ERN1","ATF6","HSPA5","ATF4","DDIT3","XBP1","TRAF2",
           "MAPK8","NFKB1","ATF6B","PPP1R15A","TRIB3","ATF3","EIF2S1","MBTPS1",
           "MBTPS2","BAX","BCL2","BCL2L1","CASP3","CASP9","BECN1","ATG5","SQSTM1",
           "MAP1LC3B","NFE2L2","HMOX1","MTOR")

s <- readRDS("/root/autodl-tmp/R_data/p3_seurat.rds")
set.seed(42)
sub <- sample(colnames(s)[s$group == "DR"], 8000)
expr <- as.matrix(GetAssayData(s, layer="data")[genes, sub])
expr <- expr[complete.cases(expr), ]
rm(s); gc()

n_genes <- nrow(expr); n_cells <- ncol(expr)
B <- 100
edge_count <- matrix(0, n_genes, n_genes, dimnames=list(rownames(expr), rownames(expr)))
t0 <- Sys.time()

for (b in 1:B) {
  set.seed(1000 + b)
  idx <- sample(n_cells, replace=TRUE)
  suffStat <- list(C = cor(t(expr[, idx])), n = length(idx))
  pc_b <- pc(suffStat=suffStat, labels=rownames(expr), indepTest=gaussCItest, alpha=0.01, skel.method="stable")
  amat_b <- as(pc_b, "amat")
  edge_count <- edge_count + (amat_b != 0)
  if (b %% 10 == 0) {
    el <- as.numeric(difftime(Sys.time(), t0, units="mins"))
    cat(sprintf("bootstrap %d/%d 完成, 累计 %.1f 分钟\n", b, B, el))
  }
}

edge_prob <- edge_count / B
saveRDS(edge_prob, file.path(fig, "bootstrap_edge_prob.rds"))

# 出图：边稳定概率热图
pdf(file.path(fig, "FigS4_bootstrap_edge_prob.pdf"), width=10, height=9)
pheatmap::pheatmap(edge_prob, cluster_rows=TRUE, cluster_cols=TRUE,
  color=colorRampPalette(c("white","#FDDBC7","#B2182B"))(100),
  main="Edge stability probability (100 bootstrap, DR n=8000)", fontsize=8)
dev.off()

# SP110 边稳定概率
sp110_out <- sort(edge_prob["SP110", ], decreasing=TRUE)
cat("\n=== SP110 出边稳定概率 ===\n")
print(round(sp110_out[sp110_out > 0.3], 3))
sp110_in <- sort(edge_prob[, "SP110"], decreasing=TRUE)
cat("\n=== SP110 入边稳定概率 ===\n")
print(round(sp110_in[sp110_in > 0.3], 3))

# 稳定边
stable <- which(edge_prob > 0.5, arr.ind=TRUE)
cat("\n=== 稳定边(概率>0.5) ===\n")
if (nrow(stable) > 0) {
  for (i in 1:nrow(stable)) {
    cat(rownames(edge_prob)[stable[i,1]], "->", colnames(edge_prob)[stable[i,2]], 
        "=", round(edge_prob[stable[i,1], stable[i,2]], 2), "\n")
  }
} else cat("无\n")

# NOTEARS 3 边的稳定概率
cat("\n=== NOTEARS 3 边稳定概率 ===\n")
for (e in list(c("CASP3","HSPA5"), c("HSPA5","XBP1"), c("EIF2S1","HSPA5"))) {
  p <- edge_prob[e[1], e[2]] + edge_prob[e[2], e[1]]
  cat(e[1], "<->", e[2], ":", round(p, 2), "\n")
}

cat("\n=== Bootstrap 100 次完成, 总耗时 %.1f 分钟 ===\n", as.numeric(difftime(Sys.time(), t0, units="mins")))
