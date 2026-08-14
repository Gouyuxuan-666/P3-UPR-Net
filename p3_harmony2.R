suppressMessages({library(Seurat); library(harmony); library(ggplot2); library(patchwork)})

fig <- "/root/autodl-tmp/P3_supp_analysis"
s <- readRDS("/root/autodl-tmp/R_data/p3_seurat.rds")

# 降采样：每 group 8000 细胞
set.seed(42)
cells <- unlist(lapply(c("HC","DR","DS"), function(g) sample(colnames(s)[s$group==g], 8000)))
s_sub <- subset(s, cells = cells)
cat("降采样:", ncol(s_sub), "细胞\n")

s_sub <- RunPCA(s_sub, npcs=30, verbose=FALSE)
cat("PCA 完成\n")
s_sub <- RunUMAP(s_sub, dims=1:30, reduction="pca", reduction.name="umap_pre", verbose=FALSE)
cat("校正前 UMAP 完成\n")
s_sub <- RunHarmony(s_sub, group.by.vars="sample", reduction.use="pca", verbose=FALSE)
cat("Harmony 完成\n")
s_sub <- RunUMAP(s_sub, reduction="harmony", dims=1:30, reduction.name="umap_harmony", verbose=FALSE)
cat("校正后 UMAP 完成\n")

p1 <- DimPlot(s_sub, reduction="umap_pre", group.by="sample", pt.size=0.3) + ggtitle("Before Harmony (by sample)") + NoLegend()
p2 <- DimPlot(s_sub, reduction="umap_harmony", group.by="sample", pt.size=0.3) + ggtitle("After Harmony (by sample)") + NoLegend()
p3 <- DimPlot(s_sub, reduction="umap_pre", group.by="group", pt.size=0.3) + ggtitle("Before (by group)") + scale_color_manual(values=c(HC="#4575B4",DR="#B2182B",DS="#F18F01"))
p4 <- DimPlot(s_sub, reduction="umap_harmony", group.by="group", pt.size=0.3) + ggtitle("After (by group)") + scale_color_manual(values=c(HC="#4575B4",DR="#B2182B",DS="#F18F01"))

ggsave(file.path(fig, "FigS5_harmony_umap.pdf"), (p1|p2)/(p3|p4), width=15, height=13)
cat("FigS5 完成\n")

sp110 <- data.frame(SP110=s_sub[["RNA"]]$data["SP110",], sample=s_sub$sample, group=s_sub$group)
p5 <- ggplot(sp110, aes(x=sample, y=SP110, fill=group)) + geom_violin(scale="width", alpha=0.7) +
  geom_boxplot(width=0.08, outlier.size=0.1) +
  scale_fill_manual(values=c(HC="#4575B4",DR="#B2182B",DS="#F18F01")) +
  labs(x="", y="SP110 expression", title="SP110 expression by sample") +
  theme_bw(base_size=10) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(fig, "FigS6_sp110_by_sample.pdf"), p5, width=10, height=5)
cat("FigS6 完成\n")

cat("=== 批次校正完成 ===\n")
