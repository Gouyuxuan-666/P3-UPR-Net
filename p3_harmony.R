suppressMessages({library(Seurat); library(harmony); library(ggplot2); library(patchwork)})

fig <- "/root/autodl-tmp/P3_supp_analysis"
s <- readRDS("/root/autodl-tmp/R_data/p3_seurat.rds")

# PCA（用现有 scale.data）
s <- RunPCA(s, npcs=30, verbose=FALSE)
cat("PCA 完成\n")

# 校正前 UMAP
s <- RunUMAP(s, dims=1:30, reduction="pca", reduction.name="umap_pre", verbose=FALSE)
cat("校正前 UMAP 完成\n")

# Harmony 批次校正
s <- RunHarmony(s, group.by.vars="sample", reduction.use="pca", verbose=FALSE)
cat("Harmony 完成\n")

# 校正后 UMAP
s <- RunUMAP(s, reduction="harmony", dims=1:30, reduction.name="umap_harmony", verbose=FALSE)
cat("校正后 UMAP 完成\n")

# ===== 出图 1: 校正前后 UMAP 对比 =====
p1 <- DimPlot(s, reduction="umap_pre", group.by="sample", pt.size=0.2) + ggtitle("Before Harmony (by sample)") + NoLegend()
p2 <- DimPlot(s, reduction="umap_harmony", group.by="sample", pt.size=0.2) + ggtitle("After Harmony (by sample)") + NoLegend()
p3 <- DimPlot(s, reduction="umap_pre", group.by="group", pt.size=0.2) + ggtitle("Before (by group)") + scale_color_manual(values=c(HC="#4575B4",DR="#B2182B",DS="#F18F01"))
p4 <- DimPlot(s, reduction="umap_harmony", group.by="group", pt.size=0.2) + ggtitle("After (by group)") + scale_color_manual(values=c(HC="#4575B4",DR="#B2182B",DS="#F18F01"))

ggsave(file.path(fig, "FigS5_harmony_umap.pdf"), (p1|p2)/(p3|p4), width=15, height=13)
cat("FigS5 完成\n")

# ===== 出图 2: SP110 表达在样本间分布（校正前后）=====
sp110 <- data.frame(SP110=s[["RNA"]]$data["SP110",], sample=s$sample, group=s$group)
p5 <- ggplot(sp110, aes(x=sample, y=SP110, fill=group)) + geom_violin(scale="width", alpha=0.7) +
  geom_boxplot(width=0.08, outlier.size=0.1) +
  scale_fill_manual(values=c(HC="#4575B4",DR="#B2182B",DS="#F18F01")) +
  labs(x="", y="SP110 expression", title="SP110 expression by sample (before batch correction)") +
  theme_bw(base_size=10) + theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave(file.path(fig, "FigS6_sp110_by_sample.pdf"), p5, width=10, height=5)
cat("FigS6 完成\n")

# ===== 保存 harmony 嵌入供因果发现 =====
harmony_emb <- Embeddings(s, "harmony")
saveRDS(harmony_emb, file.path(fig, "harmony_emb.rds"))
cat("harmony 嵌入已保存, 维度:", paste(dim(harmony_emb), collapse="x"), "\n")

cat("=== 批次校正完成 ===\n")
