suppressMessages({library(ggplot2); library(pheatmap); library(patchwork)})

fig <- "/root/autodl-tmp/radiomics_R_figs"
dir.create(fig, recursive=TRUE, showWarnings=FALSE)

img <- read.csv("/root/autodl-tmp/radiomics_features.csv", check.names=FALSE)
mol <- read.csv("/root/autodl-tmp/molecular_scores_7.csv", check.names=FALSE)

exclude <- c("original_firstorder_Energy","original_firstorder_TotalEnergy")
img_feats <- setdiff(colnames(img)[4:ncol(img)], exclude)
mol_genes <- colnames(mol)[4:ncol(mol)]

dat <- merge(img, mol, by="sample")
grp <- dat$group.x

# 排除 constant 特征（方差=0）
sds <- sapply(img_feats, function(f) sd(dat[[f]], na.rm=TRUE))
img_feats <- img_feats[sds > 1e-6]
cat("保留特征:", length(img_feats), "个\n")

# ===== 1. 关联热图 =====
corr <- matrix(0, length(img_feats), length(mol_genes),
               dimnames=list(gsub("original_firstorder_","",img_feats), mol_genes))
for (i in seq_along(img_feats))
  for (j in seq_along(mol_genes))
    corr[i,j] <- cor(dat[[img_feats[i]]], dat[[mol_genes[j]]], method="spearman")

pdf(file.path(fig,"radiomics_molecular_corr.pdf"), width=10, height=8)
pheatmap(corr, color=colorRampPalette(c("#2166AC","white","#B2182B"))(100),
         breaks=seq(-1,1,length.out=101), fontsize=9, border_color=NA,
         cluster_rows=FALSE, cluster_cols=FALSE,
         main="Imaging-molecular correlation (Spearman, n=7)")
dev.off()
cat("热图完成\n")

# ===== 2. 散点图 =====
mk_scatter <- function(feat, label) {
  d <- data.frame(x=dat[[feat]], y=dat[["ERN1"]], g=grp)
  rho <- cor(d$x, d$y, method="spearman")
  p <- cor.test(d$x, d$y, method="spearman")$p.value
  ggplot(d, aes(x=x, y=y, color=g)) + geom_point(size=3.5) +
    geom_smooth(method="lm", se=FALSE, color="grey30", linetype=2, linewidth=0.7) +
    scale_color_manual(values=c(NY="#B2182B",YM="#2166AC"), name="") +
    labs(x=label, y="ERN1 expression", title=sprintf("rho=%.2f, P=%.3f", rho, p)) +
    theme_bw(base_size=11)
}
p1 <- mk_scatter("original_firstorder_InterquartileRange","HU interquartile range")
p2 <- mk_scatter("original_firstorder_Uniformity","HU uniformity")
ggsave(file.path(fig,"radiomics_molecular_scatter.pdf"), p1+p2, width=10, height=4.8)
cat("散点图完成\n")

# ===== 3. 箱线图 =====
key <- c("original_firstorder_Mean","original_firstorder_Median","original_firstorder_Skewness",
         "original_firstorder_Kurtosis","original_firstorder_Entropy","original_firstorder_Uniformity")
ld <- do.call(rbind, lapply(key, function(f) {
  data.frame(value=dat[[f]], feature=gsub("original_firstorder_","",f), group=grp)
}))
ld$feature <- factor(ld$feature, levels=gsub("original_firstorder_","",key))
p3 <- ggplot(ld, aes(x=group, y=value, fill=group)) +
  geom_boxplot(alpha=0.6, outlier.shape=NA) +
  geom_jitter(width=0.12, size=2, alpha=0.8) +
  facet_wrap(~feature, scales="free_y", ncol=3) +
  scale_fill_manual(values=c(NY="#B2182B",YM="#2166AC"), name="") +
  labs(x="", y="", title="Whole-lung radiomics features: DR vs DS") +
  theme_bw(base_size=10) + theme(legend.position="top")
ggsave(file.path(fig,"radiomics_feature_boxplot.pdf"), p3, width=10, height=6)
cat("箱线图完成\n")

# ===== 4. PCA =====
X <- scale(as.matrix(dat[, img_feats]))
pca <- prcomp(X, center=FALSE)
pcadf <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=grp, sample=dat$sample)
var1 <- summary(pca)$importance[2,1]*100; var2 <- summary(pca)$importance[2,2]*100
p4 <- ggplot(pcadf, aes(x=PC1, y=PC2, color=group)) +
  geom_point(size=4) + geom_text(aes(label=sample), vjust=-1.2, size=3) +
  scale_color_manual(values=c(NY="#B2182B",YM="#2166AC"), name="") +
  labs(x=sprintf("PC1 (%.1f%%)",var1), y=sprintf("PC2 (%.1f%%)",var2),
       title="Whole-lung radiomics PCA (n=7)") +
  theme_bw(base_size=12) + theme(legend.position="top")
ggsave(file.path(fig,"radiomics_pca.pdf"), p4, width=6.5, height=5.5)
cat("PCA 完成\n")

cat("=== 全部完成 ===\n"); print(list.files(fig))
