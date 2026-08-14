# P3 Part 3: Network Rewiring + Pathway Enrichment
# Coulton 2024 Nat Commun fig4.R EXACT code patterns
library(Seurat)
library(dplyr)
library(ggplot2)
library(fgsea)
library(msigdbr)
library(org.Hs.eg.db)
library(ggrepel)
library(clusterProfiler)

OUT <- "/root/autodl-tmp/R_figures/Part3_Network_Rewiring"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- Load data ----
cat("Loading...\n")
cat("Loading...
")
library(Matrix)
cat("  Reading MTX...
")
m <- readMM("/root/autodl-tmp/R_data/counts_T.mtx")  # 24574 genes x 146345 cells
features <- readLines("/root/autodl-tmp/R_data/genes.csv")
barcodes <- readLines("/root/autodl-tmp/R_data/barcodes.csv")
rownames(m) <- features; colnames(m) <- barcodes
cat(sprintf("  Matrix: %d genes x %d cells
", nrow(m), ncol(m)))

meta <- read.csv("/root/autodl-tmp/R_data/obs_metadata.csv", row.names = 1)
cat(sprintf("  Metadata: %d cells, groups: %s
", nrow(meta), paste(unique(meta$group), collapse=", ")))

# Create Seurat with DR vs DS only (Coulton-style: compare two conditions)
cat("Creating Seurat (DR vs DS)...
")
dr_ds_cells <- intersect(rownames(meta)[meta$group %in% c("DR", "DS")], colnames(m))
cat(sprintf("  DR+DS cells: %d (DR=%d, DS=%d)
", length(dr_ds_cells),
           sum(meta[dr_ds_cells, "group"] == "DR"), sum(meta[dr_ds_cells, "group"] == "DS")))
so <- CreateSeuratObject(counts = m[, dr_ds_cells], meta.data = meta[dr_ds_cells, ])
so <- NormalizeData(so)
so <- FindVariableFeatures(so, nfeatures = 3000)
so <- ScaleData(so)
# ============================================================================
# Pseudobulk DEG: DR vs DS (Coulton fig4.R pattern)
# ============================================================================
cat('Computing pseudobulk DEG...\n')

sample_levels <- unique(so$sample)
pseudo_counts <- matrix(NA, nrow = nrow(so), ncol = length(sample_levels))
colnames(pseudo_counts) <- sample_levels
rownames(pseudo_counts) <- rownames(so)

sample_groups <- c()
for (i in seq_along(sample_levels)) {
    s <- sample_levels[i]
    cells <- colnames(so)[so$sample == s]
    pseudo_counts[, i] <- rowSums(GetAssayData(so, layer = 'counts')[, cells, drop = FALSE])
    sample_groups[i] <- unique(so$group[so$sample == s])
}
pseudo_counts <- pseudo_counts[rowSums(pseudo_counts) >= 10, ]

cat('  Testing DR vs DS...\n')
dr_samples <- colnames(pseudo_counts)[sample_groups == 'DR']
ds_samples <- colnames(pseudo_counts)[sample_groups == 'DS']

de_results <- data.frame(gene = rownames(pseudo_counts), logFC = NA, pval = NA)
dr_mat <- pseudo_counts[, dr_samples, drop = FALSE]
ds_mat <- pseudo_counts[, ds_samples, drop = FALSE]

for (i in seq_len(nrow(pseudo_counts))) {
    if (i %% 1000 == 0) cat(sprintf('  gene %d/%d\n', i, nrow(pseudo_counts)))
    dr_vals <- dr_mat[i, ]; ds_vals <- ds_mat[i, ]
    dr_mean <- mean(dr_vals); ds_mean <- mean(ds_vals)
    de_results$logFC[i] <- log2((ds_mean + 1) / (dr_mean + 1))
    tryCatch({
        w <- wilcox.test(dr_vals, ds_vals)
        de_results$pval[i] <- w$p.value
    }, error = function(e) { de_results$pval[i] <- NA })
}
de_results$padj <- p.adjust(de_results$pval, method = 'fdr')
de_results <- de_results[!is.na(de_results$padj), ]
de_results <- de_results[order(de_results$pval), ]
cat(sprintf('  %d DEGs (padj < 0.05)\n', sum(de_results$padj < 0.05, na.rm = TRUE)))

# ============================================================================
# Fig3A: Enhanced Volcano Plot (Burger 2024 Nat Commun standard)
# ============================================================================
cat('Fig3A: Enhanced volcano...\n')

de_results$sig <- 'NS'
de_results$sig[de_results$padj < 0.05 & abs(de_results$logFC) > 0.5] <- 'Significant'
de_results$sig[de_results$padj < 0.01 & abs(de_results$logFC) > 1.0] <- 'Highly sig'

sp110_targets <- c('SP110', 'HSPA5', 'ATF4', 'ATF6', 'DDIT3', 'XBP1', 'ERN1',
                   'EIF2AK3', 'BECN1', 'MAP1LC3B', 'SQSTM1', 'BAX', 'BCL2',
                   'CASP3', 'MAPK8', 'NFKB1', 'GBP1', 'GBP5', 'TRIB3', 'PPP1R15A')
de_results$label <- ifelse(de_results$gene %in% sp110_targets, de_results$gene, '')

n_up <- sum(de_results$sig != 'NS' & de_results$logFC > 0, na.rm = TRUE)
n_dn <- sum(de_results$sig != 'NS' & de_results$logFC < 0, na.rm = TRUE)
n_high <- sum(de_results$sig == 'Highly sig', na.rm = TRUE)

p_volcano <- ggplot(de_results, aes(x = logFC, y = -log10(pval))) +
  geom_point(aes(color = sig), size = 0.6, alpha = 0.7) +
  scale_color_manual(
    values = c('NS' = 'grey85', 'Significant' = '#377EB8', 'Highly sig' = '#E41A1C'),
    labels = c('NS' = 'NS',
               'Significant' = sprintf('p<0.05 & |FC|>0.5 (n=%d)', n_up + n_dn),
               'Highly sig' = sprintf('p<0.01 & |FC|>1 (n=%d)', n_high))) +
  geom_vline(xintercept = c(-1, -0.5, 0.5, 1),
             linetype = c('dotted', 'dashed', 'dashed', 'dotted'),
             linewidth = 0.25, color = 'grey40') +
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01)),
             linetype = c('dashed', 'dotted'), linewidth = 0.25, color = 'grey40') +
  geom_text_repel(aes(label = label), size = 3.2, max.overlaps = 25,
                  box.padding = 0.4, min.segment.length = 0.1,
                  segment.size = 0.2, fontface = 'italic') +
  annotate('text', x = max(de_results$logFC, na.rm = TRUE)*0.8,
           y = max(-log10(de_results$pval), na.rm = TRUE)*0.95,
           label = sprintf('DS up: %d', n_up), hjust = 1, size = 4.5,
           color = '#E41A1C', fontface = 'bold') +
  annotate('text', x = min(de_results$logFC, na.rm = TRUE)*0.8,
           y = max(-log10(de_results$pval), na.rm = TRUE)*0.95,
           label = sprintf('DR up: %d', n_dn), hjust = 0, size = 4.5,
           color = '#377EB8', fontface = 'bold') +
  theme_classic(base_size = 14) +
  labs(x = expression(log[2]*'(DS/DR)'), y = expression(-log[10]*'(p-value)'),
       color = '') +
  theme(
          legend.position = 'bottom',
        legend.key.size = unit(3, 'mm'))
ggsave(file.path(OUT, 'Fig3A_volcano.png'), p_volcano, width = 9, height = 8, dpi = 300)

# ============================================================================
# Fig3B: fGSEA Enhanced Bubble + Bar
# ============================================================================
cat('Fig3B: fGSEA enhanced bubble...\n')

reactome_paths <- msigdbr(species = 'Homo sapiens', collection = 'C2', subcollection = 'CP:REACTOME')
reactome_list <- split(x = reactome_paths$gene_symbol, f = reactome_paths$gs_name)
reactome_sizes <- sapply(reactome_list, length)
reactome_list <- reactome_list[reactome_sizes >= 15 & reactome_sizes <= 500]

ranks <- de_results$logFC
names(ranks) <- de_results$gene
ranks <- ranks[!is.na(ranks)]
ranks <- ranks[order(ranks, decreasing = TRUE)]
ranks <- ranks[!duplicated(names(ranks))]

cat('  Running fgsea...\n')
fgseaRes <- fgsea(pathways = reactome_list, stats = ranks, minSize = 15, maxSize = 500)
fgseaRes <- fgseaRes[order(fgseaRes$pval), ]
sig_paths <- fgseaRes[padj < 0.1]
cat(sprintf('  %d significant pathways (padj < 0.1)\n', nrow(sig_paths)))

if (nrow(sig_paths) > 0) {
    bubble_df <- as.data.frame(sig_paths)[, c('pathway', 'padj', 'NES', 'size')]
    colnames(bubble_df)[4] <- 'gene_set_size'
    bubble_df <- bubble_df[order(bubble_df$NES), ]
    bubble_df <- head(bubble_df[order(-abs(bubble_df$NES)), ], 30)
    bubble_df <- bubble_df[order(bubble_df$NES), ]

    bubble_df$pathway_label <- gsub('REACTOME_', '', bubble_df$pathway)
    bubble_df$pathway_label <- gsub('_', ' ', bubble_df$pathway_label)
    bubble_df$pathway_label <- sapply(bubble_df$pathway_label, function(x) {
        if (nchar(x) > 65) paste0(substr(x, 1, 62), '...') else x
    })
    bubble_df$pathway_label <- factor(bubble_df$pathway_label,
                                       levels = bubble_df$pathway_label)

    p_bubble <- ggplot(bubble_df, aes(x = NES, y = pathway_label)) +
      geom_point(aes(size = gene_set_size, fill = -log10(padj)),
                 shape = 21, alpha = 0.85) +
      scale_fill_gradientn(
        colors = c('#2166AC', '#92C5DE', '#F7F7F7', '#F4A582', '#B2182B'),
        values = scales::rescale(c(0, 0.5, 1, 1.5, max(-log10(bubble_df$padj), na.rm = TRUE))),
        name = expression(-log[10]*'(FDR)')) +
      scale_size_continuous(name = 'Gene set size', range = c(2, 8)) +
      geom_vline(xintercept = 0, linetype = 'dashed', linewidth = 0.3, color = 'grey50') +
      theme_classic(base_size = 11) +
      labs(x = 'Normalized Enrichment Score', y = '',
      ) +
      theme(
                  legend.position = 'right',
            axis.text.y = element_text(size = 8))
    ggsave(file.path(OUT, 'Fig3B_fgsea_bubble.png'), p_bubble, width = 13, height = 9, dpi = 300)

    top20 <- head(bubble_df[order(-abs(bubble_df$NES)), ], 20)
    top20 <- top20[order(top20$NES), ]
    top20$pathway_label <- factor(top20$pathway_label, levels = top20$pathway_label)
    top20$direction <- ifelse(top20$NES > 0, 'DS', 'DR')

    p_bar <- ggplot(top20, aes(x = NES, y = pathway_label, fill = direction)) +
      geom_col(width = 0.7, alpha = 0.9) +
      scale_fill_manual(values = c('DR' = '#2166AC', 'DS' = '#B2182B'),
                        labels = c('DR' = 'DR up', 'DS' = 'DS up')) +
      geom_vline(xintercept = 0, linewidth = 0.4) +
      theme_classic(base_size = 11) +
      labs(x = 'Normalized Enrichment Score', y = '', fill = '',
      ) +
      theme(
            legend.position = 'bottom',
            axis.text.y = element_text(size = 8.5))
    ggsave(file.path(OUT, 'Fig3B_fgsea_barplot.png'), p_bar, width = 12, height = 7, dpi = 300)
}

# ============================================================================
# Fig3C: Per-Sample Pseudobulk Heatmap (pheatmap)
# ============================================================================
cat('Fig3C: pheatmap...\n')

upr_genes <- intersect(c('SP110', 'HSPA5', 'ERN1', 'EIF2AK3', 'ATF6', 'XBP1', 'ATF4',
                          'DDIT3', 'BECN1', 'MAP1LC3B', 'SQSTM1', 'CASP3', 'BAX', 'BCL2',
                          'MAPK8', 'NFKB1', 'GBP1', 'GBP5', 'TNF', 'IL6', 'CXCL9',
                          'PPP1R15A', 'TRIB3', 'EIF2S1', 'ATF3'), rownames(so))

if (length(upr_genes) > 5) {
    library(pheatmap)

    sample_levels <- unique(so$sample)
    heat_data <- matrix(NA, nrow = length(upr_genes), ncol = length(sample_levels))
    rownames(heat_data) <- upr_genes
    colnames(heat_data) <- sample_levels

    for (i in seq_along(sample_levels)) {
        s <- sample_levels[i]
        cells <- colnames(so)[so$sample == s]
        heat_data[, i] <- rowMeans(GetAssayData(so, layer = 'data')[upr_genes, cells, drop = FALSE])
    }

    heat_data_z <- t(scale(t(heat_data)))

    ann_df <- data.frame(
      Group = ifelse(grepl('^DR', sample_levels), 'Drug Resistant', 'Drug Sensitive'),
      row.names = sample_levels
    )
    ann_colors <- list(Group = c('Drug Resistant' = '#A23B72', 'Drug Sensitive' = '#F18F01'))

    png(file.path(OUT, 'Fig3C_UPR_rewiring_heatmap.png'), width = 2800, height = 2400, res = 300)
    pheatmap(heat_data_z,
             annotation_col = ann_df,
             annotation_colors = ann_colors,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             scale = 'none',
             color = colorRampPalette(c('#2166AC', 'white', '#B2182B'))(100),
             border_color = NA,
             fontsize = 10,
             fontsize_row = 9,
             fontsize_col = 8,
             angle_col = 45,
             legend_title = 'Z-score',
             display_numbers = FALSE)
    dev.off()
}

# ============================================================================
# Fig3D: Top DEG Horizontal Barplot (with significance stars)
# ============================================================================
cat('Fig3D: Top DEG barplot...\n')

de_sig <- de_results[de_results$padj < 0.05 & abs(de_results$logFC) > 0.5, ]
top_dr <- head(de_sig[order(de_sig$logFC), ], 15)
top_ds <- head(de_sig[order(-de_sig$logFC), ], 15)
top_deg <- rbind(
    data.frame(gene = top_dr$gene, logFC = top_dr$logFC, group = 'DR',
               padj = top_dr$padj, stringsAsFactors = FALSE),
    data.frame(gene = top_ds$gene, logFC = top_ds$logFC, group = 'DS',
               padj = top_ds$padj, stringsAsFactors = FALSE)
)
top_deg <- top_deg[order(top_deg$logFC), ]
top_deg$gene <- factor(top_deg$gene, levels = top_deg$gene)
top_deg$sig_stars <- ifelse(top_deg$padj < 0.001, '***',
                     ifelse(top_deg$padj < 0.01, '**',
                     ifelse(top_deg$padj < 0.05, '*', '')))

p_topdeg <- ggplot(top_deg, aes(x = logFC, y = gene)) +
  geom_col(aes(fill = group), width = 0.7, alpha = 0.9) +
  geom_text(aes(label = sig_stars,
                x = ifelse(logFC > 0, logFC + 0.15, logFC - 0.15)),
            size = 5, vjust = 0.75, color = 'grey20') +
  geom_text(aes(label = gene,
                x = ifelse(logFC > 0, -0.05, 0.05)),
            size = 3.2, hjust = ifelse(top_deg$logFC > 0, 1, 0),
            fontface = 'italic') +
  scale_fill_manual(values = c('DR' = '#2166AC', 'DS' = '#B2182B'),
                    labels = c('DR' = 'Drug Resistant', 'DS' = 'Drug Sensitive')) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  theme_classic(base_size = 13) +
  labs(x = expression(log[2]*'(DS/DR)'), y = '',
       fill = '') +
  theme(
          legend.position = 'bottom',
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave(file.path(OUT, 'Fig3D_top_DEG_barplot.png'), p_topdeg, width = 8, height = 9, dpi = 300)

cat('\n=== Part 3 Complete (NC-standard) ===\n')
for (f in sort(list.files(OUT, pattern = '*.png'))) {
    info <- file.info(file.path(OUT, f))
    cat(sprintf('  %s (%.0f KB)\n', f, info$size / 1024))
}
