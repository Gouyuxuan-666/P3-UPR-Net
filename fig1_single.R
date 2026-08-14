suppressMessages({library(ggplot2); library(dplyr); library(viridis); library(ggsci)})

fig <- "/root/autodl-tmp/R_figures/Fig1_single"
dir.create(fig, recursive=TRUE, showWarnings=FALSE)

d <- read.csv("/root/autodl-tmp/R_data/p3_mac_clus_full.csv", row.names=1, check.names=FALSE)
d$group <- factor(d$group, levels=c("HC","DR","DS"))
cols <- c(HC="#4575B4", DR="#B2182B", DS="#F18F01")

# 精细图例主题
thm <- theme_classic(base_size=13) +
  theme(axis.title=element_text(size=12),
        axis.text=element_text(size=10),
        legend.title=element_text(size=10, face="bold"),
        legend.text=element_text(size=9),
        legend.background=element_rect(fill=alpha("white",0.7), color=NA),
        legend.key=element_blank(),
        panel.border=element_rect(fill=NA, color="grey30", linewidth=0.5))

set.seed(42); dd <- d[sample(nrow(d)), ]

# Fig1A: UMAP by group, 图例图内右上角
p1 <- ggplot(dd, aes(x=UMAP_score, y=UMAP_2, color=group)) +
  geom_point(size=0.15, alpha=0.55) +
  scale_color_manual(values=cols, name=NULL) +
  labs(x="UMAP1", y="UMAP2") +
  guides(color=guide_legend(override.aes=list(size=4, alpha=1), ncol=1)) +
  thm + theme(legend.position=c(0.98,0.98), legend.justification=c(1,1))
ggsave(file.path(fig,"Fig1A_UMAP_by_group.png"), p1, width=7.5, height=6, dpi=300)

# Fig1B: SP110 UMAP, colorbar 精细
p2 <- ggplot(dd, aes(x=UMAP_score, y=UMAP_2, color=SP110_expr)) +
  geom_point(size=0.15, alpha=0.6) +
  scale_color_viridis_c(option="magma", name="SP110", limits=c(0, quantile(d$SP110_expr,0.98)),
                        breaks=pretty(c(0, quantile(d$SP110_expr,0.98)), 4)) +
  labs(x="UMAP1", y="UMAP2") +
  thm + theme(legend.position="right", legend.key.height=unit(1.1,"cm"))
ggsave(file.path(fig,"Fig1B_SP110_UMAP.png"), p2, width=7.5, height=6, dpi=300)

# Fig1C: dotplot, size+color 双图例精细
dot <- read.csv("/root/autodl-tmp/R_data/p3_dotplot_data.csv", check.names=FALSE)
dot$short.label <- factor(dot$short.label, levels=unique(dot$short.label[order(dot$cluster)]))
p3 <- ggplot(dot, aes(x=short.label, y=gene, size=pct_expressing, color=mean_expression)) +
  geom_point() +
  scale_color_gradient(low="white", high="#B2182B", name="Mean\nexpression",
                       breaks=pretty(range(dot$mean_expression), 3)) +
  scale_size_continuous(range=c(1,8), name="% expressing", breaks=c(25,50,75,100)) +
  labs(x="Cluster", y="Gene") +
  theme_classic(base_size=13) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
        legend.position="right",
        legend.title=element_text(size=10, face="bold"),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        panel.border=element_rect(fill=NA, color="grey30", linewidth=0.5))
ggsave(file.path(fig,"Fig1C_dotplot.png"), p3, width=10, height=7, dpi=300)

# Fig1D: cluster composition, 图例精细
comp <- d %>% group_by(group, short.label) %>% summarise(n=n(), .groups="drop") %>%
  group_by(group) %>% mutate(pct=n/sum(n)*100)
p4 <- ggplot(comp, aes(x=group, y=pct, fill=short.label)) +
  geom_col(color="white", linewidth=0.3, width=0.7) +
  scale_fill_igv(name="Cluster") +
  labs(x="Group", y="Percentage of cells (%)") +
  theme_classic(base_size=13) +
  theme(legend.title=element_text(size=10, face="bold"),
        legend.text=element_text(size=9),
        legend.key.size=unit(0.5,"cm"),
        panel.border=element_rect(fill=NA, color="grey30", linewidth=0.5))
ggsave(file.path(fig,"Fig1D_cluster_composition.png"), p4, width=7.5, height=6, dpi=300)

cat("=== Fig1A-D 完成(图例精细化) ===\n")
