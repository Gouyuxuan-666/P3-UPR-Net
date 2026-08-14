suppressMessages({library(ggplot2)})
fig <- "/root/autodl-tmp/R_figures/Fig1_single"
d <- data.frame(transition=c("HC → DS","DS → DR"), cost=c(0.153,0.262), x=c(1,2))
p <- ggplot(d, aes(x=transition, y=cost)) +
  geom_col(fill=c("#4575B4","#B2182B"), width=0.5, alpha=0.85) +
  geom_text(aes(label=sprintf("%.3f", cost)), vjust=-0.5, size=5, fontface="bold") +
  annotate("segment", x=1, xend=2, y=0.29, yend=0.29, arrow=arrow(length=unit(0.2,"cm")), color="grey40", linewidth=0.6) +
  annotate("text", x=1.5, y=0.32, label="disease progression", size=4, color="grey30") +
  labs(x="", y="Optimal transport cost") +
  ylim(0, 0.35) +
  theme_classic(base_size=13) +
  theme(axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=10),
        axis.title.y=element_text(size=12),
        panel.border=element_rect(fill=NA, color="grey30", linewidth=0.5))
ggsave(file.path(fig,"Fig1E_moscot_transition.png"), p, width=6.5, height=5, dpi=300)
cat("Fig1E 完成\n")
