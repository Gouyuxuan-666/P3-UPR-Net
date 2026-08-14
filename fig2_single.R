suppressMessages({library(ggraph); library(igraph); library(ggalt); library(ComplexHeatmap); library(circlize); library(ggplot2); library(dplyr); library(tidyr); library(tibble)})

fig <- "/root/autodl-tmp/R_figures/Fig2_single"
dir.create(fig, recursive=TRUE, showWarnings=FALSE)

Nred <- "#B2182B"; Nblue <- "#2166AC"; Norange <- "#F18F01"; Npurple <- "#7B3294"; Ngrey <- "#BDBDBD"

sp110_dr <- c(ATF4=-0.36, MAPK8=0.21, BAX=0.14, HSPA5=-0.28, DDIT3=-0.22, XBP1=-0.18, ERN1=0.12, CASP3=0.08)
sp110_ds <- c(ATF6=0.40, BECN1=0.31, HSPA5=-0.30, XBP1=0.22, MAP1LC3B=0.25, ATF4=0.18, DDIT3=0.15, SQSTM1=-0.20)
notears <- data.frame(from=c("CASP3","HSPA5","EIF2S1"), to=c("HSPA5","XBP1","HSPA5"), weight=c(1.07,0.84,0.77))
all_targets <- sort(unique(c(names(sp110_dr), names(sp110_ds))))

# Fig2A: causal network
sp_edge <- full_join(tibble(to=names(sp110_dr), dr=as.numeric(sp110_dr)),
                     tibble(to=names(sp110_ds), ds=as.numeric(sp110_ds)), by="to") %>%
  mutate(dr=replace_na(dr,0), ds=replace_na(ds,0), both=dr!=0 & ds!=0,
         effect=case_when(both ~ (dr+ds)/2, dr!=0 ~ dr, TRUE ~ ds),
         abs_eff=abs(effect),
         type=case_when(both & sign(dr)!=sign(ds) ~ "Sign-flip", effect>0 ~ "Activation", TRUE ~ "Repression"),
         from="SP110") %>% select(from,to,effect,abs_eff,type)
nt_edge <- notears %>% transmute(from,to,effect=weight,abs_eff=abs(weight),type="NOTEARS")
edges <- bind_rows(sp_edge, nt_edge)
edges$disp_w <- ifelse(edges$type=="NOTEARS", edges$abs_eff*0.28, edges$abs_eff)
node_ids <- unique(c(edges$from, edges$to))
nodes <- tibble(name=node_ids, is_hub=node_ids=="SP110",
                role=case_when(node_ids=="SP110"~"Regulator", node_ids %in% notears$from~"Mediator", TRUE~"Target"))
g <- graph_from_data_frame(edges, directed=TRUE, vertices=nodes)
edge_pal <- c("Activation"=Nred,"Repression"=Nblue,"Sign-flip"=Npurple,"NOTEARS"=Norange)
set.seed(42)
p2A <- ggraph(g, layout="kk") +
  geom_edge_arc(aes(edge_color=type, edge_width=disp_w, edge_linetype=type), strength=0.16,
                arrow=arrow(type="closed", length=unit(2.4,"mm")),
                start_cap=circle(4.5,"mm"), end_cap=circle(5.5,"mm"), alpha=0.85) +
  scale_edge_color_manual(values=edge_pal, name="Edge type") +
  scale_edge_width_continuous(range=c(0.35,1.8), guide="none") +
  scale_edge_linetype_manual(values=c("Activation"="solid","Repression"="solid","Sign-flip"="solid","NOTEARS"="dashed"), guide="none") +
  geom_node_point(aes(filter=!is_hub, fill=role), shape=21, size=6, colour="grey25", stroke=0.4) +
  geom_node_point(aes(filter=is_hub, fill=role), shape=21, size=13, colour="grey15", stroke=0.6) +
  scale_fill_manual(values=c(Regulator=Nred, Mediator=Norange, Target="white"), name="Node role",
                    guide=guide_legend(override.aes=list(size=5))) +
  geom_node_text(aes(filter=!is_hub, label=name), fontface="italic", size=3.6, repel=TRUE, max.overlaps=Inf, box.padding=0.30, seed=42, colour="grey15") +
  geom_node_text(aes(filter=is_hub, label=name), fontface="bold", size=4.5, repel=TRUE, box.padding=0.55, seed=42, colour="grey10") +
  theme_graph(background="white", base_family="sans") +
  theme(legend.position="right", legend.title=element_text(face="bold", size=10), legend.text=element_text(size=9))
ggsave(file.path(fig,"Fig2A_causal_network.png"), p2A, width=8.8, height=6.6, dpi=300)

# Fig2B: DR vs DS dumbbell
db <- tibble(gene=all_targets, dr=ifelse(all_targets %in% names(sp110_dr), sp110_dr[all_targets],0),
             ds=ifelse(all_targets %in% names(sp110_ds), sp110_ds[all_targets],0)) %>%
  mutate(abs_diff=abs(dr-ds)) %>% arrange(abs_diff)
db$gene <- factor(db$gene, levels=db$gene)
p2B <- ggplot(db, aes(y=gene, x=dr, xend=ds)) +
  geom_vline(xintercept=0, colour="grey60", linewidth=0.35) +
  ggalt::geom_dumbbell(colour=Ngrey, colour_x=Nblue, colour_xend=Nred, size=0.9, size_x=3.6, size_xend=3.6, dot_guide=FALSE) +
  geom_text(aes(x=dr, label=sprintf("%.2f",dr)), hjust=1.25, size=3, colour=Nblue, fontface="bold") +
  geom_text(aes(x=ds, label=sprintf("%.2f",ds)), hjust=-0.25, size=3, colour=Nred, fontface="bold") +
  scale_x_continuous(expand=expansion(mult=c(0.18,0.18))) +
  labs(x="Jacobian causal effect", y=NULL) +
  theme_minimal(base_size=12) +
  theme(panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
        axis.text.y=element_text(face="italic", size=11), axis.text.x=element_text(size=10),
        panel.border=element_rect(fill=NA, color="grey30", linewidth=0.5))
ggsave(file.path(fig,"Fig2B_DR_vs_DS_divergent.png"), p2B, width=7.6, height=5.6, dpi=300)

# Fig2C: causal matrix heatmap
reg_rows <- c("SP110_DR","SP110_DS","HSPA5","CASP3","EIF2S1")
mat <- matrix(0, nrow=length(reg_rows), ncol=length(all_targets), dimnames=list(reg_rows, all_targets))
for(g in names(sp110_dr)) mat["SP110_DR",g] <- sp110_dr[[g]]
for(g in names(sp110_ds)) mat["SP110_DS",g] <- sp110_ds[[g]]
for(i in seq_len(nrow(notears))){ r<-notears$from[i]; ct<-notears$to[i]; if(r %in% reg_rows && ct %in% all_targets) mat[r,ct]<-notears$weight[i] }
mx <- ceiling(max(abs(mat), na.rm=TRUE)*10)/10
col_fun <- colorRamp2(c(-mx,0,mx), c(Nblue,"white",Nred))
ht <- Heatmap(mat, name="Causal\neffect", col=col_fun, cluster_rows=TRUE, cluster_columns=TRUE,
  clustering_distance_rows="euclidean", clustering_distance_columns="euclidean",
  clustering_method_rows="complete", clustering_method_columns="complete",
  show_row_dend=TRUE, show_column_dend=TRUE, row_dend_width=unit(9,"mm"), column_dend_height=unit(9,"mm"),
  row_names_side="left", column_names_side="bottom", column_names_rot=45,
  column_names_gp=gpar(fontface="italic", fontsize=10), row_names_gp=gpar(fontface="bold", fontsize=10.5),
  rect_gp=gpar(col="white", lwd=1),
  heatmap_legend_param=list(at=pretty(c(-mx,mx),n=5), labels=sprintf("%.2f",pretty(c(-mx,mx),n=5)),
                            legend_height=unit(3.6,"cm"), title_gp=gpar(fontface="bold", fontsize=10)),
  cell_fun=function(j,i,x,y,w,h,fill){ v<-mat[i,j]; if(abs(v)>=0.15) grid.text(sprintf("%.2f",v),x,y,gp=gpar(fontsize=8,col=ifelse(abs(v)>mx*0.55,"white","grey15"))) })
png(file.path(fig,"Fig2C_UPR_causal_matrix.png"), width=9.6, height=4.8, units="in", res=300)
draw(ht, heatmap_legend_side="right", annotation_legend_side="right", merge_legend=TRUE, padding=unit(c(4,4,4,4),"mm"))
dev.off()

# Fig2D/E: SP110 DR/DS lollipop (split)
make_lolli <- function(vec, x_range) {
  d <- tibble(gene=names(vec), effect=as.numeric(vec)) %>% arrange(effect)
  d$gene <- factor(d$gene, levels=d$gene); d$dir <- ifelse(d$effect>=0,"Activation","Repression")
  ggplot(d, aes(x=effect, y=gene)) +
    geom_segment(aes(x=0, xend=effect, yend=gene, colour=dir), linewidth=0.9) +
    geom_point(aes(fill=dir), shape=21, size=4.6, colour="grey20", stroke=0.4) +
    geom_text(aes(label=sprintf("%.2f",effect), hjust=ifelse(effect>=0,-0.25,1.25)), fontface="bold", size=3.2, colour="grey15") +
    geom_vline(xintercept=0, colour="grey30", linewidth=0.4) +
    scale_colour_manual(values=c(Activation=Nred, Repression=Nblue), guide="none") +
    scale_fill_manual(values=c(Activation=Nred, Repression=Nblue), name=NULL) +
    scale_x_continuous(limits=x_range, expand=expansion(mult=c(0.02,0.02))) +
    labs(x="Jacobian causal effect", y=NULL) +
    theme_classic(base_size=12) +
    theme(panel.grid.major.x=element_line(colour="grey90", linewidth=0.25),
          axis.text.y=element_text(face="italic", size=11), axis.text.x=element_text(size=10),
          legend.title=element_text(size=10, face="bold"), legend.text=element_text(size=9),
          legend.position="bottom", panel.border=element_rect(fill=NA, color="grey30", linewidth=0.5))
}
xr <- range(c(sp110_dr, sp110_ds)); xr <- xr + c(-1,1)*0.22*diff(xr)
ggsave(file.path(fig,"Fig2D_SP110_DR_effects.png"), make_lolli(sp110_dr, xr), width=6, height=5.2, dpi=300)
ggsave(file.path(fig,"Fig2E_SP110_DS_effects.png"), make_lolli(sp110_ds, xr), width=6, height=5.2, dpi=300)

cat("=== Fig2A-E 完成 ===\n")
