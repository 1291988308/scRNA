#细胞鉴定
#鉴定marker基因
setwd("/home/data/t040243/CopyOfjiangshijiu")
dir.create("cluster")
setwd("/home/data/t040243/CopyOfjiangshijiu/cluster")
sel.clust = "integrated_snn_res.0.1"
obj.combined <- SetIdent(obj.combined, value = sel.clust)
sample.markers <- FindAllMarkers(obj.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- sample.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#保存每个cluster top10的marker基因
top10_table=unstack(top10, gene ~ cluster)
names(top10_table)=gsub("X","cluster",names(top10_table))
write.csv(file="top10_marker_genes.csv",top10_table,row.names=F)
write.csv(file="sample_marker_genes.csv",sample.markers,row.names=F)
#save(obj.combined,sample.markers,top10,top10_table,file = "cluster.RData")
#细胞亚群注释——命名
if(T){
new.cluster.ids <- c("ILC2a","ILCdc","ILC2b","ILC1",
                    "ILCt","ILCdc",  "ILCdc")
names(new.cluster.ids) <- levels(obj.combined)
obj.combined <- RenameIdents(obj.combined, new.cluster.ids)
table(obj.combined@active.ident)
}

obj.combined$celltype <- Idents(obj.combined) ##保存每个细胞的细胞类型
#细胞亚群注释后的UMAP
library(RColorBrewer)   #配色
cols <- c(brewer.pal(8, "Dark2")) #配色
pdf('cluster_annotations.pdf',width = 12,height =8)                    
DimPlot(obj.combined, reduction = "umap", label = TRUE,
        label.size = 12,label.box = F,cols = cols,#shape.by = "group",
        group.by = "celltype", pt.size = 1.5,repel=T) +theme_bw()+
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size =30),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30),
        title = element_text(size = 30))+
  guides(color = guide_legend(override.aes = list(size = 6)))  ###控制图例的大小
dev.off()

pdf('cluster_annotations_group.pdf',width = 36,height = 12)
DimPlot(obj.combined, reduction = "umap", label = TRUE,
        label.size = 12,label.box = F,cols = cols,shape.by = "group",
        split.by = "group", pt.size = 3.5,repel=T,) +theme_bw()+
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size =30),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30),
        title = element_text(size = 30),
        strip.text = element_text(face = "bold", color = "#7570B3",
                                  hjust = 0, size = 50),
        strip.background = element_rect(fill = "#E6AB02", linetype = "dotted"))+
  guides(color = guide_legend(override.aes = list(size = 12)))
dev.off()

save(obj.combined,file = "rename_cluster.RData")

###########################################################################
#############################################################################

################################################################


