#加载seurat数据和包
# single-cell analysis package
setwd("/home/data/t040243/CopyOfjiangshijiu/hdWGCNA")
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# gene enrichment packages
#install.packages('enrichR')
library(enrichR)

#BiocManager::install('GeneOverlap',update = F,ask = F)
library(GeneOverlap)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# load the Zhou et al snRNA-seq dataset
seurat_obj <- readRDS('hdWGCNA_object.rds')

#GO富集分析
# enrichr databases to test
dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021')

# 富集分析，会逐个模块分析
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test
)


#记得保存上面hdWGNCA关键分析过程！！！
saveRDS(seurat_obj, file='hdWGCNA_object.rds')

# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_obj)
write.csv(enrich_df,file = "enrich_df.csv")

# make GO term plots作图，在文件夹下生成！
EnrichrBarPlot(
  seurat_obj,
  outdir = "enrichr_plots", # name of output directory
  n_terms = 5, # number of enriched terms to show (sometimes more show if there are ties!!!)
  plot_size = c(5,7), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)


#气泡图
# GO_Biological_Process_2021
pdf("EnrichrDotpoot_BP.pdf",width = 12,height = 15)
EnrichrDotPlot(
  seurat_obj,
 # mods = c("turquoise","black"), # use all modules (this is the default behavior)
  database = "GO_Biological_Process_2021", # this has to be one of the lists we used above!!!
  n_terms=3
  # number of terms for each module
)
dev.off()

pdf("EnrichrDotpoot_BP_4mods.pdf",width = 12,height = 8)
EnrichrDotPlot(
  seurat_obj,
  mods = c("blue","turquoise","pink","green"), # use all modules (this is the default behavior)
  database = "GO_Biological_Process_2021", # this has to be one of the lists we used above!!!
  n_terms=3
  # number of terms for each module
)
dev.off()


#气泡图
# GO_Cellular_Component_2021
pdf("EnrichrDotpoot_CC.pdf",width = 12,height = 15)
EnrichrDotPlot(
  seurat_obj,
 # mods = c("turquoise","black",'blue'), # use all modules (this is the default behavior)
  database = "GO_Cellular_Component_2021", # this has to be one of the lists we used above!!!
  n_terms=3
  # number of terms for each module
)
dev.off()


pdf("EnrichrDotpoot_CC_4mods.pdf",width = 12,height = 8)
EnrichrDotPlot(
  seurat_obj,
  mods = c("blue","turquoise","pink","green"), # use all modules (this is the default behavior)
  database = "GO_Cellular_Component_2021", # this has to be one of the lists we used above!!!
  n_terms=3
  # number of terms for each module
)
dev.off()


#气泡图
# GO_Biological_Process_2021
pdf("EnrichrDotpoot_MF.pdf",width = 12,height = 15)
EnrichrDotPlot(
  seurat_obj,
 # mods = c("turquoise","black",'red'), # use all modules (this is the default behavior)
  database = "GO_Molecular_Function_2021", # this has to be one of the lists we used above!!!
  n_terms=3
  # number of terms for each module
)
dev.off()

pdf("EnrichrDotpoot_MF_4mods.pdf",width = 12,height = 8)
EnrichrDotPlot(
  seurat_obj,
  mods = c("blue","turquoise","pink","green"), # use all modules (this is the default behavior)
  database = "GO_Molecular_Function_2021", # this has to be one of the lists we used above!!!
  n_terms=3
  # number of terms for each module
)
dev.off()
#差异基因重叠分析
## 这个分析帮助我们看到，哪些模块可能是相似的
# compute cell-type marker genes with Seurat:
# 常规方法计算差异基因/特征基因
Idents(seurat_obj) <- seurat_obj$celltype
markers <- Seurat::FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  logfc.threshold=1
)

# compute marker gene overlaps
overlap_df <- OverlapModulesDEGs(
  seurat_obj,
  deg_df = markers,
  fc_cutoff = 1 # log fold change cutoff for overlap analysis
)
write.csv(overlap_df,file = "overlapModulesDEGs.csv")
#条形图
# overlap barplot, produces a plot for each cell type
plot_list <- OverlapBarPlot(overlap_df,label_size = 5)

# stitch plots with patchwork
p1 <- wrap_plots(plot_list, ncol=5)
p1
ggsave("overlap barplot_ncol=5.pdf",p1,width = 24,height = 6,dpi = 300)
#气泡图
# plot odds ratio of the overlap as a dot plot
p2 <- OverlapDotPlot(
  overlap_df,
  plot_var = 'odds_ratio') +
  ggtitle('Overlap of modules & cell-type markers')
ggsave("overlap dotplot.pdf",p2,width = 6,height = 4,dpi = 300)

#----------------------------
#网络可视化
# network analysis & visualization package:
# network analysis & visualization package:
library(igraph)

#可视化每个模块的网络图
ModuleNetworkPlot(seurat_obj)

#组合网络图，在文件夹下生成
# hubgene network
pdf("Hubgenenet.pdf",width = 8,height = 6)
 HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 3, n_other=6,
  edge_prop = 0.75,
  mods = "all",
  hub.vertex.size = 8,
  other.vertex.size = 1,
  vertex.label.cex = 0.75)
dev.off()



#UMAP可视化
## 利用hub基因，重新UMAP，如此可以获得分群明显的图
g <- HubGeneNetworkPlot(seurat_obj,  return_graph=TRUE)

seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)


# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)
write.csv(umap_df,"umap_df.csv")

# plot with ggplot
pdf("ggplot_umap_df.pdf",width = 12,height = 8)
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()
dev.off()

pdf("Moduleumap.pdf",width = 12,height = 8)
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=6 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
  )
dev.off()
#####################################################################
#####################################################################
#监督UMAP
g <- ModuleUMAPPlot(seurat_obj,  return_graph=TRUE)
# run supervised UMAP:
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10,
  n_neighbors=15,
  min_dist=0.1,
  supervised=TRUE,
  target_weight=0.5
)

# get the hub gene UMAP table from the seurat object
supervised_umap_df <- GetModuleUMAP(seurat_obj)
write.csv(supervised_umap_df,"supervised_umap_df.csv")

# plot with ggplot
pdf("ggplot_supervised_umap_df.pdf",width = 12,height = 8)
ggplot(supervised_umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()
dev.off()

pdf("supervised_Moduleumap.pdf",width = 12,height = 8)
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=6 ,# how many hub genes to plot per module?
  keep_grey_edges=T,
  vertex.label.cex = 1,
  
)
dev.off()



################################################################################
################################################################################
dir.create("top6hubgene")
setwd("top6hubgene")
hub_df_top6 <- GetHubGenes(seurat_obj, n_hubs = 6)
head(hub_df_top6)
write.csv(hub_df_top6,file = "hub_df_top6.csv")

for( i in unique(hub_df_top6$module) ){
  markers_genes =  hub_df_top6[hub_df_top6$module==i,]$gene_name
  print(markers_genes)
  FeaturePlot(seurat_obj, features=markers_genes,reduction = "umap",cols = c("lightgrey","red2"))
  ggsave(filename=paste0('FeaturePlot_',i,'.pdf'),width = 12,height = 8)
}
