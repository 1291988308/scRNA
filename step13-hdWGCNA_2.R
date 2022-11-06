# 推荐本地安装
#devtools::install_local('/home/data/t040243/hdWGCNA/hdWGNCA.zip')
#BiocManager::install('harmony',update=F,ask=F)
# 正式安装
#library(devtools)
#devtools::install_github('smorabit/hdWGCNA', ref='dev')

library(hdWGCNA)
#加载单细胞分析包
library(Seurat)
#加载作图包
library(tidyverse)
library(cowplot)
library(patchwork)
#加载共表达网络分析包
library(WGCNA)
gc()

# 工作目录
setwd("/home/data/t040243/CopyOfjiangshijiu")
dir.create("hdWGCNA")
setwd("/home/data/t040243/CopyOfjiangshijiu/hdWGCNA")
#设置随机种子
set.seed(12345)

# 读取单细胞数据集，读取自己的
load("/home/data/t040243/CopyOfjiangshijiu/Doublet/after_doublet_obj.combined.RData")
scRNA=obj.combined

if (F) {###提取亚群分析
  
  # 提上皮亚群,重新降维聚类
scRNA=subset(scRNA,celltype %in% 'ILCdc')
scRNA<- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)
scRNA<- RunPCA(scRNA, features = VariableFeatures(scRNA))
DimPlot(scRNA, reduction = "pca", group.by = "orig.ident")
ElbowPlot(scRNA)
scRNA <- FindNeighbors(scRNA, dims = 1:5)
scRNA <- FindClusters(scRNA)
scRNA <- RunUMAP(scRNA, dims = 1:5)
# 瞄一眼，小亚群
DimPlot(scRNA, label=TRUE) 
}


#过滤出至少在5%的细胞中表达的基因
scRNA <- SetupForWGCNA(
  scRNA,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "Bio_com" # the name of the hdWGCNA experiment
)


#构建metacells!!这一步非常重要，WGCNA对数据的稀疏性非常敏感，与普通转录组的WGCNA分析相比
# 单细胞的稀疏矩阵的解决方法是WGCNA流程的问题核心

# construct metacells  in each group
scRNA<- MetacellsByGroups(
  seurat_obj = scRNA,k=20,
  max_shared = 10,
  # group.by一般关注的是组织类型和细胞类型!这边组织类型是orig.ident，CT正常，PR疾病
  group.by =c("celltype"), ### c("celltype",'orig.ident'), # 也可以选择别的groupby
  ident.group = 'celltype' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
scRNA <- NormalizeMetacells(scRNA)
metacell_obj <- GetMetacellObject(scRNA)


#转置表达矩阵
# 安全起见，另起一个对象，以角质细胞细胞为例
seurat_obj  <- SetDatExpr(
  scRNA,
  group_name = c("ILCdc","ILC1","ILC2a","ILC2b","ILCt"), # 选择感兴趣恶的细胞类群！
  group.by='celltype' # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
)


#选择softpower
seurat_obj <- TestSoftPowers(
  seurat_obj,
  setDatExpr = FALSE, # 这边不用转置了，前面已转置
)


# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
pdf("softpower.pdf",width = 12,height = 8)
wrap_plots(plot_list, ncol=2)
dev.off()


#查看powerTable
power_table <- GetPowerTable(seurat_obj)
head(power_table)


#构建共表达网络
softpower=6  # 根据自己的改，选5也没问题
# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power=softpower,
  group.by='celltype', group_name='ILC',setDatExpr = F)


#可视化WGCNA网络
pdf("hdWGCNA Dendrogram.pdf",width = 12,height = 8)
PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')
dev.off()

#(可选)获取TOM矩阵，可用来进行其他高级分析
TOM <- GetTOM(seurat_obj)


#计算模块协调特征
#记得scale一下 or else harmony throws an error:
seurat_obj <- Seurat::ScaleData(
  seurat_obj,
  features = GetWGCNAGenes(seurat_obj),
  
)
# 计算ME，根据组织类型分组
# harmony必须biocManager安装，不可以用github安装！！！
library(harmony)

seurat_obj <- ModuleEigengenes(
  seurat_obj,
  #group.by.vars="orig.ident" #harmony对象     #####是否需要改为celltype？
  group.by.vars="celltype" #harmony对象     #####是否需要改为celltype？
)


seurat_obj <- ModuleConnectivity(seurat_obj)
# plot genes ranked by kME for each module
#可视化每个模块中，按照kME打分的基因PlotKMEs(seurat_obj, ncol=5)



# 获取hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 25)
head(hub_df)
write.csv(hub_df,file = "hub_df.csv")


p1 <- PlotKMEs(seurat_obj,text_size =6,ncol = 5)
p1
ggsave("MEs.pdf",p1,width = 18,height = 8,dpi = 300)

#记得保存上面hdWGNCA关键分析过程！！！
saveRDS(seurat_obj, file='hdWGCNA_object.rds')

dev.off()


####------------一些可视化-----------------------
## 模块间的相关性
library(igraph)
library(qgraph)
# 载入保存的

seurat_obj=readRDS("/home/data/t040243/CopyOfjiangshijiu/hdWGCNA/hdWGCNA_object.rds")

# 画模块间相关性图
pdf("correlogram.pdf",width = 12,height = 8)
ModuleCorrelogram(seurat_obj, sig.level = 0.001, pch.cex=2)
dev.off()
# 由于识别到了每个模块的hub基因，可以去计算hub基因打分
# compute gene scoring for the top 25 hub genes by kME for each module
# (方法一)with Seurat method
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='Seurat'
)

# compute gene scoring for the top 25 hub genes by kME for each module
# (方法二)with UCell method #推荐这种方法
# 由于Ucell刚刚更新，所以4.1.x的同学请用本地安装,依赖包自行安装
#devtools::install_local("/home/data/t040243/hdWGCNA/UCell-1.3.zip")
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)

# featureplot
# 瞄一眼
DimPlot(scRNA, label=TRUE,split.by = 'celltype') 

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE ,# order so the points with highest hMEs are on top
)

# stitch together with patchwork
p <- wrap_plots(plot_list,ncol=5,widths = 12,heights = 8)
p
ggsave("hMEs.pdf",p,width = 12,height = 4,dpi = 300)

### dotplot
# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods)

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
p
ggsave("dot.pdf",p,width = 12,height = 8,dpi = 300)

