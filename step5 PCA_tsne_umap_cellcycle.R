#对整合之后的数据进行降维
#主成分分析
setwd("/home/data/t040243/CopyOfjiangshijiu")
#save(obj.combined,file = "obj.combined.RData")
saveRDS(obj.combined, "obj.combined.rds")
dir.create("PCA")
setwd("/home/data/t040243/CopyOfjiangshijiu/PCA")
DefaultAssay(obj.combined) <- "integrated"
obj.combined@meta.data[["group"]] <- factor(obj.combined@meta.data[["group"]],levels=c("Sham","MIRI","MI"))
######给分组指定因子水平Sham>MIRI>MI

obj.combined <- RunPCA(obj.combined, npcs = 50, verbose = FALSE)
#可视化PCA降维之后的结果
pdf('2_PCA.pdf', width = 10, height = 6)
DimPlot(obj.combined, reduction = "pca")
dev.off()
#热图展示前5个主成分
pdf(file="2_PCtop5_Heatmap.pdf",width=14,height=20)
DimHeatmap(obj.combined, dims = 1:5, balanced = TRUE)
dev.off()
#热图展示第6-10个主成分
pdf(file="2_PC6-10_Heatmap.pdf",width=14,height=20)
DimHeatmap(obj.combined, dims = 1:5, balanced = TRUE)
dev.off()
#确定用几个主成分做后续分析
#Elbowplot
pdf('2_ElbowPlot.pdf', width = 10, height = 6)
ElbowPlot(obj.combined, ndims = 50)
dev.off()
#对整合之后的数据进行聚类
obj.combined <- RunUMAP(obj.combined, reduction = "pca", dims = 1:30)
obj.combined <- RunTSNE(obj.combined, reduction = "pca", dims = 1:30)
 #SWNE Visualizing single-cell RNA-seq datasets with Similarity Weighted Nonnegative Embedding (tSNE可准确的捕获数据集的局部结构，但是会扭曲数据集的全局结构，比如簇之间的距离)
obj.combined <- FindNeighbors(obj.combined, reduction = "pca", dims = 1:30)
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
obj.combined <- FindClusters(obj.combined,resolution = res)
} 
#设置不同的分辨率，观察分群效果(选择哪一个？)
#可视化不同分辨率，分群效果
apply(obj.combined@meta.data[,grep("integrated_snn",colnames(obj.combined@meta.data))],2,table)
pdf('2_umap_resolution_high.pdf', width = 18)
plot_grid(ncol = 3, DimPlot(obj.combined, reduction = "umap", group.by = "integrated_snn_res.0.8") + 
                ggtitle("louvain_0.8"), DimPlot(obj.combined, reduction = "umap", group.by = "integrated_snn_res.1") + 
                 ggtitle("louvain_1"), DimPlot(obj.combined, reduction = "umap", group.by = "integrated_snn_res.0.3") + 
                  ggtitle("louvain_0.3"))
dev.off()
pdf('2_umap_resolution_low.pdf', width = 18)
plot_grid(ncol = 3, DimPlot(obj.combined, reduction = "umap", group.by = "integrated_snn_res.0.01") + 
                 ggtitle("louvain_0.01"), DimPlot(obj.combined, reduction = "umap", group.by = "integrated_snn_res.0.1") + 
                 ggtitle("louvain_0.1"), DimPlot(obj.combined, reduction = "umap", group.by = "integrated_snn_res.0.2") + 
                 ggtitle("louvain_0.2"))
dev.off()
pdf('2_Tree_diff_resolution.pdf', width = 10,height = 10)
clustree(obj.combined@meta.data, prefix = "integrated_snn_res.",layout = "sugiyama")
dev.off()
#接下来分析，按照分辨率为0.1进行
sel.clust = "integrated_snn_res.0.1"
obj.combined <- SetIdent(obj.combined, value = sel.clust)
table(obj.combined@active.ident)
#绘制UMAP聚类图
pdf('3_UMAP_multi_samples_combined.pdf',width = 11,height = 6)
p1 <- DimPlot(obj.combined, reduction = "umap", group.by = "group")
p2 <- DimPlot(obj.combined, reduction = "umap", label = TRUE)
p1 + p2
dev.off()
#分别展示两组的UMAP聚类图
pdf('3_UMAP_multi_samples_split.pdf',width = 10,height = 6)
DimPlot(obj.combined, reduction = "umap", split.by = "group",label = TRUE)
dev.off()

#查看细胞周期影响
pdf(file="3_Pcna_Mki67_expression.pdf",width=12)
FeaturePlot(obj.combined,features = c("Pcna","Mki67"),reduction = "umap")  #s期基因PCNA，G2/M期基因MKI67均高表达，说明细胞均处于细胞周期
dev.off()
#计算细胞周期分值，判断这些细胞分别处于什么期
cc_gene=unlist(cc.genes)  #cc.genes 为细胞周期相关的基因，为list，有两个元素s.genes，g2m.genes
s.genes=cc.genes$s.genes
s.genes<-tolower(s.genes)
s.genes<-capitalize(s.genes)
g2m.genes=cc.genes$g2m.genes
g2m.genes<-tolower(g2m.genes)
g2m.genes<-capitalize(g2m.genes)
obj.combined <- CellCycleScoring(obj.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
pdf(file="3_cellcycle_score.pdf",width=8)
obj.combined@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+theme_minimal()
dev.off()
#去除细胞周期影响
obj.combined1 <- ScaleData(obj.combined, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(obj.combined))
save(obj.combined,obj.combined1,file = "PCA_cellcycle.RData")

