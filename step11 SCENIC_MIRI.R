#基因调控网络 SCENIC 和DoRothEA
#DoRothEA
#获取包自带数据库
set.seed(123)  #设置随机数种子，使结果可重复
setwd("/home/shpc_100537/MIRI_SCENIC")
dir.create("/home/shpc_100537/MIRI_SCENIC/DoRothEA")
setwd("/home/shpc_100537/MIRI_SCENIC/DoRothEA")
load("/home/shpc_100537/CopyOfjiangshijiu/Doublet/after_doublet_obj.combined.RData")


dorothea_regulon_mouse <- get(data("dorothea_mm", package = "dorothea"))
#推断regulons基于levelA,B,C
regulon <- dorothea_regulon_mouse %>%
  dplyr::filter(confidence %in% c("A","B","C"))
#viper得分矩阵获取
obj.combined <- run_viper(obj.combined, regulon,
                          options = list(method = "scale", minsize = 4, 
                                         eset.filter = FALSE, cores = 1,  
                                         verbose = FALSE))
obj.combined@assays$dorothea@data[1:4,1:4]
dim(obj.combined@assays$dorothea@data)
#根据细胞亚群汇总
viper_scores_df <-GetAssayData(obj.combined,slot="data",assay="dorothea")%>%data.frame(check.names = F) %>%t()
CellsClusters <- data.frame(cell = names(Idents(obj.combined)), 
                            cell_type=as.character(Idents(obj.combined)),
                            check.names = F)
viper_scores_clusters <- viper_scores_df  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)
summarized_viper_scores <- viper_scores_clusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))
head(summarized_viper_scores)
#各个单细胞亚群特异性的转录因子
DefaultAssay(object = obj.combined) <- "dorothea"
markers <- FindAllMarkers(object = obj.combined, only.pos = TRUE,min.pct = 0.25, thresh.use = 0.25)
pro='dorothea-markers-for-mouse'
write.csv(markers,file=paste0(pro,'_.markers.csv'))
#挑选前40TFs可视化
highly_variable_tfs <- summarized_viper_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(200, var) %>% ##40*5群
  distinct(tf)
highly_variable_tfs
#准备数据可视化
summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) 
palette_length = 100
my_color = colorRampPalette(c("Darkblue","white","red"))(palette_length)
colnames(summarized_viper_scores_df)
#pheatmap热图可视化
my_breaks <- c(seq(min(summarized_viper_scores_df), 0, 
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(summarized_viper_scores_df)/palette_length, 
                   max(summarized_viper_scores_df), 
                   length.out=floor(palette_length/2)))
pdf("DoRothEA_heatmap.pdf",width=12)
pheatmap(t(summarized_viper_scores_df),fontsize=14, 
         fontsize_row = 10, 
         color=my_color, breaks = my_breaks, 
         main = "DoRothEA (ABC)", angle_col = 45,
         treeheight_col = 0,  border_color = NA)
dev.off()






#SCENIC
#参考数据库和结果输出文件夹
#数据库下载https://resources.aertslab.org/cistarget/
#dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather")
#for(featherURL in dbFiles)
#{
#  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
#}

dir.create("../SCENIC")
dir.create("../SCENIC/int")
setwd("../SCENIC")
#准备细胞meta信息和表达矩阵（挑选疾病组的细胞亚群）
Idents(obj.combined) <- "group"
obj.combined <- subset(obj.combined,idents="MIRI")
Idents(obj.combined) <- "celltype"
cellInfo <- data.frame(obj.combined@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="integrated_snn_res.0.1")] <- "cluster"
colnames(cellInfo)[which(colnames(cellInfo)=="celltype")] <- "celltype"
cellInfo <- cellInfo[,c("sample","cluster","celltype")]
#随机抽取1000个细胞的数据
#subcell <- sample(colnames(obj.combined),1000)
combine.sub <- obj.combined
saveRDS(combine.sub, "combine.sub.rds")
#准备meta信息
exprMat <- as.matrix(combine.sub@assays$RNA@counts) #准备表达矩阵
#设置分析环境
mydbs <- c("mm9-500bp-upstream-7species.mc9nr.feather",
           "mm9-tss-centered-10kb-7species.mc9nr.feather")
names(mydbs) <- c("500bp", "10kb")
scenicOptions <- initializeScenic(org="mgi",
                                  nCores=14,
                                  dbs = mydbs,
                                  dbDir="cisTarget_databases")
saveRDS(scenicOptions, "int/scenicOptions.rds")
#转录调控网络推断
#基因过滤
genesKept <- geneFiltering(exprMat, scenicOptions, 
                           minCountsPerGene = 3 * 0.01 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ] 
#过滤标准是基因表达量之和>细胞数*3%，且在1%的细胞中表达
#计算相关性矩阵
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) #TF-Targets相关性分析
#共表达模块和regulon推断
runGenie3(exprMat_filtered_log, scenicOptions, nParts = 5)   #这一步需要根据自己的处理器核心进行，时间比较长
runSCENIC_1_coexNetwork2modules(scenicOptions)#推断共表达模块
runSCENIC_2_createRegulons(scenicOptions,coexMethods = "top5perTarget")  #推断转录调控网络（regulon）
#以上代码可增加参数coexMethod=c("w001", "w005", "top50", "top5perTarget", "top10perTarget", "top50perTarget"))
#默认6种方法的共表达网络都计算，可以少选几种方法以减少计算量
#regulon活性评分与可视化
exprMat_all <- as.matrix(obj.combined@assays$RNA@counts)
exprMat_all <- log2(exprMat_all+1)
runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat_all)
#活性评分转换为二进制
#使用shiny互动调整阈值
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_all)
savedSelections <- shiny::runApp(aucellApp)
#保存调整后的阈值
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_all)
#分析结果可视化
#利用Seurat可视化AUC
#转录因子富集结果
scenicOptions=readRDS(file="int/scenicOptions.Rds")
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes") 
as.data.frame(sort(table(motifEnrichment_selfMotifs_wGenes$highlightedTFs),decreasing = T)) #每个基因的motif数量
#可视化某个基因的motif序列特征
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Cebpb"]
viewMotifs(tableSubset) 
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC #整合原始regulonAUC矩阵
obj.combined.auc <- AddMetaData(obj.combined, AUCmatrix)
obj.combined.auc@assays$integrated <- NULL
saveRDS(obj.combined.auc,'obj.combined.auc.rds')
#整合二进制regulonAUC矩阵
BINmatrix <- readRDS("int/4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
obj.combined.bin <- AddMetaData(obj.combined, BINmatrix)
obj.combined.bin@assays$integrated <- NULL
saveRDS(obj.combined.bin, 'obj.combined.bin.rds')





###FeaturePlot&idgePlot&VlnPlot
dir.create("3group_plot")
setwd("3group_plot")

for (i in 1:length(RegulonName_BIN)) {
  pdf(paste0(RegulonName_BIN[i],".pdf"),  width=18 ,height=8)
  p1 = FeaturePlot(obj.combined.auc, features=RegulonName_BIN[i], label=T, reduction = 'umap')
  p2 = FeaturePlot(obj.combined.bin, features=RegulonName_BIN[i], label=T, reduction = 'umap')
  p3 = DimPlot(obj.combined, reduction = 'umap', group.by = "celltype", label=T)
  plotc = p1|p2|p3
  print(plotc)
  dev.off()
}

dir.create("../3group_plot_auc")
setwd("../3group_plot_auc")
for (i in 1:length(RegulonName_AUC)) {
  pdf(paste0(RegulonName_AUC[i],".pdf"),  width=18 ,height=8)
  p1 = FeaturePlot(obj.combined.auc, features=RegulonName_AUC[i], label=T, reduction = 'umap')
  # p2 = FeaturePlot(obj.combined.bin, features=RegulonName_BIN[i], label=T, reduction = 'umap')
  p3 = DimPlot(obj.combined, reduction = 'umap', group.by = "celltype", label=T)
  plotc = p1|p3
  print(plotc)
  dev.off()
}


dir.create("../3group_plot_bin")
setwd("../3group_plot_bin")
for (i in 1:length(RegulonName_BIN)) {
  pdf(paste0(RegulonName_BIN[i],".pdf"),  width=18 ,height=8)
  #p1 = FeaturePlot(obj.combined.auc, features=RegulonName_AUC[i], label=T, reduction = 'umap')
  p2 = FeaturePlot(obj.combined.bin, features=RegulonName_BIN[i], label=T, reduction = 'umap')
  p3 = DimPlot(obj.combined, reduction = 'umap', group.by = "celltype", label=T)
  plotc = p2|p3
  print(plotc)
  dev.off()
}



#RidgePlot&VlnPlot
dir.create("../3group_RidgePlot_VlnPlot_auc")
setwd("../3group_RidgePlot_VlnPlot_auc")
for (i in 1:length(RegulonName_AUC)) {
  pdf(file=paste0(RegulonName_AUC[i],".pdf"),  width=10 ,height=12)
  p1 = RidgePlot(obj.combined.auc, features = RegulonName_AUC[i], group.by="celltype") + 
    theme(legend.position='none')
  p2 = VlnPlot(obj.combined.auc, features = RegulonName_AUC[i], pt.size = 0, group.by="celltype") + 
    theme(legend.position='none')
  plotc = p1 + p2
  print(plotc)
  dev.off()
}


dir.create("../3group_RidgePlot_VlnPlot_bin")
setwd("../3group_RidgePlot_VlnPlot_bin")
for (i in 1:length(RegulonName_BIN)) {
  pdf(file=paste0(RegulonName_BIN[i],".pdf"),  width=10 ,height=12)
  p1 = RidgePlot(obj.combined.bin, features = RegulonName_BIN[i], group.by="celltype") + 
    theme(legend.position='none')
  p2 = VlnPlot(obj.combined.bin, features = RegulonName_BIN[i], pt.size = 0, group.by="celltype") + 
    theme(legend.position='none')
  plotc = p1 + p2
  print(plotc)
  dev.off()
}





################################################################################
##################################################三组合并一块出图
###FeaturePlot&idgePlot&VlnPlot
dir.create("../3group_plot_split")
setwd("../3group_plot_split")

for (i in 1:length(RegulonName_BIN)) {
  pdf(paste0(RegulonName_BIN[i],".pdf"),  width=18 ,height=12)
  p1 = FeaturePlot(obj.combined.auc, features=RegulonName_BIN[i], label=T, reduction = 'umap')
  p2 = FeaturePlot(obj.combined.bin, features=RegulonName_BIN[i], label=T, reduction = 'umap')
  p3 = DimPlot(obj.combined, reduction = 'umap', group.by = "celltype", label=T)
  p4 = FeaturePlot(obj.combined.bin, features=RegulonName_BIN[i], label=T, reduction = 'umap',split.by = "group")
  plotc = (p1|p2|p3)/
    p4
  #plotc= p1 + p2 + p3 + p4 + plot_layout(nrow = 2)
  print(plotc)
  dev.off()
}



dir.create("../3group_RidgePlot_VlnPlot_bin")
setwd("../3group_RidgePlot_VlnPlot_bin")
for (i in 1:length(RegulonName_BIN)) {
  pdf(file=paste0(RegulonName_BIN[i],".pdf"),  width=10 ,height=12)
  p1 = RidgePlot(obj.combined.bin, features = RegulonName_BIN[i], group.by="celltype") + 
    theme(legend.position='none')
  p2 = VlnPlot(obj.combined.bin, features = RegulonName_BIN[i], pt.size = 0, group.by="celltype") + 
    theme(legend.position='none')
  p3 = RidgePlot(obj.combined.bin, features = RegulonName_BIN[i],group.by = "group",) + 
    theme(legend.position='none')
  p4 = VlnPlot(obj.combined.bin, features = RegulonName_BIN[i], pt.size = 0, split.by = "group") + 
    theme(legend.position='none')
  plotc = (p1 + p2)/(p3+p4)
  print(plotc)
  dev.off()
}

###################################################出版图
dir.create("../Published")
setwd("../Published")
for (i in 1:length(RegulonName_BIN)) {
  pdf(file=paste0(RegulonName_BIN[i],".pdf"),  width=8 ,height=6)
  p1=FeaturePlot(obj.combined.bin, features=RegulonName_BIN[i], label=T,label.size = 10, 
                 cols = c("lightgrey","red2"),#shape.by = "celltype",
                 pt.size = 3,repel=T,
                 reduction = 'umap')+theme_bw()+
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size =20),
          legend.text = element_text(size = 30),
          legend.title = element_text(size = 30),
          title = element_text(size = 40))+
    guides(color = guide_colorsteps(override.aes = list(size = 12)))
  print(p1)
  dev.off()
}
##############################################################
#利用pheatmap可视化
dir.create("../3group_pheatmap")
setwd("../3group_pheatmap")
celltype = subset(cellInfo,select = 'celltype')
AUCmatrix <- t(AUCmatrix)
BINmatrix <- t(BINmatrix)
#my.regulons <- c('Gata3_57g','Ets1_126g','Crem_24g','Nfkb1_64g','Eomes_extended_33g','Tbx21_14g','Irf5_27g') #挑选部分感兴趣的regulons
my.regulons <- RegulonName_BIN
myAUCmatrix <- AUCmatrix[rownames(AUCmatrix)%in%my.regulons,]
myBINmatrix <- BINmatrix[rownames(BINmatrix)%in%my.regulons,]
#使用regulon原始AUC值绘制热图
pheatmap(myAUCmatrix, show_colnames=F, annotation_col=celltype,
         cluster_cols = T,
         filename = 'myAUCmatrix_heatmap.pdf',
         width = 6, height = 10)
#使用regulon二进制AUC值绘制热图
pheatmap(myBINmatrix, show_colnames=F, annotation_col=celltype,
         filename = 'myBINmatrix_heatmap.pdf',
         color = colorRampPalette(colors =c("#330066","#FFCC33"))(100),   ###c("blue","white","red"))(100),
         # cutree_rows=5,
         #cutree_cols=10,
         cluster_rows = T,
         cluster_cols = T,
         fontsize = 20, 
         width =12, height = 12)


###################################################################################

###Cell-type specific regulators (RSS)
##*Cell-type specific regulators (based on the  Regulon Specificity Score (RSS) proposed by  Suo et al.* for the Mouse Cell Atlas in 2018).Useful for big analysis with many cell types, to identify the cell-type specific regulons.
##当细胞种类比较多时可以用RSS来识别细胞类型特异性regulons
##得到热点图，颜色深浅代表z-score值，点的大小代表RSS评分。
library(AUCell)
scenicOptions=readRDS(file="int/scenicOptions.Rds")
regulonAUC <-loadInt(scenicOptions, "aucell_regulonAUC")
cellInfo1 <- cellInfo
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo1[colnames(regulonAUC), "celltype"])
rssPlot <-plotRSS(rss) #大小 rss评分，颜色 Z-score
rssPlot[["plot"]][["data"]][["cellType"]] <- factor(rssPlot[["plot"]][["data"]][["cellType"]],levels = c("ILC2a","ILC2b","ILCdc","ILC1","ILCt"))
rssPlot$plot
ggsave("rssPlot.pdf",rssPlot$plot,width = 8,height =6)


