#两组细胞数目和比例对比
#整体细胞数
setwd("/home/data/t040243/CopyOfjiangshijiu")
dir.create("number and DEG")
setwd("number and DEG")
load("/home/data/t040243/CopyOfjiangshijiu/Doublet/after_doublet_obj.combined.RData")
cell.num <- table(Idents(obj.combined))
cell.freq <- round(prop.table(table(Idents(obj.combined)))*100,2)
cell.combined <- rbind(cell.num, cell.freq)
write.csv(file="combined_cell_counts_freq.csv",cell.combined)
#分组计算细胞数目和比例
cell.num.group <- table(Idents(obj.combined), obj.combined$group) 
colnames(cell.num.group)<-paste0(colnames(cell.num.group),'_cell_counts')
cell.freq.group <- round(prop.table(table(Idents(obj.combined), obj.combined$group), margin = 2) *100,2)
colnames(cell.freq.group)<-paste0(colnames(cell.freq.group),'_cell_Freq')
cell.group <- cbind(cell.num.group, cell.freq.group)
write.csv(file="group_cell_counts_freq.csv",cell.group)
#表格和UMAP展示细胞数目变化
pdf('UMAP_multi_samples_split_anno.pdf',width = 25,height = 8)
p<-DimPlot(obj.combined, reduction = "umap", split.by = "group",label=T,repel=T) + NoLegend()
tb <- tableGrob(cell.group)
plot_grid(p, tb,ncol=2,rel_widths=c(0.6,0.4))
dev.off()    
#火山图展示 文献：Integrated single cell analysis of blood and cerebrospinal fluid leukocytes in multiple sclerosis
#两组间差异基因
#将每个细胞的identity转换成celltype.group
DefaultAssay(obj.combined) <- "RNA"
obj.combined$celltype <- Idents(obj.combined) ##保存每个细胞的细胞类型
obj.combined$celltype.group <- paste(Idents(obj.combined), obj.combined$group, sep = "_")   ##细胞类型前面加上分组信息

######################################################################################
#######只需Ctrl+F替换ILC类型即可
########################################ILC2a细胞
##鉴定ILC2a和ILCreg细胞在MI和Health组中差异表达的基因
Idents(obj.combined) <- "celltype.group"
ILC2a.diff.MIvsSham <- FindMarkers(obj.combined, ident.1 = "ILC2a_MI", ident.2 = "ILC2a_Sham", verbose = FALSE)
write.csv(file="ILC2a_Diff_MIvsSham.csv",ILC2a.diff.MIvsSham)
ILC2a.diff.MIRIvsSham <- FindMarkers(obj.combined, ident.1 = "ILC2a_MIRI", ident.2 = "ILC2a_Sham", verbose = FALSE)
write.csv(file="ILC2a_Diff_MIRIvsSham.csv",ILC2a.diff.MIRIvsSham)

##提取出ILC2a细胞
Idents(obj.combined) <- 'celltype'
ILC2a <- subset(obj.combined, idents = "ILC2a")
Idents(ILC2a) <- "group"
avg.ILC2a <- log1p(AverageExpression(ILC2a, verbose = FALSE)$RNA)
avg.ILC2a <- data.frame(avg.ILC2a ,gene=rownames(avg.ILC2a))       ############可以用于Wgcna分析
##expr <- AverageExpression(obj.combined, assays = "RNA", slot = "data")[[1]]   ############可以用于Wgcna分析
##火山图可视化差异基因
low<-floor(range(ILC2a.diff.MIvsSham$avg_log2FC)[1])
high<-ceiling(range(ILC2a.diff.MIvsSham$avg_log2FC)[2])
pdf('DEG_ILC2a_MIvsSham.pdf',width = 12,height = 15)
print(EnhancedVolcano(ILC2a.diff.MIvsSham,
                      title = 'ILC2a MI versus Sham',
                      lab = rownames(ILC2a.diff.MIvsSham),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      FCcutoff = 1,
                      pCutoff = 1e-05,
                      ##设置标题、点和标签、轴的大小
                      titleLabSize = 45,
                      pointSize = 8.0,
                      labSize = 12,
                      axisLabSize = 50,
                      col=c('grey', 'blue', 'green','red3'),  ##颜色分别代表NS，log2FC，pvalue
                      colAlpha = 1, ##透明度（0-1之间）
                      legendPosition = 'top',  ##（'bottom',"right","left"
                      legendLabSize = 25,   ###图例标签大小
                      legendIconSize = 10,   ###图例点的大小 
                      #添加Connectors
                      drawConnectors = TRUE,  ##添加连线
                      widthConnectors = 0.5, ##连线的宽度
                      colConnectors = "red",  ##连线的颜色
                      max.overlaps = 10, ##连线达到最大重叠数目后就不显示标签了
                      xlim = c(low, high)))
dev.off()

low<-floor(range(ILC2a.diff.MIRIvsSham$avg_log2FC)[1])
high<-ceiling(range(ILC2a.diff.MIRIvsSham$avg_log2FC)[2])
pdf('DEG_ILC2a_MIRIvsSham.pdf',width = 12,height = 15)
print(EnhancedVolcano(ILC2a.diff.MIRIvsSham,
                      title = 'ILC2a MIRI versus Sham',
                      lab = rownames(ILC2a.diff.MIRIvsSham),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      FCcutoff = 1,
                      pCutoff = 1e-05,
                      ##设置标题、点和标签、轴的大小
                      titleLabSize = 45,
                      pointSize = 8.0,
                      labSize = 12,
                      axisLabSize = 50,
                      col=c('grey', 'blue', 'green','red3'),  ##颜色分别代表NS，log2FC，pvalue
                      colAlpha = 1, ##透明度（0-1之间）
                      legendPosition = 'top',  ##（'bottom',"right","left"
                      legendLabSize = 25,   ###图例标签大小
                      legendIconSize = 10,   ###图例点的大小 
                      #添加Connectors
                      drawConnectors = TRUE,  ##添加连线
                      widthConnectors = 0.5, ##连线的宽度
                      colConnectors = "red",  ##连线的颜色
                      max.overlaps = 10, ##连线达到最大重叠数目后就不显示标签了
                      xlim = c(low, high)))
dev.off()
##################################################################################
###########################################################################






###########################ILC2b
Idents(obj.combined) <- "celltype.group"
##鉴定ILC2b和ILCreg细胞在MI和Health组中差异表达的基因
ILC2b.diff.MIvsSham <- FindMarkers(obj.combined, ident.1 = "ILC2b_MI", ident.2 = "ILC2b_Sham", verbose = FALSE)
write.csv(file="ILC2b_Diff_MIvsSham.csv",ILC2b.diff.MIvsSham)
ILC2b.diff.MIRIvsSham <- FindMarkers(obj.combined, ident.1 = "ILC2b_MIRI", ident.2 = "ILC2b_Sham", verbose = FALSE)
write.csv(file="ILC2b_Diff_MIRIvsSham.csv",ILC2b.diff.MIRIvsSham)

##提取出ILC2b细胞
Idents(obj.combined) <- 'celltype'
ILC2b <- subset(obj.combined, idents = "ILC2b")
Idents(ILC2b) <- "group"
avg.ILC2b <- log1p(AverageExpression(ILC2b, verbose = FALSE)$RNA)
avg.ILC2b <- data.frame(avg.ILC2b ,gene=rownames(avg.ILC2b))       ############可以用于Wgcna分析

##火山图可视化差异基因
low<-floor(range(ILC2b.diff.MIvsSham$avg_log2FC)[1])
high<-ceiling(range(ILC2b.diff.MIvsSham$avg_log2FC)[2])
pdf('DEG_ILC2b_MIvsSham.pdf',width = 12,height = 15)
print(EnhancedVolcano(ILC2b.diff.MIvsSham,
                      title = 'ILC2b MI versus Sham',
                      lab = rownames(ILC2b.diff.MIvsSham),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      FCcutoff = 1,
                      pCutoff = 1e-05,
                      ##设置标题、点和标签、轴的大小
                      titleLabSize = 45,
                      pointSize = 8.0,
                      labSize = 12,
                      axisLabSize = 50,
                      col=c('grey', 'blue', 'green','red3'),  ##颜色分别代表NS，log2FC，pvalue
                      colAlpha = 1, ##透明度（0-1之间）
                      legendPosition = 'top',  ##（'bottom',"right","left"
                      legendLabSize = 25,   ###图例标签大小
                      legendIconSize = 10,   ###图例点的大小 
                      #添加Connectors
                      drawConnectors = TRUE,  ##添加连线
                      widthConnectors = 0.5, ##连线的宽度
                      colConnectors = "red",  ##连线的颜色
                      max.overlaps = 10, ##连线达到最大重叠数目后就不显示标签了
                      xlim = c(low, high)))
dev.off()

##############################################
low<-floor(range(ILC2b.diff.MIRIvsSham$avg_log2FC)[1])
high<-ceiling(range(ILC2b.diff.MIRIvsSham$avg_log2FC)[2])
pdf('DEG_ILC2b_MIRIvsSham.pdf',width = 12,height = 15)
print(EnhancedVolcano(ILC2b.diff.MIRIvsSham,
                      title = 'ILC2b MIRI versus Sham',
                      lab = rownames(ILC2b.diff.MIRIvsSham),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      FCcutoff = 1,
                      pCutoff = 1e-05,
                      ##设置标题、点和标签、轴的大小
                      titleLabSize = 45,
                      pointSize = 8.0,
                      labSize = 12,
                      axisLabSize = 50,
                      col=c('grey', 'blue', 'green','red3'),  ##颜色分别代表NS，log2FC，pvalue
                      colAlpha = 1, ##透明度（0-1之间）
                      legendPosition = 'top',  ##（'bottom',"right","left"
                      legendLabSize = 25,   ###图例标签大小
                      legendIconSize = 10,   ###图例点的大小 
                      #添加Connectors
                      drawConnectors = TRUE,  ##添加连线
                      widthConnectors = 0.5, ##连线的宽度
                      colConnectors = "red",  ##连线的颜色
                      max.overlaps = 10, ##连线达到最大重叠数目后就不显示标签了
                      xlim = c(low, high)))
dev.off()




############################################################
######ILC1
Idents(obj.combined) <- "celltype.group"
##鉴定ILC1和ILCreg细胞在MI和Health组中差异表达的基因
ILC1.diff.MIvsSham <- FindMarkers(obj.combined, ident.1 = "ILC1_MI", ident.2 = "ILC1_Sham", verbose = FALSE)
write.csv(file="ILC1_Diff_MIvsSham.csv",ILC1.diff.MIvsSham)
ILC1.diff.MIRIvsSham <- FindMarkers(obj.combined, ident.1 = "ILC1_MIRI", ident.2 = "ILC1_Sham", verbose = FALSE)
write.csv(file="ILC1_Diff_MIRIvsSham.csv",ILC1.diff.MIRIvsSham)

##提取出ILC1细胞
Idents(obj.combined) <- 'celltype'
ILC1 <- subset(obj.combined, idents = "ILC1")
Idents(ILC1) <- "group"
avg.ILC1 <- log1p(AverageExpression(ILC1, verbose = FALSE)$RNA)
avg.ILC1 <- data.frame(avg.ILC1 ,gene=rownames(avg.ILC1))       ############可以用于Wgcna分析

##火山图可视化差异基因
low<-floor(range(ILC1.diff.MIvsSham$avg_log2FC)[1])
high<-ceiling(range(ILC1.diff.MIvsSham$avg_log2FC)[2])
pdf('DEG_ILC1_MIvsSham.pdf',width = 12,height = 15)
print(EnhancedVolcano(ILC1.diff.MIvsSham,
                      title = 'ILC1 MI versus Sham',
                      lab = rownames(ILC1.diff.MIvsSham),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      FCcutoff = 1,
                      pCutoff = 1e-05,
                      ##设置标题、点和标签、轴的大小
                      titleLabSize = 45,
                      pointSize = 8.0,
                      labSize = 12,
                      axisLabSize = 50,
                      col=c('grey', 'blue', 'green','red3'),  ##颜色分别代表NS，log2FC，pvalue
                      colAlpha = 1, ##透明度（0-1之间）
                      legendPosition = 'top',  ##（'bottom',"right","left"
                      legendLabSize = 25,   ###图例标签大小
                      legendIconSize = 10,   ###图例点的大小 
                      #添加Connectors
                      drawConnectors = TRUE,  ##添加连线
                      widthConnectors = 0.5, ##连线的宽度
                      colConnectors = "red",  ##连线的颜色
                      max.overlaps = 10, ##连线达到最大重叠数目后就不显示标签了
                      xlim = c(low, high)))
dev.off()

##############################################
low<-floor(range(ILC1.diff.MIRIvsSham$avg_log2FC)[1])
high<-ceiling(range(ILC1.diff.MIRIvsSham$avg_log2FC)[2])
pdf('DEG_ILC1_MIRIvsSham.pdf',width = 12,height = 15)
print(EnhancedVolcano(ILC1.diff.MIRIvsSham,
                      title = 'ILC1 MIRI versus Sham',
                      lab = rownames(ILC1.diff.MIRIvsSham),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      FCcutoff = 1,
                      pCutoff = 1e-05,
                      ##设置标题、点和标签、轴的大小
                      titleLabSize = 45,
                      pointSize = 8.0,
                      labSize = 12,
                      axisLabSize = 50,
                      col=c('grey', 'blue', 'green','red3'),  ##颜色分别代表NS，log2FC，pvalue
                      colAlpha = 1, ##透明度（0-1之间）
                      legendPosition = 'top',  ##（'bottom',"right","left"
                      legendLabSize = 25,   ###图例标签大小
                      legendIconSize = 10,   ###图例点的大小 
                      #添加Connectors
                      drawConnectors = TRUE,  ##添加连线
                      widthConnectors = 0.5, ##连线的宽度
                      colConnectors = "red",  ##连线的颜色
                      max.overlaps = 10, ##连线达到最大重叠数目后就不显示标签了
                      xlim = c(low, high)))
dev.off()







###########################################################
####DcILC
Idents(obj.combined) <- "celltype.group"
##鉴定ILCdc和ILCreg细胞在MI和Health组中差异表达的基因
ILCdc.diff.MIvsSham <- FindMarkers(obj.combined, ident.1 = "ILCdc_MI", ident.2 = "ILCdc_Sham", verbose = FALSE)
write.csv(file="ILCdc_Diff_MIvsSham.csv",ILCdc.diff.MIvsSham)
ILCdc.diff.MIRIvsSham <- FindMarkers(obj.combined, ident.1 = "ILCdc_MIRI", ident.2 = "ILCdc_Sham", verbose = FALSE)
write.csv(file="ILCdc_Diff_MIRIvsSham.csv",ILCdc.diff.MIRIvsSham)

##提取出ILCdc细胞
Idents(obj.combined) <- 'celltype'
ILCdc <- subset(obj.combined, idents = "ILCdc")
Idents(ILCdc) <- "group"
avg.ILCdc <- log1p(AverageExpression(ILCdc, verbose = FALSE)$RNA)
avg.ILCdc <- data.frame(avg.ILCdc ,gene=rownames(avg.ILCdc))       ############可以用于Wgcna分析

##火山图可视化差异基因
low<-floor(range(ILCdc.diff.MIvsSham$avg_log2FC)[1])
high<-ceiling(range(ILCdc.diff.MIvsSham$avg_log2FC)[2])
pdf('DEG_ILCdc_MIvsSham.pdf',width = 12,height = 15)
print(EnhancedVolcano(ILCdc.diff.MIvsSham,
                      title = 'ILCdc MI versus Sham',
                      lab = rownames(ILCdc.diff.MIvsSham),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      FCcutoff = 1,
                      pCutoff = 1e-05,
                      ##设置标题、点和标签、轴的大小
                      titleLabSize = 45,
                      pointSize = 8.0,
                      labSize = 12,
                      axisLabSize = 50,
                      col=c('grey', 'blue', 'green','red3'),  ##颜色分别代表NS，log2FC，pvalue
                      colAlpha = 1, ##透明度（0-1之间）
                      legendPosition = 'top',  ##（'bottom',"right","left"
                      legendLabSize = 25,   ###图例标签大小
                      legendIconSize = 10,   ###图例点的大小 
                      #添加Connectors
                      drawConnectors = TRUE,  ##添加连线
                      widthConnectors = 0.5, ##连线的宽度
                      colConnectors = "red",  ##连线的颜色
                      max.overlaps = 10, ##连线达到最大重叠数目后就不显示标签了
                      xlim = c(low, high)))
dev.off()

##############################################
low<-floor(range(ILCdc.diff.MIRIvsSham$avg_log2FC)[1])
high<-ceiling(range(ILCdc.diff.MIRIvsSham$avg_log2FC)[2])
pdf('DEG_ILCdc_MIRIvsSham.pdf',width = 12,height = 15)
print(EnhancedVolcano(ILCdc.diff.MIRIvsSham,
                      title = 'ILCdc MIRI versus Sham',
                      lab = rownames(ILCdc.diff.MIRIvsSham),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      FCcutoff = 1,
                      pCutoff = 1e-05,
                      ##设置标题、点和标签、轴的大小
                      titleLabSize = 45,
                      pointSize = 8.0,
                      labSize = 12,
                      axisLabSize = 50,
                      col=c('grey', 'blue', 'green','red3'),  ##颜色分别代表NS，log2FC，pvalue
                      colAlpha = 1, ##透明度（0-1之间）
                      legendPosition = 'top',  ##（'bottom',"right","left"
                      legendLabSize = 25,   ###图例标签大小
                      legendIconSize = 10,   ###图例点的大小 
                      #添加Connectors
                      drawConnectors = TRUE,  ##添加连线
                      widthConnectors = 0.5, ##连线的宽度
                      colConnectors = "red",  ##连线的颜色
                      max.overlaps = 10, ##连线达到最大重叠数目后就不显示标签了
                      xlim = c(low, high)))
dev.off()






################ILCt
Idents(obj.combined) <- "celltype.group"
##鉴定ILCt和ILCreg细胞在MI和Health组中差异表达的基因
ILCt.diff.MIvsSham <- FindMarkers(obj.combined, ident.1 = "ILCt_MI", ident.2 = "ILCt_Sham", verbose = FALSE)
write.csv(file="ILCt_Diff_MIvsSham.csv",ILCt.diff.MIvsSham)
ILCt.diff.MIRIvsSham <- FindMarkers(obj.combined, ident.1 = "ILCt_MIRI", ident.2 = "ILCt_Sham", verbose = FALSE)
write.csv(file="ILCt_Diff_MIRIvsSham.csv",ILCt.diff.MIRIvsSham)

##提取出ILCt细胞
Idents(obj.combined) <- 'celltype'
ILCt <- subset(obj.combined, idents = "ILCt")
Idents(ILCt) <- "group"
avg.ILCt <- log1p(AverageExpression(ILCt, verbose = FALSE)$RNA)
avg.ILCt <- data.frame(avg.ILCt ,gene=rownames(avg.ILCt))       ############可以用于Wgcna分析

##火山图可视化差异基因
low<-floor(range(ILCt.diff.MIvsSham$avg_log2FC)[1])
high<-ceiling(range(ILCt.diff.MIvsSham$avg_log2FC)[2])
pdf('DEG_ILCt_MIvsSham.pdf',width = 12,height = 15)
print(EnhancedVolcano(ILCt.diff.MIvsSham,
                      title = 'ILCt MI versus Sham',
                      lab = rownames(ILCt.diff.MIvsSham),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      FCcutoff = 1,
                      pCutoff = 1e-05,
                      ##设置标题、点和标签、轴的大小
                      titleLabSize = 45,
                      pointSize = 8.0,
                      labSize = 12,
                      axisLabSize = 50,
                      col=c('grey', 'blue', 'green','red3'),  ##颜色分别代表NS，log2FC，pvalue
                      colAlpha = 1, ##透明度（0-1之间）
                      legendPosition = 'top',  ##（'bottom',"right","left"
                      legendLabSize = 25,   ###图例标签大小
                      legendIconSize = 10,   ###图例点的大小 
                      #添加Connectors
                      drawConnectors = TRUE,  ##添加连线
                      widthConnectors = 0.5, ##连线的宽度
                      colConnectors = "red",  ##连线的颜色
                      max.overlaps = 10, ##连线达到最大重叠数目后就不显示标签了
                      xlim = c(low, high)))
dev.off()

##############################################
low<-floor(range(ILCt.diff.MIRIvsSham$avg_log2FC)[1])
high<-ceiling(range(ILCt.diff.MIRIvsSham$avg_log2FC)[2])
pdf('DEG_ILCt_MIRIvsSham.pdf',width = 12,height = 15)
print(EnhancedVolcano(ILCt.diff.MIRIvsSham,
                      title = 'ILCt MIRI versus Sham',
                      lab = rownames(ILCt.diff.MIRIvsSham),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      FCcutoff = 1,
                      pCutoff = 1e-05,
                      ##设置标题、点和标签、轴的大小
                      titleLabSize = 45,
                      pointSize = 8.0,
                      labSize = 12,
                      axisLabSize = 50,
                      col=c('grey', 'blue', 'green','red3'),  ##颜色分别代表NS，log2FC，pvalue
                      colAlpha = 1, ##透明度（0-1之间）
                      legendPosition = 'top',  ##（'bottom',"right","left"
                      legendLabSize = 25,   ###图例标签大小
                      legendIconSize = 10,   ###图例点的大小 
                      #添加Connectors
                      drawConnectors = TRUE,  ##添加连线
                      widthConnectors = 0.5, ##连线的宽度
                      colConnectors = "red",  ##连线的颜色
                      max.overlaps = 10, ##连线达到最大重叠数目后就不显示标签了
                      xlim = c(low, high)))
dev.off()

###############################################################################
#################################################################################
###################################################################################
#############################################ILC亚群之间整体对比
####DcILC Vs ILC2a
Idents(obj.combined) <- "celltype"
ILCdc_VS_ILC2a.diff <- FindMarkers(obj.combined, ident.1 = "ILCdc", ident.2 = "ILC2a", verbose = FALSE)
write.csv(file="ILCdc_VS_ILC2a.csv",ILCdc_VS_ILC2a.diff)


##火山图可视化差异基因
low<-floor(range(ILCdc_VS_ILC2a.diff$avg_log2FC)[1])
high<-ceiling(range(ILCdc_VS_ILC2a.diff$avg_log2FC)[2])
pdf('ILCdc_VS_ILC2a.pdf',width = 12,height = 15)
print(EnhancedVolcano(ILCdc_VS_ILC2a.diff,
                      title = 'ILCdc_VS_ILC2a',
                      lab = rownames(ILCdc_VS_ILC2a.diff),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      FCcutoff = 1,
                      pCutoff = 1e-05,
                      ##设置标题、点和标签、轴的大小
                      titleLabSize = 45,
                      pointSize = 8.0,
                      labSize = 12,
                      axisLabSize = 50,
                      col=c('grey', 'blue', 'green','red3'),  ##颜色分别代表NS，log2FC，pvalue
                      colAlpha = 1, ##透明度（0-1之间）
                      legendPosition = 'top',  ##（'bottom',"right","left"
                      legendLabSize = 25,   ###图例标签大小
                      legendIconSize = 10,   ###图例点的大小 
                      #添加Connectors
                      drawConnectors = TRUE,  ##添加连线
                      widthConnectors = 0.5, ##连线的宽度
                      colConnectors = "red",  ##连线的颜色
                      max.overlaps = 10, ##连线达到最大重叠数目后就不显示标签了
                      xlim = c(low, high)))
dev.off()
########################################################################
########################################################################
########################################################################
####ILCdc Vs ILC2b
Idents(obj.combined) <- "celltype"
ILCdc_VS_ILC2b.diff <- FindMarkers(obj.combined, ident.1 = "ILCdc", ident.2 = "ILC2b", verbose = FALSE)
write.csv(file="ILCdc_VS_ILC2b.csv",ILCdc_VS_ILC2b.diff)


##火山图可视化差异基因
low<-floor(range(ILCdc_VS_ILC2b.diff$avg_log2FC)[1])
high<-ceiling(range(ILCdc_VS_ILC2b.diff$avg_log2FC)[2])
pdf('ILCdc_VS_ILC2b.pdf',width = 12,height = 15)
print(EnhancedVolcano(ILCdc_VS_ILC2b.diff,
                      title = 'ILCdc_VS_ILC2b',
                      lab = rownames(ILCdc_VS_ILC2b.diff),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      FCcutoff = 1,
                      pCutoff = 1e-05,
                      ##设置标题、点和标签、轴的大小
                      titleLabSize = 45,
                      pointSize = 8.0,
                      labSize = 12,
                      axisLabSize = 50,
                      col=c('grey', 'blue', 'green','red3'),  ##颜色分别代表NS，log2FC，pvalue
                      colAlpha = 1, ##透明度（0-1之间）
                      legendPosition = 'top',  ##（'bottom',"right","left"
                      legendLabSize = 25,   ###图例标签大小
                      legendIconSize = 10,   ###图例点的大小 
                      #添加Connectors
                      drawConnectors = TRUE,  ##添加连线
                      widthConnectors = 0.5, ##连线的宽度
                      colConnectors = "red",  ##连线的颜色
                      max.overlaps = 10, ##连线达到最大重叠数目后就不显示标签了
                      xlim = c(low, high)))
dev.off()




###############################################################################
#################################################################################
###################################################################################
#############################################ILC亚群之间整体对比
####ILCdc Vs ILC1
Idents(obj.combined) <- "celltype"
ILCdc_VS_ILC1.diff <- FindMarkers(obj.combined, ident.1 = "ILCdc", ident.2 = "ILC1", verbose = FALSE)
write.csv(file="ILCdc_VS_ILC1.csv",ILCdc_VS_ILC1.diff)


##火山图可视化差异基因
low<-floor(range(ILCdc_VS_ILC1.diff$avg_log2FC)[1])
high<-ceiling(range(ILCdc_VS_ILC1.diff$avg_log2FC)[2])
pdf('ILCdc_VS_ILC1.pdf',width = 12,height = 15)
print(EnhancedVolcano(ILCdc_VS_ILC1.diff,
                      title = 'ILCdc_VS_ILC1',
                      lab = rownames(ILCdc_VS_ILC1.diff),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      FCcutoff = 1,
                      pCutoff = 1e-05,
                      ##设置标题、点和标签、轴的大小
                      titleLabSize = 45,
                      pointSize = 8.0,
                      labSize = 12,
                      axisLabSize = 50,
                      col=c('grey', 'blue', 'green','red3'),  ##颜色分别代表NS，log2FC，pvalue
                      colAlpha = 1, ##透明度（0-1之间）
                      legendPosition = 'top',  ##（'bottom',"right","left"
                      legendLabSize = 25,   ###图例标签大小
                      legendIconSize = 10,   ###图例点的大小 
                      #添加Connectors
                      drawConnectors = TRUE,  ##添加连线
                      widthConnectors = 0.5, ##连线的宽度
                      colConnectors = "red",  ##连线的颜色
                      max.overlaps = 10, ##连线达到最大重叠数目后就不显示标签了
                      xlim = c(low, high)))
dev.off()





###############################################################################
#################################################################################
###################################################################################
#############################################ILC亚群之间整体对比
####ILCdc Vs ILCt
Idents(obj.combined) <- "celltype"
ILCdc_VS_ILCt.diff <- FindMarkers(obj.combined, ident.1 = "ILCdc", ident.2 = "ILCt", verbose = FALSE)
write.csv(file="ILCdc_VS_ILCt.csv",ILCdc_VS_ILCt.diff)


##火山图可视化差异基因
low<-floor(range(ILCdc_VS_ILCt.diff$avg_log2FC)[1])
high<-ceiling(range(ILCdc_VS_ILCt.diff$avg_log2FC)[2])
pdf('ILCdc_VS_ILCt.pdf',width = 12,height = 15)
print(EnhancedVolcano(ILCdc_VS_ILCt.diff,
                      title = 'ILCdc_VS_ILCt',
                      lab = rownames(ILCdc_VS_ILCt.diff),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      FCcutoff = 1,
                      pCutoff = 1e-05,
                      ##设置标题、点和标签、轴的大小
                      titleLabSize = 45,
                      pointSize = 8.0,
                      labSize = 12,
                      axisLabSize = 50,
                      col=c('grey', 'blue', 'green','red3'),  ##颜色分别代表NS，log2FC，pvalue
                      colAlpha = 1, ##透明度（0-1之间）
                      legendPosition = 'top',  ##（'bottom',"right","left"
                      legendLabSize = 25,   ###图例标签大小
                      legendIconSize = 10,   ###图例点的大小 
                      #添加Connectors
                      drawConnectors = TRUE,  ##添加连线
                      widthConnectors = 0.5, ##连线的宽度
                      colConnectors = "red",  ##连线的颜色
                      max.overlaps = 10, ##连线达到最大重叠数目后就不显示标签了
                      xlim = c(low, high)))
dev.off()