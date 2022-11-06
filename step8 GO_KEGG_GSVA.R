#富集分析GO和KEGG
#Mac细胞两组的差异基因GO富集
setwd("/home/shpc_100537/CopyOfjiangshijiu")
dir.create("GO_KEGG")
setwd("/home/shpc_100537/CopyOfjiangshijiu/GO_KEGG")
###有意义的差异基因太少了，就使用全部的差异基因吧！！！
#sig.ILC2a_VS_DcILC.diff.MIvsSham <- subset(ILC2a_VS_DcILC.diff.MIvsSham, p_val_adj<0.01&abs(avg_log2FC)>1)
#sig.ILC2a_VS_DcILC.diff.MIRIvsSham <- subset(ILC2a_VS_DcILC.diff.MIRIvsSham, p_val_adj<0.01&abs(avg_log2FC)>1)
#sig.Myeloid_DC.diff.MIRIvsSham <- subset(Myeloid_DC.diff.MIRIvsSham, p_val_adj<0.01&abs(avg_log2FC)>1)


ego_ALL <- enrichGO(gene        = row.names(ILCt.diff.MIvsSham),
                    OrgDb         = 'org.Mm.eg.db',
                    keyType       = 'SYMBOL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
ego_all <- data.frame(ego_ALL)
write.csv(file="ILCt_MI_vs_Sham_GO_result.csv",data.frame(ego_all))
#GO中CC,MF,BP
ego_CC <- enrichGO(gene          = row.names(ILCt.diff.MIvsSham),
                   OrgDb         = 'org.Mm.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_MF <- enrichGO(gene          = row.names(ILCt.diff.MIvsSham),
                   OrgDb         = 'org.Mm.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_BP <- enrichGO(gene          = row.names(ILCt.diff.MIvsSham),
                   OrgDb         = 'org.Mm.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
#可视化GO富集结果
ego_CC@result$Description <- substring(ego_CC@result$Description,1,70)
ego_MF@result$Description <- substring(ego_MF@result$Description,1,70)
ego_BP@result$Description <- substring(ego_BP@result$Description,1,70)   #截取结果中Description的1-70个字符，方便图形显示
pdf(file="ILCt_MI_vs_Sham_GO_barplot.pdf",width=20,height = 20)
p_BP <- barplot(ego_BP, showCategory = 5, label_format = 80) + ggtitle("barplot for Biological process")
p_CC <- barplot(ego_CC, showCategory = 5,label_format = 80) + ggtitle("barplot for Cellular component")
p_MF <- barplot(ego_MF, showCategory = 5,label_format = 80) + ggtitle("barplot for Molecular function")
plotc <- p_BP/p_CC/p_MF
plotc
dev.off()
pdf(file="ILCt_MI_vs_Sham_GO_dot.pdf",width=20,height = 20)
enrichplot::dotplot(ego_ALL,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
dev.off()

#KEGG富集
genelist <- bitr(row.names(ILCt.diff.MIvsSham), fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Mm.eg.db')
genelist <- pull(genelist, ENTREZID)                          
ekegg <- enrichKEGG(gene = genelist, organism = 'mmu')
write.csv(data.frame(ekegg),'ILCt_MI_vs_Sham_enrichKEGG.csv',row.names = F)
#KEGG可视化
pdf(file="ILCt_MI_vs_Sham_KEGG.pdf",width=20,height = 15)
p1 <- barplot(ekegg, showCategory=10,label_format=80)
p2 <- enrichplot::dotplot(ekegg, showCategory=10,label_format=80)
plotc = p1/p2
plotc
dev.off()


#GSEA和GSVA
#GSEA 挑选两组ILCt_MI_vs_Sham细胞的差异基因
geneList= ILCt.diff.MIvsSham$avg_log2FC ##制作基因list
names(geneList)= toupper(rownames(ILCt.diff.MIvsSham))  ##name为基因名字
geneList=sort(geneList,decreasing = T)  ###按FC从大到小排序
head(geneList)
#从GSEA官网下载GSEA分析需要的基因集 http://www.gsea-msigdb.org/gsea/index.jsp
#下载免疫相关的基因集，c7: immunologic signature gene sets
##gmtfile ='c7.immunesigdb.v7.5.1.symbols.gmt'  也可以选这个
gmtfile ='c2.cp.kegg.v7.5.1.symbols.gmt'
pathway<-read.gmt(gmtfile) ##读取gmt文件中的pathway信息
y <- GSEA(geneList,TERM2GENE =pathway) ##进行GSEA分析
write.csv(file="ILCt_MI_vs_Sham_GSEA_result.csv",data.frame(y))
#气泡图展示显著富集的前30条通路
pdf(file="ILCt_MI_vs_Sham_GSEA_dotplot.pdf",width=18)
clusterProfiler::dotplot(y,showCategory=10,label_format=80)
dev.off()
#绘制具体通路的GSEA图
pdf(file="ILCt_MI_vs_Sham_GSEA_.pdf",width=18)
gseaplot2(y,geneSetID = 'KEGG_RIBOSOME',pvalue_table=T)
dev.off()



#######################################################################################################
############################GSVA
#Idents(obj.combined) <- 'celltype.group'
Idents(obj.combined) <- 'celltype'
Idents(obj.combined)   #获取细胞类型
expr <- AverageExpression(obj.combined, assays = "RNA", slot = "data")[[1]]
View(expr) #计算每个基因在每个细胞亚群中的平均表达值
expr <- expr[rowSums(expr)>0,]  #选取非零基因
rownames(expr) <- toupper(rownames(expr))  #大写基因名
expr <- as.matrix(expr)  #转换成矩阵
#从GSEA官网下载GSEA分析需要的基因集
gmtfile ='c2.cp.kegg.v7.5.1.symbols.gmt'
pathway<-read.gmt(gmtfile)[,c(2,1)]  #读取gmt文件中的pathway信息
genesets=unstack(pathway)  #去堆叠,转换成list
gsva.res <- gsva(expr, genesets, method="ssgsea")   #进行GSVA分析 
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "gsva_res.csv", row.names = F)
##绘制热图可视化
#######################################################
######################################################
#自定义渐变色；
mycol = colorRampPalette(c("blue", "white","red"))(100)
#创建分组颜色条所需的数据框；
celltype = factor(c(rep("ILC2a",1), rep("ILC2b",1), rep("ILCdc",1),rep("ILC1",1),rep("ILCt",1)))
anndf <- data.frame(celltype)
rownames(anndf) <- c("ILC2a" ,"ILC2b" ,"ILCdc","ILC1"  ,"ILCt" )
#自定义分组颜色条的颜色；
anncol = list(celltype=c(ILC2a="#c77cff",ILC2b="#FF9999",ILCdc="#99CC00",ILC1="#FF9900",ILCt="#FFCD71"))

#绘制热图；
set.seed(123)  #设置随机数种子，使结果可重复
#tmp <- sample(rownames(gsva.res),100) %>% sort()  ###随机抽取100个
#scRNA <- gsva.res[tmp,]
scRNA <- gsva.res[c(1,	2,	3,	6,	8,	10,	11,	14,	15,	16,	17,	18,	21,	24,	25,	26,	27,	33,	34,	35,	36,	37,	38,	42,	43,	44,	47,	48,	49,	50,	51,	52,	53,	54,	57,	58,	59,	60,	61,	65,	66,	67,	68,	69,	70,	71,	72,	73,	85,	86,	87,	88,	89,	90,	91,	92,	93,	94,	95,	96,	97,	98,	99,	101,	102,	103,	106,	107,	109,	110,	111,	113,	114,	115,	116,	118,	120,	122,	124,	125,	126,	127,	128,	129,	130,	131,	132,	133,	135,	137,	138,	139,	140,	141,	142,	148,	153,	157,	158,	186),]
###挑选中意的通路

#对热图中的文字字号，聚类树高度，渐变颜色条做了自定义；
#同时，根据聚类结果在行和列方向都添加了“gap”；
pdf(file="GSVA_heatmap.pdf",width=14,height=18)
htmap <- pheatmap(scRNA, scale = "row",
                  border="white",
                  cellwidth =60,
                  cellheight = 10,
                  color = mycol,
                  cutree_rows=5,
                  cutree_cols=5,
                  fontsize = 15,
                  fontsize_col = 20,
                  fontsize_row = 12,
                  treeheight_row=30,
                  treeheight_col=30,
                  annotation_col=anndf,
                  annotation_colors=anncol,
                  annotation_legend=T,
                  legend = T,
                  show_colnames = T,
                  show_row_dend = F,
                  row_names_side = "left"
                  #legend_breaks=c(-1,-0.5,0,0.5,1),
                  # legend_labels=c(-1,-0.5,0,0.5,1),
                  #main = "GSVA"
                  #angle_col = 90
)
htmap
dev.off()
##########################################################
#######################################################
if(F){############绘制指定位置通路
ha = rowAnnotation(foo = anno_mark(at = c(1:50), 
                                   labels = rownames(gsva.res)[1:50]))
pdf(file="GSVA_heatmap.pdf",width=12,height=40)
pheatmap(gsva.res,show_colnames = T, scale = "row",show_rownames = F,
         display_numbers = F,fontsize_number = 20, cluster_cols = T, # 去掉横向、纵向聚类
         cluster_rows = T,  cutree_cols = 5, cutree_rows =5,right_annotation = ha,
        # row_names_side = "left", 
         cellheight=10,fontsize_row=6,fontsize_col = 16)
dev.off()
}

if (F) {
  #绘制热图可视化,展示前10个。
pdf(file="GSVA_heatmap_50.pdf",width=12,height=10)
pheatmap(gsva.res[c(1:50),], show_colnames = T, scale = "row",cellheight=12,fontsize_row=14)
dev.off()
}






#####################
if (F) {######绘制GSVA_FeaturePlot图
library(Seurat)
library(GSVA)
library(tidyverse)
##创建gmt文件转list函数
gmt2list <- function(gmtfile){
  sets <- as.list(read_lines(gmtfile))
  for(i in 1:length(sets)){
    tmp = str_split(sets[[i]], '\t')
    n = length(tmp[[1]])
    names(sets)[i] = tmp[[1]][1]
    sets[[i]] = tmp[[1]][3:n]
    rm(tmp, n)
  }
  return(sets)
}
#读取基因集数据库
s.sets = gmt2list("c2.cp.kegg.v7.5.1.symbols.gmt")
#读取表达矩阵
scRNA <- obj.combined
# 随机提取1000个细胞演示GSVA，非常规操作
# tmp <- sample(colnames(scRNA),1000) %>% sort()
# scRNA <- scRNA[,tmp]
expr <- as.matrix(scRNA@assays$RNA@counts)
meta <- scRNA@meta.data[,c("seurat_clusters", "celltype")]
es.matrix = gsva(expr, s.sets, kcdf="Poisson")
write.table(es.matrix, 'GSVA.xls', row.names=T, col.names=NA, sep='\t')

library(pheatmap)
library(patchwork)
#绘制热图
pheatmap(es.matrix, show_rownames=1, show_colnames=0, annotation_col=meta,
         fontsize_row=15, filename='GSVA.png', width=15, height=12)
#挑选感兴趣的基因集绘制featureplot
es <- data.frame(t(es.matrix),stringsAsFactors=F)
scRNA <- AddMetaData(scRNA, es)
p1 <- FeaturePlot(scRNA, features = "KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS", reduction = 'umap')
p2 <- FeaturePlot(scRNA, features = "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES", reduction = 'umap')
p3 <- FeaturePlot(scRNA, features = "KEGG_LEISHMANIA_INFECTION", reduction = 'umap')
plotc = (p1|p2)/(p3)
ggsave('GSVA_featureplot_demo.png', plotc, width = 10, height = 8)
}
