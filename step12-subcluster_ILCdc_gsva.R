#Idents(obj.combined) <- 'celltype.group'
obj.combined <- sce
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


mycol = colorRampPalette(c("blue", "white","red"))(100)
#创建分组颜色条所需的数据框；
celltype = factor(c(rep("ILCdc_cluster0",1), rep("ILCdc_cluster1",1), rep("ILCdc_cluster2",1),rep("ILCdc_cluster3",1),rep("ILCdc_cluster4",1)))
anndf <- data.frame(celltype)
rownames(anndf) <- c("ILCdc_cluster0" ,"ILCdc_cluster1" ,"ILCdc_cluster2","ILCdc_cluster3"  ,"ILCdc_cluster4" )
#自定义分组颜色条的颜色；
anncol = list(celltype=c(ILCdc_cluster0="#c77cff",ILCdc_cluster1="#FF9999",ILCdc_cluster2="#99CC00",ILCdc_cluster3="#FF9900",ILCdc_cluster4="#FFCD71"))

#绘制热图；
set.seed(123)  #设置随机数种子，使结果可重复
tmp <- sample(rownames(gsva.res),100) %>% sort()
#scRNA <- gsva.res[tmp,]
scRNA <- gsva.res[c(1,	3,	4,	5,	10,	11,	17,	18,	24,	25,	26,	33,	34,	35,	36,	43,	48,	49,	56,	57,	58,	59,	60,	61,	69,	70,	87,	88,	89,	90,	91,	92,	95,	96,	97,	98,	99,	102,	103,	109,	110,	111,	114,	115,	116,	118,	120,	125,	128,	130,	131,	133,	137,	138,	139,	142,	148,	157,	158,	186),]
#对热图中的文字字号，聚类树高度，渐变颜色条做了自定义；
#同时，根据聚类结果在行和列方向都添加了“gap”；
pdf(file="GSVA_heatmap.pdf",width=14,height=12)
htmap <- pheatmap(scRNA, scale = "row",
                  border="white",
                  cellwidth =60,
                  cellheight = 10,
                  color = mycol,
                  cutree_rows=5,
                  cutree_cols=5,
                  fontsize = 15,
                  fontsize_col = 20,
                  fontsize_row = 10,
                  treeheight_row=30,
                  treeheight_col=30,
                  annotation_col=anndf,
                  annotation_colors=anncol,
                  annotation_legend=T,
                  legend = T,
                  show_colnames = T,
                  show_row_dend = F,
                  row_names_side = "left", 
                  #legend_breaks=c(-1,-0.5,0,0.5,1),
                  # legend_labels=c(-1,-0.5,0,0.5,1),
                  #main = "GSVA"
                  #angle_col = 90
)
htmap
dev.off()
##################