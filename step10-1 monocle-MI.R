###伪时序分析
#Monocle进行伪时间分析的核心技术是一种机器学习算法——反向图形嵌入 (Reversed Graph Embedding)。
#它分析的前提需要一张展现细胞转录特征相似性关系的图，
#Monocle2使用DDTree降维图，Monocle3使用UMAP降维图。
#Monocle的机器学习算法可以依据上述降维图形，学习描述细胞如何从一种状态过渡到另一种状态的轨迹。
#Monocle假设轨迹是树状结构，一端是“根”，另一端是“叶”。一个细胞在生物过程的开始，
#从根开始沿着主干进行，直到它到达第一个分支。然后，该细胞必须选择一条路径，并沿着树移动越来越远，
#直到它到达一片叶子。一个细胞的假时间值是它返回根所需的距离。
#降维方面monocle与seurat的过程大同小异，首先进行数据标准化，其次选择部分基因代表细胞转录特征，
#最后选用适当的算法降维。对Monocle原理感兴趣的同学可以登录官网查看：http://cole-trapnell-lab.github.io/monocle-release/


#数据导入与处理
#轨迹分析的前提是待分析的细胞有紧密的发育关系，PBMC细胞不是很好的的示例数据，
#我们选择T细胞群体演示一下。Monocle建议导入原始表达矩阵，由它完成数据标准化和其他预处理。

rm(list=ls())
set.seed(123)
setwd("/home/data/t040243/CopyOfjiangshijiu")
dir.create("monocle10-1")
setwd("/home/data/t040243/CopyOfjiangshijiu/monocle10-1")
load("/home/data/t040243/CopyOfjiangshijiu/Doublet/after_doublet_obj.combined.RData")
########################################################################
#挑选疾病组的所有细胞亚群进行拟时序分析(obj.combine为上一步中挑选MI组作为分析,obj.combined为全部组作为分析)
DefaultAssay(obj.combined) <- "RNA"   #更改默认数组,最好选择原始表达矩阵
scRNA<-obj.combined
scRNAsub = obj.combined[,obj.combined@meta.data$orig.ident %in% c("MI")]

dir.create("MI")
setwd("MI")
#挑选疾病组的所有细胞亚群进行拟时序分析(obj.combine为上一步中挑选MI组作为分析)

data <- as(as.matrix(scRNAsub@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNAsub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())

##expressionFamily参数用于指定表达矩阵的数据类型，有几个选项可以选择：
#稀疏矩阵用negbinomial.size()，
#FPKM值用tobit()，
#logFPKM值用gaussianff()
#mycds是Monocle为我们的数据生成的对象，相当于我们在seurat使用的scRNA对象


#数据导入后需要进行标准化和其他预处理：
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)
#mycds <- detectGenes(mycds, min_expr = 2)  #很多教程不用

#与seurat把标准化后的表达矩阵保存在对象中不同，monocle只保存一些中间结果在对象中，
#需要用时再用这些中间结果转化。经过上面三个函数的计算，
#mycds对象中多了SizeFactors、Dipersions、num_cells_expressed和num_genes_expressed等信息。

#提选择代表性基因
################################
#完成数据导入和预处理后，就可以考虑选择哪些基因代表细胞的发育特征
#Monocle官网教程提供了4个选择方法：
#（1）选择发育差异表达基因
#（2）选择clusters差异表达基因
#（3）选择离散程度高的基因
#（4）自定义发育marker基因
#前三种都是无监督分析方法，细胞发育轨迹生成完全不受人工干预；
#最后一种是半监督分析方法，可以使用先验知识辅助分析。
#第一种方法要求实验设计有不同的时间点，对起点和终点的样本做基因表达差异分析，
#挑选显著差异的基因进行后续分析。
#对于没有时序设计的实验样本，可以使用第2、3种方法挑选基因。
#第2种方法要先对细胞降维聚类，然后用clusters之间差异表达的基因开展后续分析。
#Monocle有一套自己的降维聚类方法，与seurat的方法大同小异，很多教程直接使用seurat的差异分析结果。
#第3种方法使用离散程度高的基因开展分析，seurat有挑选高变基因的方法，monocle也有自己选择的算法。
#本案例数据不具备使用第1、4种方法的条件，因此这里只演示2、3种方法的使用。
if (F) {
  
  ##使用clusters差异表达基因
  diff.wilcox = FindAllMarkers(scRNAsub)
  
  all.markers = diff.wilcox %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)  ##用dplyr::select替换select
  top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)   ##注意log角标有2
  #diff.genes <- read.csv('subcluster/diff_genes_wilcox.csv')
  diff.genes <- subset(diff.wilcox,p_val_adj<0.01)$gene
  diff.genes <- subset(diff.wilcox,p_val_adj<0.0001&abs(avg_log2FC)>0.75)$gene   ##标准更严格了
  mycds <- setOrderingFilter(mycds, diff.genes)
  p1 <- plot_ordering_genes(mycds);p1
  ##使用seurat选择的高变基因
  var.genes <- VariableFeatures(scRNAsub)
  mycds <- setOrderingFilter(mycds, var.genes)
  p2 <- plot_ordering_genes(mycds);p2
}
##使用monocle选择的高变基因
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
p3 <- plot_ordering_genes(mycds);p3
##结果对比
p3

#选择不同的基因集，拟时分析的结果不同，实践中可以几种方法都试一下。
#####降维及细胞排序
####使用disp.genes开展后续分析

#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#排序
mycds <- orderCells(mycds)


##可视化
#State轨迹分布图
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080")
plot1 <- plot_cell_trajectory(mycds, color_by = "State") + scale_color_manual(values = colour)+theme_bw()+
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size =30),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30),
        title = element_text(size = 30))+
  guides(color = guide_legend(override.aes = list(size = 12)))  ;plot1 ###控制图例的大小  
ggsave("State.pdf", plot = plot1, width = 8, height = 6)


##Cluster轨迹分布图
plot2 <- plot_cell_trajectory(mycds, color_by = "celltype")+ scale_color_manual(values = colour)+theme_bw()+
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size =30),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30),
        title = element_text(size = 30))+
  guides(color = guide_legend(override.aes = list(size = 12)));plot2
ggsave("Cluster.pdf", plot = plot2, width = 8, height = 6)


##Pseudotime轨迹图
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")+theme_bw()+
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size =30),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30),
        title = element_text(size = 30))+
  guides(color = guide_legend(override.aes = list(size = 12)));plot3 
ggsave("Pseudotime.pdf", plot = plot3, width =8, height = 6)

##合并作图
plotc <- plot1|plot2|plot3
ggsave("Combination.pdf", plot = plotc, width = 24, height = 6)
#ggsave("pseudotime/Combination.png", plot = plotc, width = 10, height = 3.5)
##保存结果
write.csv(pData(mycds), "pseudotime.csv")


####轨迹图分面显示
p1 <- plot_cell_trajectory(mycds, color_by = "State") + facet_wrap(~State, nrow = 1)+  scale_color_nejm() +theme_bw()+
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size =30),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30),
        title = element_text(size = 30),
        strip.text = element_text(face = "bold", color = "#7570B3",
                                  hjust = 0, size = 30),
        strip.background = element_rect(fill = "#E6AB02", linetype = "dotted"))+
  guides(color = guide_legend(override.aes = list(size = 12)));p1

p2 <- plot_cell_trajectory(mycds, color_by = "celltype") + facet_wrap(~celltype, nrow = 1)+ scale_color_manual(values = colour)+theme_bw()+
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size =30),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30),
        title = element_text(size = 30),
        strip.text = element_text(face = "bold", color = "#7570B3",
                                  hjust = 0, size = 30),
        strip.background = element_rect(fill = "#E6AB02", linetype = "dotted"))+
  guides(color = guide_legend(override.aes = list(size = 12)));p2


p3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime") +theme_bw()+
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size =30),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30),
        title = element_text(size = 30),
        strip.text = element_text(face = "bold", color = "#7570B3",
                                  hjust = 0, size = 30),
        strip.background = element_rect(fill = "#E6AB02", linetype = "dotted"))+
  guides(color = guide_legend(override.aes = list(size = 12)));p3

#密度图
library(ggpubr)
df <- pData(mycds) 
## pData(mycds)取出的是mycds对象中mycds@phenoData@data的内容
View(df)
p4 <- ggplot(df, aes(Pseudotime, colour = celltype, fill=celltype)) +
  geom_density(bw=1,size=2,alpha = 0.5)+theme_classic2()+theme_bw()+
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size =30),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30),
        title = element_text(size = 30),
        strip.text = element_text(face = "bold", color = "#7570B3",
                                  hjust = 0, size = 30),
        strip.background = element_rect(fill = "#E6AB02", linetype = "dotted"))+
  guides(color = guide_legend(override.aes = list(size = 12)));p4
plotc <- p2/p3/p4
ggsave("trajectory_facet.pdf", plot = plotc, width = 20, height = 20)


##########树形图
p2 <- plot_complex_cell_trajectory(mycds, color_by = "celltype") + facet_wrap(~celltype, nrow = 1)+ scale_color_manual(values = colour)+theme_bw()+
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size =30),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30),
        title = element_text(size = 30),
        strip.text = element_text(face = "bold", color = "#7570B3",
                                  hjust = 0, size = 30),
        strip.background = element_rect(fill = "#E6AB02", linetype = "dotted"))+
  guides(color = guide_legend(override.aes = list(size = 12)));p2


p3 <- plot_complex_cell_trajectory(mycds, color_by = "Pseudotime") +theme_bw()+
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size =30),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30),
        title = element_text(size = 30),
        strip.text = element_text(face = "bold", color = "#7570B3",
                                  hjust = 0, size = 30),
        strip.background = element_rect(fill = "#E6AB02", linetype = "dotted"))+
  guides(color = guide_legend(override.aes = list(size = 12)));p3

#密度图
library(ggpubr)
df <- pData(mycds) 
## pData(mycds)取出的是mycds对象中mycds@phenoData@data的内容
View(df)
p4 <- ggplot(df, aes(Pseudotime, colour = celltype, fill=celltype)) +
  geom_density(bw=1,size=2,alpha = 0.5)+theme_classic2()+theme_bw()+
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size =30),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30),
        title = element_text(size = 30),
        strip.text = element_text(face = "bold", color = "#7570B3",
                                  hjust = 0, size = 30),
        strip.background = element_rect(fill = "#E6AB02", linetype = "dotted"))+
  guides(color = guide_legend(override.aes = list(size = 12)));p4
plotc <- p2/p3/p4
ggsave("trajectory_facet_tree.pdf", plot = plotc, width = 20, height = 20)



###Monocle基因可视化
#s.genes <- c("Cd74","Fscn1","Gata3","Rora","Pclaf","Birc5","Xcl1","Ifng","Trac","Lef1")
s.genes <- c("Cd74","Gata3","Pclaf","Xcl1","Trac")
#p1 <- plot_genes_jitter(mycds[s.genes,], grouping = "State", color_by = "State")
p2 <- plot_genes_violin(mycds[s.genes,], grouping = "State", color_by = "State");p2
#p2 <- plot_genes_violin(mycds[s.genes,], grouping = "celltype", color_by = "celltype");p2
#p3 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = "State")
p3 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = "Pseudotime")
plotc <- p2|p3
ggsave("genes_visual.pdf", plot = plotc, width = 20, height = 15)


#拟时序展示单个基因表达量
colnames(pData(mycds))
pData(mycds)$Gata3 = log2( exprs(mycds)['Gata3',]+1)
p1=plot_cell_trajectory(mycds, color_by = "Gata3")  + scale_color_gsea()
pData(mycds)$Cd74 = log2(exprs(mycds)['Cd74',]+1)
p2=plot_cell_trajectory(mycds, color_by = "Cd74")    + scale_color_gsea()
pData(mycds)$Xcl1 = log2(exprs(mycds)['Xcl1',]+1)
p3=plot_cell_trajectory(mycds, color_by = "Xcl1")    + scale_color_gsea()
pData(mycds)$Trac = log2(exprs(mycds)['Trac',]+1)
p4=plot_cell_trajectory(mycds, color_by = "Trac")    + scale_color_gsea()
pData(mycds)$Pclaf = log2(exprs(mycds)['Pclaf',]+1)
p5=plot_cell_trajectory(mycds, color_by = "Pclaf")    + scale_color_gsea()

library(patchwork)
pc <- p1+p2+p3+p4+p5
ggsave("single_genes_visual.pdf", plot = pc, width = 15, height = 12)



###拟时相关基因聚类热图
##Monocle中differentialGeneTest()函数可以按条件进行差异分析，
#将相关参数设为fullModelFormulaStr = "~sm.ns(Pseudotime)"时，可以找到与拟时相关的差异基因。
#我们可以按一定的条件筛选基因后进行差异分析，全部基因都输入会耗费比较长的时间。
#建议使用cluster差异基因或高变基因输入函数计算。分析结果主要依据qval区分差异的显著性，
#筛选之后可以用plot_pseudotime_heatmap函数绘制成热图。

if (F) {
  #cluster差异基因
  #diff.wilcox <- FindAllMarkers(object = obj.combine, only.pos = TRUE, MIn.pct = 0.25,
  #                             thresh.use = 0.25)
  #diff.wilcox已经运行过来
  
  sig_diff.genes <- subset(diff.wilcox,p_val_adj<0.0001&abs(avg_log2FC)>0.75)$gene
  sig_diff.genes <- unique(as.character(sig_diff.genes))
  diff_test <- differentialGeneTest(mycds[sig_diff.genes,], cores = 1, 
                                    fullModelFormulaStr = "~sm.ns(Pseudotime)")
  sig_gene_names <- row.names(subset(diff_test, qval < 0.01))
  p1 = plot_pseudotime_heatmap(mycds[sig_gene_names,], num_clusters=5,
                               show_rownames=T, return_heatmap=T)
  ggsave("pseudotime_heatmap1.pdf", plot = p1, width = 10, height = 12)
}

#monocle高变基因
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1&dispersion_empirical >= 1*dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
diff_test <- differentialGeneTest(mycds[disp.genes,], cores = 4,fullModelFormulaStr = "~sm.ns(Pseudotime)")
#diff_test <- differentialGeneTest(mycds[disp.genes,], cores = 4,fullModelFormulaStr = "~celltype")
#~后面是表示对谁做差异分析的变量，理论上可以为p_data的任意列名
diff_test<- subset(diff_test, qval < 1e-02)   ###基因太多，可以根据qval < 1e-100标准筛选部分基因展示
diff_test <- diff_test[order(diff_test$qval,decreasing=F),] ##差异表达基因作为轨迹构建的基因,差异基因的选择标准是qval<0.01,decreasing=F表示按数值增加排序
##差异基因的结果文件保存
write.csv(diff_test,file="train.monocle.DEG.csv")
sig_gene_names <- row.names(diff_test)[1:50]   ###取100个基因展示

p2 = plot_pseudotime_heatmap(mycds[sig_gene_names,], num_clusters=5,
                             show_rownames=T, return_heatmap=T)
ggsave("pseudotime_heatmap.pdf", plot = p2, width = 8, height =8)

#################################################################################################
#################################################################################################
#################################################################################################
#前面通过设置num_clusters将热图聚成了5个cluster，如过想要把每个cluster的基因单独提出来富集做分析
p2$tree_row
# Call:
# hclust(d = d, method = method)
# Cluster method   : ward.D2 
# Number of objects: 2829 
clusters <- cutree(p2$tree_row, k = 5)
gene_group <- data.frame(clusters)
gene_group[,1] <- as.character(gene_group[,1])
colnames(gene_group) <- "Cluster"
gene_group$gene <- rownames(gene_group)
table(gene_group)
# 1    2    3    4 
# 570 1031  506  722 
write.csv(gene_group, "Time_clustering_all.csv", row.names = F)
library(clusterProfiler)
library(org.Mm.eg.db)
allcluster_go=data.frame()
for (i in unique(gene_group$Cluster)) {
  small_gene_group=filter(gene_group,gene_group$Cluster==i)
  df_name=bitr(small_gene_group$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
  go <- enrichGO(gene         = unique(df_name$ENTREZID),
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.2,
                 readable      = TRUE)
  go_res=go@result
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    allcluster_go=rbind(allcluster_go,go_res)
  }
}

head(allcluster_go[,c("ID","Description","qvalue","cluster")])

###################################################################################################
###################################################################################################
#########################绘制每一群的GO富集分析
for (i in unique(allcluster_go$cluster)) {
  subcluster <- subset(allcluster_go,allcluster_go$cluster==i)  #提取出第一群的cluster富集
  #提取结果前Top20绘图(或自定义所需pathway绘图)：
  top20 <- subcluster[1:10,]
  #指定绘图顺序（转换为因子）：
  top20$pathway <- factor(top20$Description,levels = rev(top20$Description))
  
  ##计算Rich Factor（富集因子）：
  top20 <- mutate(top20,
                  RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
  ##计算Fold Enrichment（富集倍数）：
  top20$FoldEnrichment <- apply(top20,1,function(x){
    GeneRatio <- eval(parse(text = x["GeneRatio"]))
    BgRatio <- eval(parse(text = x["BgRatio"]))
    foldEnrichment <- round(GeneRatio/BgRatio,2)
    foldEnrichment
  })
  head(top20$RichFactor)
  head(top20$FoldEnrichment)
  colnames(top20)
  ######
  #1.常规画法：
  #自定义主题：
  mytheme <- theme(axis.title = element_text(size = 13),
                   axis.text = element_text(size = 11),
                   plot.title = element_text(size = 14,
                                             hjust = 0.5,
                                             face = "bold"),
                   legend.title = element_text(size = 13),
                   legend.text = element_text(size = 11))
  ######
  #2.论文画法：
  mytheme2 <- mytheme + theme(axis.text.y = element_blank()) #先在自定义主题中隐去y轴文本标签显示
  top20$text_x <- rep(0,10) #新增一列重复数组，使绘图时文本标签能从固定位置开始；
  p2 <- ggplot(data = top20,
               aes(x = RichFactor, y = pathway)) +
    geom_bar(aes(fill = -log10(pvalue)), stat = "identity", width = 0.8, alpha = 0.7) +
    scale_fill_distiller(palette = "YlOrRd", direction = 1) +
    labs(x = "RichFactor", y = "pathway", title = "Go BP enrichment barplot") +
    geom_text(aes(x = text_x, #用新增的重复数组控制文本标签起始位置
                  label = pathway),
              hjust = 0)+ #hjust=0，左对齐
    theme_classic() + mytheme2
  p2
  ggsave(paste0("goenrich_cluster",i,".pdf"),p2,width = 8,height = 3)
}

#######################################################################################################
#######################################################################################################
########拟时序的top6基因展示
keygenes <- head(sig_gene_names,5)
mycds_subset <- mycds[keygenes,]
##可视化：以state/celltype/pseudotime进行
#p1 <- plot_genes_in_pseudotime(mycds_subset, color_by = "State")
p2 <- plot_genes_in_pseudotime(mycds_subset, color_by = "celltype")
p3 <- plot_genes_in_pseudotime(mycds_subset, color_by = "Pseudotime")
plotc <- p2|p3
ggsave("Genes_top5_pseudotimeplot.pdf", plot = plotc, width = 16, height = 8)




##BEAM分析
##单细胞轨迹中通常包括分支，它们的出现是因为细胞的表达模式不同。
#当细胞做出命运选择时，或者遗传、化学或环境扰动时，就会表现出不同的基因表达模式。
#BEAM(Branched expression analysis modeling)是一种统计方法，
#用于寻找以依赖于分支的方式调控的基因。
#monocle高变基因
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
mycds_sub <- mycds[disp.genes,]
plot_cell_trajectory(mycds_sub, color_by = "State")
beam_res <- BEAM(mycds_sub, branch_point = 1, cores = 8)
beam_res <- subset(beam_res,qval<1e-02)
beam_res <- beam_res[order(beam_res$qval),]   ####默认的升序排列decreasing=F
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
write.csv(beam_res,file = "beam_res.csv")
mycds_sub_beam <- mycds_sub[row.names(beam_res)[1:50],]  #提取beam前top50基因画图


p <- plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = 1, num_clusters = 5, show_rownames = T,return_heatmap = T)
ggsave("BEAM_heatmap.pdf",p$ph_res,width = 8,height = 12)

#####可以提取基因进行后续富集分析，步骤同上，此步略，参考https://www.jianshu.com/p/5d6fd4561bc0
##参考https://blog.csdn.net/qq_38774801/article/details/117631707
gene_group=p$annotation_row
gene_group$gene=rownames(gene_group)
table(gene_group)

#显著差异基因(top100)按热图结果排序并保存
#genes <- rownames(head(beam_res,3))
genes <- c("Cd74","Gata3","Pclaf","Xcl1","Trac")
pdf("genes_branched_pseudotime_5.pdf",width = 9,height = 6)
plot_genes_branched_pseudotime(mycds[genes,],
                               branch_point = 1,
                               color_by = "celltype",
                               ncol = 2)
dev.off()

genes <- rownames(head(beam_res,4))
pdf("genes_branched_pseudotime_top4.pdf",width = 9,height = 4)
plot_genes_branched_pseudotime(mycds[genes,],
                               branch_point = 1,
                               color_by = "Pseudotime",
                               ncol = 2)
dev.off()

save(mycds,file = "monocle.RData")
