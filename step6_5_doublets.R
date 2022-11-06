
setwd( "/home/data/t040243/CopyOfjiangshijiu")
dir.create("Doublet")
setwd("Doublet")
##devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
##########################################################################################################################
##实际分析中主要有以下四步：

##使用单细胞数据创建seurat对象，并进行数据标准化、降维，最好进行聚类和细胞类型鉴定；
##使用BCmvn（均值-方差标准化双峰系数）寻找计算pANN的最优pK值；
##根据泊松分布统计原理估计样本中doublets的数量，并排除DoubletFinder不能检出的同源doublets，得到优化后的预估doublets数量；
##使用前两步得到的参数运行函数鉴定doublets数据。

#这是一个测试最佳参数的过程，运行速度慢
load("~/CopyOfjiangshijiu/cluster/rename_cluster.RData")
scRNA_harmony <- obj.combined      ###最好选注释后有celltype的数据
DefaultAssay(scRNA_harmony) <- "RNA"  ##如果是intergrated数据可以用这一步
scRNA_harmony <- NormalizeData(scRNA_harmony) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)


sweep.res.list <- paramSweep_v3(scRNA_harmony, PCs = 1:10, sct = F)
#使用log标准化，sct参数设置为 sct = F（默认 ）,如使用SCT标准化方法，设置为T
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats) #可以看到最佳参数的点
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() #提取最佳pk值

## 排除不能检出的同源doublets，优化期望的doublets数量
#DoubletRate = 0.076                     # 5000细胞对应的doublets rate是7.6%
DoubletRate = ncol(scRNA_harmony)*8*1e-6 #更通用
#估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞
scRNA_harmony@meta.data$celltype <-scRNA_harmony@active.ident 
homotypic.prop <- modelHomotypic(scRNA_harmony$celltype) 
# 计算双细胞比例
nExp_poi <- round(DoubletRate*ncol(scRNA_harmony)) 
# 使用同源双细胞比例对计算的双细胞比例进行校正 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## 使用确定好的参数鉴定doublets
scRNA_harmony <- doubletFinder_v3(scRNA_harmony, PCs = 1:10, pN = 0.25, pK = pK_bcmvn, 
                                  nExp = nExp_poi.adj, reuse.pANN = F, sct = T)

## 结果展示，分类结果在scRNA_harmony@meta.data中
colnames(scRNA_harmony@meta.data)
plot1=DimPlot(scRNA_harmony, reduction = "umap", group.by = "DF.classifications_0.25_0.29_45")
plot2=DimPlot(scRNA_harmony, reduction = "umap",label = T) 
#combinate
plotc <- plot1+plot2
pdf("doublets.pdf",plotc,width = 10,height = 5)
plotc
dev.off()
##seu  这是一个经过充分处理的Seurat对象（即，在NormalizeData，FindVariableGenes，ScaleData，RunPCA和RunTSNE全部运行之后）。 
#pN  定义生成的人工双峰的数量 (variable numbers of artificial doublets)，表示为合并的真实人工数据的一部分。基于DoubletFinder在很大程度上是pN不变，默认设置为25％ 
#pK  定义用于计算pANN (proportion of artificial nearest neighbors) 的PC邻域大小，表示为合并的真实人工数据的一部分。没有设置默认值，应该根据每个scRNA-seq数据集调整pK。 
#nExp  定义用于进行最终双峰/单峰预测的pANN阈值。可以从10X / Drop-Seq中的细胞密度估计该值，并根据同型双峰的估计比例进行调整。



Cells.singlet <- subset(scRNA_harmony@meta.data, DF.classifications_0.25_0.29_45=="Singlet")
scRNA_harmony_singlet <- subset(scRNA_harmony, cells=row.names(Cells.singlet))
plot_singlet=DimPlot(scRNA_harmony_singlet, reduction = "umap",split.by="group",label = T) ;plot_singlet

Cells.doublet <- subset(scRNA_harmony@meta.data, DF.classifications_0.25_0.29_45=="Doublet")
scRNA_harmony_doublet <- subset(scRNA_harmony, cells=row.names(Cells.doublet))
plot_doublet=DimPlot(scRNA_harmony_doublet, reduction = "umap",label = T) ;plot_doublet



#markers <- FindAllMarkers(object = scRNA_harmony_singlet, test.use="wilcox" ,
#                          only.pos = TRUE,
#                          logfc.threshold = 0.25)   
#all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
#top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#write.csv(all.markers,file = "scRNA_harmony_diffMarkers.csv")
#write.csv(top10,file = "scRNA_harmony_diffMarkers_top10.csv")

###此处应该将scRNA_harmony_singlet再用scTYPE跑一遍，确定celltype
obj.combined <- scRNA_harmony_singlet
sel.clust = "integrated_snn_res.0.1"
obj.combined <- SetIdent(obj.combined, value = sel.clust)
if(T){  #########去除双细胞后重新经过scTYPE命名 
  new.cluster.ids <- c("ILC2a","ILCdc","ILC2b","ILC1",
                       "ILCt","ILCdc",  "ILCdc")
  names(new.cluster.ids) <- levels(obj.combined)
  obj.combined <- RenameIdents(obj.combined, new.cluster.ids)
  table(obj.combined@active.ident)
}


sample.markers <- FindAllMarkers(obj.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- sample.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#保存每个cluster top10的marker基因
top10_table=unstack(top10, gene ~ cluster)
names(top10_table)=gsub("X","cluster",names(top10_table))
write.csv(file="top10_marker_genes.csv",top10_table,row.names=F)
write.csv(file="sample_marker_genes.csv",sample.markers,row.names=F)


obj.combined@meta.data[["celltype"]] <- obj.combined@active.ident
save(obj.combined,scRNA_harmony,file = "after_doublet_obj.combined.RData")


##########################################################################################
##########################################################################################
####正式出版图Marker   正式出版图Marker   正式出版图Marker
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
  guides(color = guide_legend(override.aes = list(size = 12)))  ###控制图例的大小
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



##########################################################################################
##########################################################################################
###正式出版图Marker
ILC=c("Ptprc","Il7r","Id2")
Lin = c("Itgax","Itgam","Cd8a","Cd3e","Cd3d","Cd4","Cd19","Cma1","Ly6c1","Adgre1","Pdgfra")
ILC2a=c("Klrg1","Gata3","Rora","Il5","Areg","Hs3st1")
ILC2b=c("Pclaf","Birc5","Stmn1")
#ILC1=c("Xcl1","Ctsw","Ccl4","Ifng","Klrc1")
#ILC1=c("Tbx21", "Eomes", "Ifng","Tnf","Prf1","Gzma")
ILC1=c("Ifng","Tnf","Tbx21","Klrc1","Xcl1")
ILCdc=c("Cd74","H2-Eb1","H2-Aa","H2-Ab1","Cd83")
ILCt=c("Trac","Cd2","Lef1","Klf2","Cd3g")
ILC3=c("Rorc","Il22","Il17a")   ####"Ahr","Ccr6"

pdf('ILC_markers.pdf',width = 12,height = 8)
DotPlot(obj.combined, assay = "RNA", 
        features = c(ILC,Lin,ILC3,ILC2a,ILC2b,ILCdc,ILC1,ILCt), 
        cols = c("blue", "red","green"), split.by = "group") + 
  coord_flip()+ RotatedAxis()
dev.off()

#####################
######heatmap
#top10 <- sample.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 <- sample.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
p=DoHeatmap(obj.combined, features = top5$gene, 
            size = 6,assay="RNA",group.bar.height = 0.05);p
ggsave(filename='top5.markers_heatmap.pdf',width = 12,height = 8)

#########################
#####Feature图
p1 <- FeaturePlot(obj.combined,features = c("Gata3","Pclaf","Cd74","Xcl1","Trac"),cols = c("lightgrey","red"),reduction = "umap")
p1
ggsave("Featureplot.pdf", plot = p1, width = 12, height = 8)
p2 <- FeaturePlot(obj.combined,features = c("Gata3"),pt.size = 0.1,cols = c("lightgrey","red"),reduction = "umap")
p2
ggsave("Featureplot_Gata3.pdf", plot = p2, width = 4, height = 3)
p3 <- FeaturePlot(obj.combined,features = c("Pclaf"),pt.size = 0.1,cols = c("lightgrey","red"),reduction = "umap")
p3
ggsave("Featureplot_Pclaf.pdf", plot = p3, width = 4, height = 3)
p4 <- FeaturePlot(obj.combined,features = c("Cd74"),pt.size = 0.1,cols = c("lightgrey","red"),reduction = "umap")
p4
ggsave("Featureplot_Cd74.pdf", plot = p4, width = 4, height = 3)
p5 <- FeaturePlot(obj.combined,features = c("Xcl1"),pt.size = 0.1,cols = c("lightgrey","red"),reduction = "umap")
p5
ggsave("Featureplot_Xcl1.pdf", plot = p5, width = 4, height = 3)
p6 <- FeaturePlot(obj.combined,features = c("Trac"),pt.size = 0.1,cols = c("lightgrey","red"),reduction = "umap")
p6
ggsave("Featureplot_Trac.pdf", plot = p6, width = 4, height = 3)
p7 <- FeaturePlot(obj.combined,features = c("Mki67"),pt.size = 0.1,cols = c("lightgrey","red"),reduction = "umap")
p7
ggsave("Featureplot_Mki67.pdf", plot = p7, width = 4, height = 3)
p8 <- FeaturePlot(obj.combined,features = c("Pcna"),pt.size = 0.1,cols = c("lightgrey","red"),reduction = "umap")
p8
ggsave("Featureplot_Pcna.pdf", plot = p8, width = 4, height = 3)
p9 <- FeaturePlot(obj.combined,features = c("Top2a"),pt.size = 0.1,cols = c("lightgrey","red"),reduction = "umap")
p9
ggsave("Featureplot_Top2a.pdf", plot = p9, width = 4, height = 3)
p10 <- FeaturePlot(obj.combined,features = c("Ifng"),pt.size = 0.1,cols = c("lightgrey","red"),reduction = "umap")
p10
ggsave("Featureplot_Ifng.pdf", plot = p10, width = 4, height = 3)

##########################################################################
######堆叠小提琴图
pdf('stack_Vlnplot_markers.pdf',width = 24,height = 8)
VlnPlot(obj.combined,features = c(ILC,ILC1,ILC2a,ILC2b,ILCdc,ILCt),stack = T)
dev.off()

##############ILC特异性小提琴图
pdf('Vlnplot_Zbtb46.pdf',width = 12,height = 8)
VlnPlot(obj.combined,features = "Zbtb46")
dev.off()
pdf('Vlnplot_Batf3.pdf',width = 12,height = 8)
VlnPlot(obj.combined,features = "Batf3")
dev.off()
pdf('Vlnplot_Id3.pdf',width = 12,height = 8)
VlnPlot(obj.combined,features = "Id3")
dev.off()
pdf('Vlnplot_Il1b.pdf',width = 12,height = 8)
VlnPlot(obj.combined,features = "Il1b")
dev.off()
pdf('Vlnplot_Lef1.pdf',width = 12,height = 8)
VlnPlot(obj.combined,features = "Lef1")
dev.off()
pdf('Vlnplot_S1pr1.pdf',width = 12,height = 8)
VlnPlot(obj.combined,features = "S1pr1")
dev.off()
pdf('Vlnplot_Trac.pdf',width = 12,height = 8)
VlnPlot(obj.combined,features = "Trac")
dev.off()
pdf('Vlnplot_Xcl1.pdf',width = 12,height = 8)
VlnPlot(obj.combined,features = "Xcl1")
dev.off()
pdf('Vlnplot_Ifng.pdf',width = 12,height = 8)
VlnPlot(obj.combined,features = "Ifng")
dev.off()
pdf('Vlnplot_Tbx21.pdf',width = 12,height = 8)
VlnPlot(obj.combined,features = "Tbx21")
dev.off()
pdf('Vlnplot_Hs3st1.pdf',width = 12,height = 8)
VlnPlot(obj.combined,features = "Hs3st1")
dev.off()
pdf('Vlnplot_Il5.pdf',width = 12,height = 8)
VlnPlot(obj.combined,features = "Il5")
dev.off()
pdf('Vlnplot_Areg.pdf',width = 12,height = 8)
VlnPlot(obj.combined,features = "Areg")
dev.off()
pdf('Vlnplot_Gata3.pdf',width = 12,height = 8)
VlnPlot(obj.combined,features = "Gata3")
dev.off()
pdf('Vlnplot_Pclaf.pdf',width = 12,height = 8)
VlnPlot(obj.combined,features = "Pclaf")
dev.off()
pdf('Vlnplot_Birc5.pdf',width = 12,height = 8)
VlnPlot(obj.combined,features = "Birc5")
dev.off()
pdf('Vlnplot_Bcl2a1b.pdf',width = 12,height = 8)
VlnPlot(obj.combined,features = "Bcl2a1b")
dev.off()





###########################################################################################
##经典出图marker基因   经典出图marker基因   经典出图marker基因   经典出图marker基因
ILC=c("Ptprc","Il7r")
ILC2=c("Id2","Rora" ,"Gata3","Bcl11b","Ets1", "Il5","Il13","Areg","Il1rl1","Hs3st1","Pdcd1","Il17rb","Klrg1"
       ,"Pclaf","Birc5","Stmn1")
ILC1=c("Tbx21", "Ifng","Klrc1","Klrb1c","Ccl3","Il21r","Ccl5","Ccl4"," Xcl1","Ctsw","Nkg7","Klf2")
ILCdc=c("Cd74","H2-Eb1","H2-Aa","H2-Ab1","Ccr7","Irf8","Batf3","Zbtb46","Cd86","Cd83","Fscn1","Il1b","Irf5","Id3","Sox4","Bcl2a1b")
ILCt=c("Trac","Gm8369","S1pr1","Cd2","Ms4a4b","Cd3g","Tcf7","Lef1")   ##"Trbc2"
Lin <- c("Itgax","Itgam","Cd8a","Cd3e","Cd3d","Cd4","Cd19","Cma1","Ly6c1","Adgre1","Pdgfra")

pdf('umap_markers.pdf',width = 12,height = 8)
DotPlot(obj.combined, assay = "RNA", 
        features = c(ILC,Lin,ILC2,ILCdc,ILC1,ILCt), 
        cols = c("blue", "red","green"), split.by = "group") + 
        coord_flip()+ RotatedAxis()
dev.off()

######################################################################################
########################堆叠小提琴图   堆叠小提琴图   堆叠小提琴图   堆叠小提琴图
ILC=c("Ptprc","Il7r")
ILC2=c("Rora" ,"Gata3","Bcl11b", "Il5","Il13","Areg","Il1rl1","Hs3st","Pdcd1","Il17rb","Klrg1"
       ,"Pclaf","Birc5","Stmn1")
ILC1=c("Tbx21", "Ifng","Klrc1","Klrb1c","Ccl3","Il21r","Ccl4"," Xcl1","Ctsw","Nkg7")
ILCdc=c("Cd74","Batf3","Zbtb46","Cd86","Cd83","Fscn1","H2-Eb1","H2-Aa","Il1b","Ccr7","Irf5","Id3","Sox4")
ILCt=c("Trac","Lef1","S1pr1","Cd3g","Tcf7")

pdf('stack_Vlnplot_markers.pdf',width = 12,height = 8)
VlnPlot(obj.combined,features = c(ILC,ILC1,ILC2,ILCdc,ILCt),stack = T)
dev.off()
#########################################################################


genes2=c("Zbtb16","Id2","Ly6a","Tcf7","Tox","Nfil3",   #ILCP,CLP,EILPs
         "Rora" ,"Gata3", "Bcl11b","Gfi1","Ets1",# ILC2
         "Il5","Il13","Areg","Il4","Il9", "Bmp7",# ILC2
         "Il1rl1", #The IL-33 receptor IL1RL1 (T1/ST2)
         "Klrg1", #ILC2
         "Thy1",  #CD90
         "Il10ra","Vipr2","Icos","Tnfsf4",  #ILC2其他的受体
         "Il17rb","Klrg1",  #IL25和IL17B受体，判断iILC2
         "Arg1",
         "Ptgdr2", # Prostaglandin D2 Receptor 2, CRTH2, 人ILC2表达
         "Tbx21", "Eomes", #ILC1细胞
         "Ifng","Tnf","Prf1","Gzma", #ILC1细胞
         "Klrb1c",  #NK1.1
         "Il22","Il17a", # ILC3细胞
         "Rorc","Ahr", #ILC3
         "Ncr1", "Ccr6", #ILC3细胞再分类用
         "Kit",   
         "Il2ra",   #CD25
         "Mki67","Top2a",  #细胞增殖相关基因
         "Cd68", "Cd163", "Cd14"  #Monocytes and macrophages,有待进一步认定
)
