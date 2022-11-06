#质控QC，计算每种特征在每个细胞里面的百分比
dir.create("QC")
setwd("QC")
for (x in names(obj.list)){
obj.list[[x]][["percent.MT"]] <- PercentageFeatureSet(obj.list[[x]], pattern = "^mt-")
obj.list[[x]][["percent.RP"]] <- PercentageFeatureSet(obj.list[[x]], pattern = "^Rp[sl]")
obj.list[[x]][["percent.HB"]] <- PercentageFeatureSet(obj.list[[x]], pattern = "^Hb[^(p)]")
}

#循环绘制每个样本的QC图
qc_feature <- c("nFeature_RNA", "nCount_RNA", "percent.HB", "percent.MT",  "percent.RP")
for(x in names(obj.list)){
pdf(file=paste0("1_",x,"_quality_control.pdf"),width = 15,height=7)
print(VlnPlot(obj.list[[x]], features = qc_feature, ncol = 5, pt.size = 0.5))
dev.off()
}

#循环绘制特征之间的相互关系图
for(x in names(obj.list)){
pdf(file=paste0("1_",x,"_feature_relationship.pdf"),width = 12,height=10)
p1 <- FeatureScatter(obj.list[[x]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p2 <- FeatureScatter(obj.list[[x]], feature1 = "nCount_RNA", feature2 = "percent.HB")
p3 <- FeatureScatter(obj.list[[x]], feature1 = "nCount_RNA", feature2 = "percent.MT")
p4 <- FeatureScatter(obj.list[[x]], feature1 = "nCount_RNA", feature2 = "percent.RP")
print(p1 + p2 + p3 + p4)
dev.off()
}

#过滤细胞
obj.list <- lapply(obj.list, function(x) {
subset(x, subset = nFeature_RNA > 200 & 
            nFeature_RNA < 6500 &
           percent.MT < 20)})
obj.list


###阳选cd127和cd45,剔除cd11c和cd11b
for(x in names(obj.list)){
  sc <- obj.list[[x]]
  count=sc@assays$RNA@counts  
  ###阳选
  sc[["cd127"]] <- count["Il7r",]
  sc[["cd45"]] <- count["Ptprc",]
  sc[["cd11c"]] <- count["Itgax",]  #DC CD11c
  sc[["cd11b"]] <- count["Itgam",]  #髓系细胞 CD11b
  sc[["cd8a"]] <- count["Cd8a",]    #T细胞
  sc[["cd3e"]] <- count["Cd3e",]    #T细胞   
  sc[["cd19"]] <- count["Cd19",]    #B细胞  
  sc[["cma1"]] <- count["Cma1",]   #肥大细胞
  sc[["ly6c1"]] <- count["Ly6c1",]   #单核细胞
  sc[["adgre1"]] <- count["Adgre1",]   #巨噬细胞 F4/80
  sc[["Pdgfra"]] <- count["Pdgfra",]                #成纤维细胞
  sc[["Cd3d"]] <- count["Cd3d",]                #T细胞
  #sc[["ly6g"]] <- count["Ly6g",]   # 中性粒细胞（Gr1）
  #sc[["Cd3"]] <- count["Cd3",]
  #sc[["Cd45r"]] <- count["",]
  sc[["Cd4"]] <- count["Cd4",]
  sc_dp <- subset(sc, subset = cd127 >0 & cd45 >0 & cd11c<1 & cd11b<1 
                  &cd8a<1&cd3e<1&cd19<1&cma1<1&ly6c1<1&adgre1<1&Pdgfra<1&Cd3d<1
                  &Cd4<1) 
  obj.list[[x]]=sc_dp
}

#####################################################
#循环绘制QC后每个样本的QC图
qc_feature <- c("nFeature_RNA", "nCount_RNA", "percent.HB", "percent.MT",  "percent.RP")
for(x in names(obj.list)){
  pdf(file=paste0("2_afterQC",x,"_quality_control.pdf"),width = 15,height=7)
  print(VlnPlot(obj.list[[x]], features = qc_feature, ncol = 5, pt.size = 0.5))
  dev.off()
}

#循环绘制QC后特征之间的相互关系图
for(x in names(obj.list)){
  pdf(file=paste0("2_afterQC",x,"_feature_relationship.pdf"),width = 12,height=10)
  p1 <- FeatureScatter(obj.list[[x]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p2 <- FeatureScatter(obj.list[[x]], feature1 = "nCount_RNA", feature2 = "percent.HB")
  p3 <- FeatureScatter(obj.list[[x]], feature1 = "nCount_RNA", feature2 = "percent.MT")
  p4 <- FeatureScatter(obj.list[[x]], feature1 = "nCount_RNA", feature2 = "percent.RP")
  print(p1 + p2 + p3 + p4)
  dev.off()
}




