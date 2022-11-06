#多样本整合
#归一化
obj.list <- lapply(obj.list,function(x) {
 NormalizeData(x)
})

#寻找高变异度基因
obj.list <- lapply(obj.list, function(x) {
  FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
  })
obj.list
#整合成一个对象
#寻找锚点
obj.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:20)
#根据锚点来整合
obj.combined <- IntegrateData(anchorset = obj.anchors, dims = 1:20)
DefaultAssay(obj.combined) <- "integrated"   #更改默认数组
#对整合后的数据进行尺度变换
all.genes <- rownames(obj.combined[["RNA"]]@data)    #针对所有基因
length(all.genes)
obj.combined <- ScaleData(obj.combined, features = all.genes)
