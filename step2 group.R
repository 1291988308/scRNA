##分别读入每个样本的原始数据
setwd("/home/data/t040243/CopyOfjiangshijiu")
set.seed(123)  #设置随机数种子，使结果可重复
dir = c('data/MI/', 
        'data/MIRI/',
        'data/Sham/')
names(dir) = c('MI', 'MIRI', 'Sham')

obj.list <- list()
#以下代码会把每个样本的数据创建一个seurat对象，并存放到列表scRNAlist里
for(i in 1:length(dir))
{
  counts <- Read10X(data.dir = dir[i])
  obj.list[[i]] <- CreateSeuratObject(counts, min.cells=3,min.features = 200)
}     #保留至少在三个细胞里面表达的基因 ;#保留至少表达200个基因的细胞
names(obj.list) <- names(dir)
obj.list
obj.list$MI@meta.data$group <- 'MI'
obj.list$MIRI@meta.data$group <- 'MIRI'
obj.list$Sham@meta.data$group <- 'Sham'



#########################################################下面不用运行
