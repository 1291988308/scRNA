#2. 读入数据并查看分类变量
rm(list=ls())
setwd("/home/data/t040243/CopyOfjiangshijiu")
dir.create("Figure")
setwd("Figure")
load("/home/data/t040243/CopyOfjiangshijiu/Doublet/after_doublet_obj.combined.RData")
Idents(obj.combined) <- obj.combined@meta.data[["group"]]
obj.combined_Sham <- subset(obj.combined,idents="Sham")
scRNA <- obj.combined_Sham

unique(scRNA$group)
## [1] "LN"       "BRONCHO"  "NS"       "EBUS"     "EFFUSION" "LUNG"
table(scRNA$celltype)
## 
##     B lymphocytes Endothelial cells  Epithelial cells       Fibroblasts 
##              7403               462             10241              1142 
##        MAST cells     Myeloid cells          NK cells  Oligodendrocytes 
##               785             12016              3397               457 
##     T lymphocytes      UndeterShamned 
##             21896               201
#3. 饼图
table(scRNA$celltype)
## 
##     B lymphocytes Endothelial cells  Epithelial cells       Fibroblasts 
##              7403               462             10241              1142 
##        MAST cells     Myeloid cells          NK cells  Oligodendrocytes 
##               785             12016              3397               457 
##     T lymphocytes      UndeterShamned 
##             21896               201

#BiocManager::install("plotrix")
library(plotrix)
library(dplyr)
## 
## Attaching package: 'dplyr'
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
library(ggsci)
mynames <-   table(scRNA$celltype) %>% names()
myratio <-  table(scRNA$celltype) %>% as.numeric()
pielabel <- paste0(mynames," (", round(myratio/sum(myratio)*100,2), "%)")

cols <-pal_npg("nrc")(10)#
cols
##  [1] "#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF" "#8491B4FF"
##  [7] "#91D1C2FF" "#DC0000FF" "#7E6148FF" "#B09C85FF"
pdf("pie_Sham.pdf",width = 12,height = 8)
pie(myratio, labels=pielabel,
    radius = 1.0,clockwise=T,
    main = "celltype_Sham",col = cols)
dev.off()


# 绘制 3D 图，faShamly 要设置你系统支持的中文字体库
pdf("3Dpie_Sham.pdf",width = 8,height = 6)
pie3D(myratio,labels = pielabel,explode = 0.1, radius=1,theta=pi/4,
      main = "ILCs Proportion of Sham",
      height = 0.3,
      labelcex  = 2)
dev.off()

#4. 甜甜圈图
# 并没有直接画甜甜圈图的R包，所以在饼图源代码的基础上改改
doughnut <- function (x, labels = names(x), edges = 200, outer.radius = 0.8,
                      inner.radius=0.6, clockwise = FALSE,
                      init.angle = if (clockwise) 90 else 0, density = NULL,
                      angle = 45, col = NULL, border = FALSE, lty = NULL,
                      main = NULL, ...)
{
  if (!is.numeric(x) || any(is.na(x) | x < 0))
    stop("'x' values must be positive.")
  if (is.null(labels))
    labels <- as.character(seq_along(x))
  else labels <- as.graphicsAnnot(labels)
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
  plot.new()
  pin <- par("pin")
  xlim <- ylim <- c(-1, 1)
  if (pin[1L] > pin[2L])
    xlim <- (pin[1L]/pin[2L]) * xlim
  else ylim <- (pin[2L]/pin[1L]) * ylim
  plot.window(xlim, ylim, "", asp = 1)
  if (is.null(col))
    col <- if (is.null(density))
      palette()
  else par("fg")
  col <- rep(col, length.out = nx)
  border <- rep(border, length.out = nx)
  lty <- rep(lty, length.out = nx)
  angle <- rep(angle, length.out = nx)
  density <- rep(density, length.out = nx)
  twopi <- if (clockwise)
    -2 * pi
  else 2 * pi
  t2xy <- function(t, radius) {
    t2p <- twopi * t + init.angle * pi/180
    list(x = radius * cos(t2p),
         y = radius * sin(t2p))
  }
  for (i in 1L:nx) {
    n <- max(2, floor(edges * dx[i]))
    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n),
              outer.radius)
    polygon(c(P$x, 0), c(P$y, 0), density = density[i],
            angle = angle[i], border = border[i],
            col = col[i], lty = lty[i])
    Pout <- t2xy(mean(x[i + 0:1]), outer.radius)
    lab <- as.character(labels[i])
    if (!is.na(lab) && nzchar(lab)) {
      lines(c(1, 1.05) * Pout$x, c(1, 1.05) * Pout$y)
      text(1.1 * Pout$x, 1.1 * Pout$y, labels[i],
           xpd = TRUE, adj = ifelse(Pout$x < 0, 1, 0),
           ...)
    }      
    Pin <- t2xy(seq.int(0, 1, length.out = n*nx),
                inner.radius)
    polygon(Pin$x, Pin$y, density = density[i],
            angle = angle[i], border = border[i],
            col = "white", lty = lty[i])
  }
  
  title(main = main, ...)
  invisible(NULL)
}

# 绘图
df <- table(scRNA$celltype) %>% as.data.frame()
labs <- paste0(df$Var1," (", round(df$Freq/sum(df$Freq)*100,2), "%)")


pdf("circle_Sham.pdf",width = 12,height = 8)
p.circle <- doughnut(
  df$Freq,
  labels=labs, 
  init.angle=90,     # 设置初始角度
  col = cols , # 设置颜色 
  border="white",    # 边框颜色 
  inner.radius= 0.6, # 内环大小
  outer.radius = 1,  #外环大小
  cex = 1.5,           # 字体大小
  main = "ILC Proportion of Sham")          
dev.off()
## Warning in rep(lty, length.out = nx): 'x' is NULL so the result will be NULL
## Warning in rep(density, length.out = nx): 'x' is NULL so the result will be NULL


p.circle
## NULL
#5. 堆积柱状图
library(ggplot2)
scRNA <- obj.combined
cellnum <- table(scRNA$celltype,scRNA$group);cellnum 
cell.prop<-as.data.frame(prop.table(cellnum));cell.prop
colnames(cell.prop)<-c("celltype","group","Proportion")



p.bar <- ggplot(cell.prop,aes(group,Proportion,fill=celltype))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values=cols[1:10])+#自定义fill的颜色
  ggtitle("cell Proportion")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(size=20,title="celltype"))+
  theme(axis.text = element_text(size = 30),
    axis.title = element_text(size =30),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30),
        title = element_text(size = 30))+
  guides(color = guide_legend(override.aes = list(size = 6)))

pdf("stack_bar_3group.pdf",width = 12,height = 8)
p.bar
dev.off()


#6. 箱线图
scRNA <- obj.combined
cellnum <- table(scRNA$orig.ident,scRNA$celltype)
cellnum[1:3,1:5]
##             
##              B lymphocytes Endothelial cells Epithelial cells Fibroblasts
##   BRONCHO_11           421                 0               67           2
##   BRONCHO_58            51                 0              190          37
##   EBUS_06               52                 0              520           2
##   EBUS_10              435                 0               70           1
for (i in 1:nrow(cellnum)) {
  cellnum[i,] <- cellnum[i,]/sum(cellnum[i,])  
}
cellnum <- as.data.frame(cellnum)

library(reshape2)
colnames(cellnum) <- c('Sample','celltype','Freq')

for(i in 1:nrow(cellnum)){
  cellnum$group[i] <- strsplit(as.character(cellnum$Sample[i]),'_')[[1]][1]
}
unique(cellnum$group)
## [1] "BRONCHO"  "EBUS"     "EFFUSION" "LN"       "LUNG"     "NS"
library(tidyverse)
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──
## ✔ tibble  3.1.7     ✔ purrr   0.3.4
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()

cellnum$group <- factor(cellnum$group,levels = c("MI","MIRI","Sham" ))

library(ggplot2)
library(ggsci)
library(ggsignif)
unique(cellnum$group)
## [1] BRONCHO  EBUS     EFFUSION LN       LUNG     NS      
## Levels: LUNG BRONCHO EBUS EFFUSION LN NS
compaired <- list(c('MI','Sham'),
                  c('MIRI',"Sham"),
                  c("MI","MIRI"))

myplot<- list()
for (i in 1:length(unique(cellnum$celltype))) {
  myplot[[i]] <-ggplot(data=cellnum[cellnum$celltype==unique(cellnum$celltype)[i],],
                       aes(x=group,y=Freq,fill=group))+
    geom_violin()+
    geom_boxplot(width=0.2,
                 position = position_dodge(0.9))+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1))+
    scale_fill_manual(values = pal_npg("nrc")(10))+facet_wrap(~celltype)+
    # stat_compare_means(aes(group=group))+
    geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = 't.test')#t.test, wilcox.test
}

p.box <- Seurat::CombinePlots(myplot,ncol = 3,legend='right')
## Warning: CombinePlots is being deprecated. Plots should now be combined using
## the patchwork system.
## Warning: Computation failed in `stat_signif()`:
## missing value where TRUE/FALSE needed
## Computation failed in `stat_signif()`:
## missing value where TRUE/FALSE needed
pdf("boxplot_celltype.pdf",width = 12,height = 8)
p.box
dev.off()



if(F){
#BiocManager::install("ggpubr")
library("ggpubr")
####
for (i in 1:length(unique(cellnum$celltype))) {
  myplot[[1]] <-ggplot(data=cellnum[cellnum$celltype==unique(cellnum$celltype)[1],],
                       aes(x=group,y=Freq,fill=group))+
    geom_violin()+
    geom_boxplot(width=0.2,
                 position = position_dodge(0.9))+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1))+
    scale_fill_manual(values = pal_npg("nrc")(10))+facet_wrap(~celltype)+
stat_compare_means(aes(group=group),label = "p.signif",method = "wilcox.test",hide.ns = T)+
geom_signif(comparisons = compaired,step_increase = 0.5,map_signif_level = F,test = 't.test')#t.test, wilcox.test
}
}

#####################################################################
##################################单细胞数目条形图
color_cluster=c("#e5192c","#3a77b7","#3cac4c","#813c93","#f36c24",
                "#37b8c3","#a54922","#6b7627","#28996b",
                "#965b6a","#e9148f","#595b5e","#76c3ad",
                "#80d08a","#d29099","#f2e010")
#将每个细胞的identity转换成celltype.group
DefaultAssay(obj.combined) <- "RNA"
Idents(obj.combined) <-"celltype" ##保存每个细胞的细胞类型
obj.combined$celltype.group <- paste(Idents(obj.combined), obj.combined$group, sep = "_")   ##细胞类型前面加上分组信息

sce <- obj.combined
bar.df=sce@meta.data
text.df=as.data.frame(table(bar.df$celltype.group))

p1=bar.df%>%ggplot(aes(x=celltype.group))+geom_bar(aes(fill=celltype.group))+
  scale_x_discrete("")+
  scale_y_continuous("cell number",expand = c(0.02,0))+
  scale_fill_manual(values = color_cluster)+
  geom_text(data = text.df,aes(x=Var1,y=Freq,label=Freq),size=8)+
  theme_bw()+
  theme(#panel.grid = element_blank(),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1))+
theme(axis.text = element_text(size = 20),
      axis.title = element_text(size =20),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20),
      title = element_text(size = 20))+
  guides(color = guide_legend(override.aes = list(size = 4)))



pdf("barplot_3group.pdf",width = 12,height = 8)
p1
dev.off()
