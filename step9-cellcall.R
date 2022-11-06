setwd("/home/data/t040243/CopyOfjiangshijiu")
dir.create("cellcall")
setwd("cellcall")
library(devtools)
#devtools::install_github("ShellyCoder/cellcall")
library("cellcall")
if(!require(networkD3)){
  BiocManager::install("networkD3")
}

#devtools::install_github("wch/webshot")
library(webshot)
if(!is_phantomjs_installed()){
  install_phantomjs()
}
is_phantomjs_installed()
## [1] TRUE

#f.tmp <- system.file("extdata", "example_Data.Rdata", package="cellcall")
#load(f.tmp)
## gene expression stored in the variable in.content
#dim(in.content)
#in.content[1:4, 1:4]
#table(str_split(colnames(in.content), "_", simplify = T)[,2])
#mt <- CreateNichConObject(data=in.content, Shamn.feature = 3,
#                          names.field = 2,
#                          names.delim = "_",
#                          source = "TPM",
#                          scale.factor = 10^6,
#                          Org = "Mus musculus",
#                          project = "Shamcroenvironment")

set.seed(123)
dir.create("Sham")
setwd("Sham")
load("/home/data/t040243/CopyOfjiangshijiu/Doublet/after_doublet_obj.combined.RData")
Idents(obj.combined) <- "group"
obj.combine <- subset(obj.combined,idents = c("Sham"))
Idents(obj.combine) <- "celltype"
obj.combine@active.assay <- "RNA"


mt <- CreateObject_fromSeurat(Seurat.object=obj.combine, 
                                slot="counts", 
                                cell_type="celltype",
                                data_source="UMI",
                                scale.factor = 10^6, 
                                Org = "Mus musculus")

mt <- TransCommuProfile(object = mt,
                        pValueCor = 0.05,
                        CorValue = 0.1,
                        topTargetCor=1,
                        p.adjust = 0.05,
                        use.type="median",
                        probs = 0.9,
                        method="weighted",
                        IS_core = TRUE,
                        Org = 'Mus musculus')


n <- mt@data$expr_l_r_log2_scale
n[1:3,1:3]


pathway.hyper.list <- lapply(colnames(n), function(i){
  print(i)
  tmp <- getHyperPathway(data = n, object = mt, cella_cellb = i, Org="Mus musculus")
  return(tmp)
})


myPub.df <- getForBubble(pathway.hyper.list, cella_cellb=colnames(n))
p <- plotBubble(myPub.df)+theme_bw()+
  theme(axis.text.x  = element_text(size = 30,angle = 270),
        axis.text.y  = element_text(size = 30),
        axis.title = element_text(size =30),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30),
        title = element_text(size = 30),
        strip.text = element_text(face = "bold", color = "#7570B3",
                                  hjust = 0, size = 50),
        strip.background = element_rect(fill = "#E6AB02", linetype = "dotted"))+
  guides(color = guide_legend(override.aes = list(size = 12)))
p
ggsave("Sham_pathway_Bubble.pdf",p,width = 25,height = 10,dpi = 300)
write.csv(myPub.df,file = "Sham_pathway_Bubble.csv")


#先保存一下
save(mt,obj.combine,pathway.hyper.list,file = "mt.RData")

#彩虹色 
#cols <- rainbow(9)  ##选取彩虹前9个颜色
#cell_color <- data.frame(color=c(rainbow(5)), stringsAsFactors = FALSE)
library(RColorBrewer)   #配色
cell_color <- data.frame(color=c(brewer.pal(5, "Dark2")), stringsAsFactors = FALSE)

rownames(cell_color) <- c("ILC2a","ILC2b","ILC1","ILCdc",  "ILCt")

ViewInterCircos(object = mt, font = 2, cellColor = cell_color, 
                lrColor = c("#F16B6F", "#84B1ED"),
                arr.type = "big.arrow",arr.length = 0.04,
                trackhight1 = 0.05, slot="expr_l_r_log2_scale",
                linkcolor.from.sender = TRUE,
                linkcolor = NULL, gap.degree = 2,
                order.vector=c("ILC2a","ILC2b","ILC1","ILCdc",  "ILCt"),
                trackhight2 = 0.032, track.margin2 = c(0.01,0.12), DIY = FALSE)


pdf("Sham_circle_L_R.pdf",width = 8,height = 6)
 ViewInterCircos(object = mt@data$expr_l_r_log2_scale, font =2, 
                cellColor = cell_color,
                lrColor = c("#F16B6F", "#84B1ED"),
                arr.type = "big.arrow",arr.length = 0.04,
                trackhight1 = 0.05, slot="expr_l_r_log2_scale",
                linkcolor.from.sender = TRUE,
                linkcolor = NULL, gap.degree = 2,
                order.vector=c("ILC2a","ILC2b","ILC1","ILCdc",  "ILCt"),
                trackhight2 = 0.032, track.margin2 = c(0.01,0.12), DIY = T)
dev.off()


expr_l_r_log2_scale <- mt@data[["expr_l_r_log2_scale"]]
write.csv(expr_l_r_log2_scale,file = "expr_l_r_log2_scale.csv")
pdf("Sham_Pheatmap_L_R.pdf",width = 24,height = 18)
 viewPheatmap(object = mt, slot="expr_l_r_log2_scale", show_rownames = T,
             show_colnames = T,treeheight_row=0, treeheight_col=50,
             cluster_rows = T,cluster_cols = F,fontsize = 36,angle_col = "315",  
             main="score",color = NULL)

dev.off()


#########################################################################################
mt <- LR2TF(object = mt, sender_cell="ILC2a", recevier_cell="ILCdc",
            slot="expr_l_r_log2_scale", org="Mus musculus")
head(mt@reductions$sankey)


sank <- LRT.Dimplot(mt, 
                    fontSize = 30, 
                    nodeWidth = 30, 
                    height = NULL, 
                    width = 1200,    
                    sinksRight=FALSE, 
                    DIY.color = FALSE)
networkD3::saveNetwork(sank, 
                       "ILC2a_ILCdc_full.html")

#webshot("sankey.html" , "sankey.png")
webshot("ILC2a_ILCdc_full.html" , "ILC2a_ILCdc_full.pdf")


#########################################################################
#循环批量出图
dir.create("sank")
setwd("sank")
celltype <- unique(mt@meta.data[["celltype"]])

for (j in 1:length(celltype)) {
for (i in 1:length(celltype)) {
    
  b<- tryCatch(
    LR2TF(object = mt, sender_cell=celltype[j], recevier_cell=celltype[i],
          slot="expr_l_r_log2_scale", org="Mus musculus") ,
    # warning=function(w){message('Waring @ ',x) ; return(NA)},
    error=function(e){message('error @ ') ;  return(NA)},
    finally = {message('next...')}
  )
  
  sank <- tryCatch(LRT.Dimplot(b, 
                      fontSize = 15, 
                   nodeWidth = 30, 
                   height = NULL, 
                   width = 1200,    
                   sinksRight=FALSE, 
                  DIY.color = FALSE) ,
                  # warning=function(w){message('Waring @ ',x) ; return(NA)},
                  error=function(e){message('error @ ') ;  return(NA)},
                  finally = {message('next...')}
  )
  
  tryCatch(networkD3::saveNetwork(sank, 
                                  paste0(celltype[j],"_",celltype[i],".html"))  ,
           # warning=function(w){message('Waring @ ',x) ; return(NA)},
           error=function(e){message('error @ ') ;  return(NA)},
           finally = {message('next...')}
  )
   tryCatch(webshot(paste0(celltype[j],"_",celltype[i],".html") , paste0(celltype[j],"_",celltype[i],".pdf")),
            # warning=function(w){message('Waring @ ',x) ; return(NA)},
            error=function(e){message('error @ ') ;  return(NA)},
            finally = {message('next...')}
   )
  
   }
}

##########################################################################################

###转录因子富集图，有时候，我们会关注到TF下游富集到的TGs。
#此工具提供显示或导出TG列表详细信息的选项，该列表存储在NichConObject@datacell_type@geneSets.
dir.create("../GSEA")
setwd("../GSEA")

w <- mt@data$gsea.list$ILCt@geneSets
ILCt.tf <- names(mt@data$gsea.list$ILCt@geneSets)
ILCt.tf

#########################################################批量出图
pdf("ILCt_GSEA.pdf",width = 12,height = 9)
getGSEAplot(gsea.list=mt@data$gsea.list, geneSetID=c(ILCt.tf), 
            myCelltype="ILCt", fc.list=mt@data$fc.list,  
            selectedGeneID = mt@data$gsea.list$ILCt@geneSets$"Nfkb1",    ###标出TF调控的基因
            mycol = NULL)
dev.off()
#######################################################
## gsea object
egmt <- mt@data$gsea.list$ILCt
## filter TF
egmt.df <- data.frame(egmt)
head(egmt.df[,1:6])
flag.index <- which(egmt.df$p.adjust < 0.05)

plot <- ridgeplot.DIY(x=egmt, fill="p.adjust", showCategory=flag.index, core_enrichment = T,
                      orderBy = "NES", decreasing = T)
plot
ggsave(paste0("ILCt_ridgeplot",".pdf"),plot,width = 12,height = 9,dpi = 300)


################################################################################
####桑基图调整与配色
library(magrittr)
library(dplyr)
tmp <- mt@reductions$sankey
tmp1 <- dplyr::filter(tmp, weight1>0) ## filter triple relation with weight1 (LR score)
tmp.df <- trans2tripleScore(tmp1)  ## transform weight1 and weight2 to one value (weight)
head(tmp.df)

## set the color of node in sankey graph
mycol.vector = c('#5d62b5','#29c3be','#f2726f','#62b58f','#bc95df', '#67cdf2', '#ffc533', '#5d62b5', '#29c3be')  
elments.num <-  tmp.df %>% unlist %>% unique %>% length()
mycol.vector.list <- rep(mycol.vector, times=ceiling(elments.num/length(mycol.vector)))

sankey_graph(df = tmp.df, axes=1:3, mycol = mycol.vector.list[1:elments.num], nudge_x = NULL, font.size = 4, boder.col="white", isGrandSon = F)



library(magrittr)
library(dplyr)
tmp <- mt@reductions$sankey
tmp1 <- dplyr::filter(tmp, weight1>0)  ## filter triple relation with weight1 (LR score)
tmp.df <- trans2tripleScore(tmp1)  ## transform weight1 and weight2 to one value (weight)

## set the color of node in sankey graph
mycol.vector = c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')
elments.num <-  length(unique(tmp.df$Ligand))
mycol.vector.list <- rep(mycol.vector, times=ceiling(elments.num/length(mycol.vector)))

sankey_graph(df = tmp.df, 
             axes=1:3,
             mycol = mycol.vector.list[1:elments.num], 
             isGrandSon = TRUE, 
             nudge_x = nudge_x, 
             font.size = 2, 
             boder.col="white",    
             set_alpha = 0.8)

