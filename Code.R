##"Integration of single-cell transcriptomics and bulk transcriptomics to explore prognostic and immunotherapeutic characteristics of nucleotide metabolism in lung adenocarcinoma"Main R codes analyzed
##Author：Zhang et al. 
##Date:2024-10-6

#Fig 2A-Single-cell analysis
library(Seurat)
folders=list.files('./',pattern='[0123456789t]$')
folders
scList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder,
                     min.cells = 3, min.features = 200)
})
BM <- merge(scList[[1]], 
            y = c(scList[[2]],scList[[3]]), 
            add.cell.ids = c("sample1", "sample2","sample3"), 
            project = "BM") 

GM <- merge(scList[[4]], 
            y = c(scList[[5]],scList[[6]]), 
            add.cell.ids = c("sample4", "sample5","sample6"), 
            project = "GM")
GM[["percent.mt"]] <- PercentageFeatureSet(GM,pattern = "^MT-")
BM[["percent.mt"]] <- PercentageFeatureSet(BM,pattern = "^MT-")

preQC_GM <- VlnPlot(GM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                    ncol = 3, 
                    group.by = "orig.ident", 
                    pt.size = 0)
preQC_BM <- VlnPlot(BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                    ncol = 3, 
                    group.by = "orig.ident", 
                    pt.size = 0)
GM <- subset(GM, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
BM <- subset(BM, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
postQC_GM <- VlnPlot(GM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                     ncol = 3, 
                     group.by = "orig.ident", 
                     pt.size = 0)
postQC_BM <- VlnPlot(BM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                     ncol = 3, 
                     group.by = "orig.ident", 
                     pt.size = 0)
library(Seurat)
library(dplyr)
sce <- merge(BM, GM)
sce <- NormalizeData(sce, normalization.method = "LogNormalize") 
sce <- FindVariableFeatures(sce , nfeatures = 2000)
top10 <- head(VariableFeatures(sce), 10)
plot1 <- VariableFeaturePlot(sce)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = FALSE)
plot2
sce <- ScaleData(sce, vars.to.regress = c("S.Score", "G2M.Score","percent.mt", "nCount_RNA"), verbose = T)
sce <- FindVariableFeatures(sce)
sce <- RunPCA(sce,npcs = 50, verbose = FALSE)
library(harmony)
sce <- RunTSNE(sce,  dims = 1:25, reduction = "harmony")
DimPlot(sce,reduction = "tsne",label=T,split.by  = 'orig.ident') 
sce <- FindNeighbors(sce, reduction = "harmony",dims = 1:25) 
sce  <- FindClusters(object = sce , resolution = 1, verbose = FALSE) 
DimPlot(sce, label = T)
Macrodata <- sce
new.cluster.ids <- c("0"="T_cell", 
                     "1"="Epithelial_cells", 
                     "2"="Mono/Macro", 
                     "3"="T_cell", 
                     "4"="Mono/Macro", 
                     "5"="Mono/Macro", 
                     "6"="T_cell", 
                     "7"="B_cell", 
                     "8"="T_cell", 
                     "9"="Mono/Macro", 
                     "10"="B_cell", 
                     "11"="T_cell", 
                     "12"="Mono/Macro",
                     "13"="B_cell",
                     "14"="Mono/Macro",
                     "15"="Mono/Macro",
                     "16"="Epithelial_cells",
                     "17"="Fibroblasts",
                     "18"="Epithelial_cells",
                     "19"="B_cell",
                     "20"="Endothelial_cells",
                     "21"="Epithelial_cells",
                     "22"="T_cell",
                     "23"="Mono/Macro",
                     "24"="Epithelial_cells",
                     "25"="Epithelial_cells")
Macrodata <- RenameIdents(Macrodata, new.cluster.ids)                        
levels(Macrodata) <- c("Epithelial_cells",
                       "Endothelial_cells",
                       "Fibroblasts",
                       "T_cell",
                       "B_cell",
                       "Mono/Macro")
Macrodata$celltype <- Macrodata@active.ident
DimPlot(Macrodata, group.by = "celltype")
save(Macrodata, file = "Macrodata.RData")
#Fig 2B-single-cell thermogram
library(ClusterGVis)
library(org.Hs.eg.db)
load(Macrodata.RData)
pbmc=Macrodata
# find markers for every cluster compared to all remaining cells
# report only the positive ones
pbmc.markers.all <- Seurat::FindAllMarkers(pbmc,
                                           only.pos = TRUE,
                                           min.pct = 0.25,
                                           logfc.threshold = 0.5)
#get top 10 genes
pbmc.markers <- pbmc.markers.all %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 20, wt = avg_log2FC)
# check
head(pbmc.markers)
# prepare data from seurat object
st.data <- prepareDataFromscRNA(object = pbmc,
                                diffData = pbmc.markers,
                                showAverage = TRUE)
# check
str(st.data)
# enrich for clusters
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 5,
                        seed = 5201314)
# check
head(enrich)
# add gene name
markGenes = unique(pbmc.markers$gene)[sample(1:length(unique(pbmc.markers$gene)),40,
                                             replace = F)]
# line plot
visCluster(object = st.data,
           plot.type = "line")
# heatmap plot
pdf('sc1.pdf',height = 10,width = 6,onefile = F)
visCluster(object = st.data,
           plot.type = "heatmap",
           column_names_rot = 45,
           markGenes = markGenes,
           cluster.order = c(1:9))
dev.off()
# add bar plot
pdf('sc2.pdf',height = 10,width = 14,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           cluster.order = c(1:6),
           go.col = rep(jjAnno::useMyCol("stallion",n = 6),each = 5),
           add.bar = T)
dev.off()
#Fig 2C-single-celled volcano diagram (geology)
load("Macrodata.RData")
library(scRNAtoolVis)
library(Seurat)
pbmc.markers <- FindAllMarkers(Macrodata, only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = 0)
library(scRNAtoolVis)

head(pbmc.markers,3)

jjVolcano(diffData = pbmc.markers)

jjVolcano(diffData = pbmc.markers,
          aesCol = c('purple','orange'))

jjVolcano(diffData = pbmc.markers,
          tile.col = corrplot::COL2('RdBu', 15)[4:12])

p=jjVolcano(diffData = pbmc.markers,
            tile.col = corrplot::COL2('RdBu', 15)[4:12],
            size  = 3.5,
            fontface = 'bold')
p
ggsave(p,file="plot1.pdf",width = 10,height = 5)
#Fig 2D-AddModuleScore score
library(Seurat)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(msigdbr)
load("Macrodata.RData")
sce.all<-Macrodata
#signature gene
paper<-"signature gene"
papermarker<-str_to_upper(trimws(strsplit(paper,',')[[1]]))
papermarker
#list 
features <- list(papermarker)
#AddModuleScore
sce_score<- AddModuleScore(sce.all,
                           features = features,
                           ctrl = 100,
                           name = "features")
head(sce_score@meta.data)
colnames(sce_score@meta.data)[8] <- 'Score'
sce_score@meta.data[1:2,1:8]
library(ggpubr)
ggviolin(sce_score@meta.data, 
         x="celltype", y="Score", 
         width = 0.8,color = "black",
         fill="celltype",
         xlab = F, 
         add = 'mean_sd', 
         bxp.errorbar=T,      
         bxp.errorbar.width=0.05,          
         size=0.5, 
         palette = "npg", 
         legend = "right")
ggsave(filename="violin.pdf",width = 14,height = 8)
#Fig 2E-G-Copykat analysis and Aucell scores
library(tidyverse)
library(Seurat)
library(copykat)
library(ggplot2)
EP <- load("data.RData")
counts <- as.matrix(EP@assays$RNA@counts)
cnv <- copykat(rawmat=counts, ngene.chr=5, sam.name="test", n.cores=1)
save(EP,file = "RE.RData")
saveRDS(cnv, "cnv_1107.rds") 
table(rownames(cnv$CNAmat))
sc_cnv <-readRDS('cnv_1107.rds') 
single_RNA<- EP
all_cells <- sc_cnv$prediction[sc_cnv$prediction$copykat.pred != 'not.defined',]
sc_counts <- as.matrix(single_RNA@assays$RNA@counts)[,all_cells$cell.names]
CRC.data <- Matrix::Matrix(sc_counts)
single_RNA<- CreateSeuratObject(counts = CRC.data, project = "CRC_Eds", min.cells = 3, min.features = 200)
single_RNA[['group_copykat']] <- all_cells[match(rownames(single_RNA@meta.data),all_cells$cell.names),2]
single_RNA<- NormalizeData(single_RNA, normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(single_RNA)) 
single_RNA<- RunPCA(single_RNA, features = VariableFeatures(object = single_RNA)) %>%
  RuntSNE(dims = 1:20,reduction = 'pca',) %>%
  RunTSNE(dims = 1:20,reduction = 'pca') %>%
  FindNeighbors(dims = 1:20,reduction = 'pca') %>%
  FindClusters(resolution = 0.8)
plot1 <- DimPlot(single_RNA, reduction = "tsne",label = TRUE,label.box = TRUE,repel = TRUE,pt.size = 1)
plot2 <- DimPlot(single_RNA, reduction = "tsne",group.by = 'group_copykat',label = TRUE,label.box = TRUE,pt.size = 1,cols=c("#e3716e",'#54beaa'))
ggsave(plot2,file="1.pdf",width = 10,height = 8)
plot1 + plot2
plot3 <- DimPlot(single_RNA, reduction = "tSNE",label = TRUE,label.box = TRUE,repel = TRUE,pt.size = 1)
plot4 <- DimPlot(single_RNA, reduction = "tSNE",group.by = 'group_copykat',label = TRUE,label.box = TRUE,pt.size = 1)
plot3 + plot4
save(single_RNA,file = "single_RNA.RData")
EPI=as.data.frame(single_RNA$group_copykat)
write.csv(plotdata,file = "EP.csv")
library(tidyverse)
library(ggplot2)
library(cowplot)
library(AUCell) 
library(clusterProfiler)
Geneset <- read.gmt("gene.gmt") 
Geneset2 <- lapply(unique(Geneset$term), function(x){print(x);Geneset$gene[Geneset$term == x]})
names(Geneset2) <- unique(Geneset$term)
# step1
cells_rankings <- AUCell_buildRankings(scedata@assays$RNA@data)
save(cells_rankings,file = "AUCellrank.rda")
# step2
cells_AUC <- AUCell_calcAUC(Geneset2, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
set.seed(123)
par(mfrow=c(1,1)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 
geneSetID <- "gene"
aucs <- as.numeric(getAUC(cells_AUC)[geneSetID, ])
scedata$AUC  <- aucs
plotdata <- data.frame(scedata@meta.data, scedata@reductions$tsne@cell.embeddings)
library(tidyverse)
median_df <- plotdata %>% group_by(group_copykat) %>%
  summarise(median.1 = median(tSNE_1),
            median.2 = median(tSNE_2)) 
head(median_df) 
ggplot(plotdata,aes(tSNE_1,tSNE_2,col= AUC))+
  geom_point(size=0.1)+
  scale_colour_gradientn(colours = c("#92afe3","#1E90FF","orange","#E03131","#C92A2A"))+
  ggrepel::geom_text_repel(data = median_df,
                           aes(median.1,median.2, label =group_copykat),
                           size=4,
                           color= "gray20",
                           min.segment.length = 0)+
  theme_cowplot()+
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "right",
        axis.title = element_text(size = 12),
        axis.text = element_text(colour = 'black',size = 10))
ggsave("AUCell_tSNE.pdf",height = 8,width = 10)

#Fig 3A-B、D-G-Model construction and KM curve, timeROX curve
library(ggpubr)
module.gene=read.delim('gene.txt',sep='\t',header = F)
module.gene=module.gene$V1
tcga_dat_m<-read.delim('500tpms01A_log2.txt',sep='\t',header = T,row.names = 1,check.names = F)
tcga_cli<-read.delim('clinical01A.txt',sep='\t',header = T)
tcga_dat_m=t(tcga_dat_m[module.gene,])
tcga_dat_m=cbind.data.frame(OS.time=tcga_cli$OS.time,
                            OS=tcga_cli$OS,
                            tcga_dat_m[tcga_cli$Samples,])
library(survival)
tcga.cox=cox_batch(dat = t(tcga_dat_m[,-c(1,2)]),
                   time = tcga_dat_m$OS.time,
                   event = tcga_dat_m$OS)
head(tcga.cox)
fitp=0.05
tcga.sig.fit=tcga.cox[which(tcga.cox$p.value<fitp),]
dim(tcga.sig.fit)
sig_gene<-rownames(tcga.sig.fit)
length(sig_gene) 
write.csv(tcga.sig.fit,'tcga.sig.fit.csv',quote = F,row.names = T)
options(ggrepel.max.overlaps=Inf)
tcga.lasso<-lasso_cox(dat=tcga_dat_m[,sig_gene],
                      time=tcga_dat_m$OS.time,
                      event=tcga_dat_m$OS,
                      labels=c('A','B'))
length(tcga.lasso$gene)
tcga.lasso$plot
ggsave(plot = tcga.lasso$plot,filename = 'tcga.lasso.pdf',height = 5,width = 10)
tcga.lasso$lambda
tcga_risk<-get_riskscore(dat=tcga_dat_m[,-c(1,2)],
                         time=tcga_dat_m$OS.time,
                         event=tcga_dat_m$OS,
                         #cv.gene=rownames(tcga.sig.fit),
                         cv.gene=tcga.lasso$gene,
                         step=T,
                         cutoff='median')
length(tcga_risk$module_gene)#9
tcga_risk$module_gene
tcga_risk$formula
tcga_roc_km<-risk_plot(time=tcga_risk$result$time/365,
                       event=tcga_risk$result$status,
                       riskscore=tcga_risk$result$riskscore,
                       group=tcga_risk$result$risktype,
                       mk=c(1,2,3,4,5),labs=c('High','Low'),
                       palette=ggsci::pal_aaas()(10)[c(2,1)])
tcga_roc_km
ggsave(plot = tcga_roc_km,filename = 'tcga.roc.km.pdf',height = 5,width = 10)
#GSE31210
gse31210_dat_m<-read.delim('GSE31210.txt',sep='\t',header = T,row.names = 1,check.names = F)
gse31210_dat_m<-t(gse31210_dat_m[intersect(rownames(gse31210_dat_m),module.gene),])
range(gse31210_dat_m)
gse31210_cli<-read.delim('Cli.txt',sep='\t',header = T)
#GSE72094
gse72094_dat_m<-read.delim('GSE72094.txt',sep='\t',header = T,row.names = 1,check.names = F)
gse72094_dat_m<-t(gse72094_dat_m[intersect(rownames(gse72094_dat_m),module.gene),])
gse72094_cli<-read.delim('Cli.txt',sep='\t',header = T)
#GSE50081
gse50081_dat_m<-read.delim('GSE50081.txt',sep='\t',header = T,row.names = 1,check.names = F)
gse50081_dat_m<-t(gse50081_dat_m[intersect(rownames(gse50081_dat_m),module.gene),])
gse50081_cli<-read.delim('Cli.txt',sep='\t',header = T)
gse31210_dat_m<-cbind.data.frame(OS.time=gse31210_cli$OS.time,
                                 OS=gse31210_cli$OS,
                                 gse31210_dat_m[gse31210_cli$Accession,])
gse72094_dat_m<-cbind.data.frame(OS.time=gse72094_cli$OS.time,
                                 OS=gse72094_cli$OS,
                                 gse72094_dat_m[gse72094_cli$Accession,])

gse50081_dat_m<-cbind.data.frame(OS.time=gse50081_cli$OS.time,
                                 OS=gse50081_cli$OS,
                                 gse50081_dat_m[gse50081_cli$Accession,])
#GEO dataset validation
gse31210_risk<-get_riskscore(dat=gse31210_dat_m[,-c(1,2)],
                             time=gse31210_dat_m$OS.time,
                             event=gse31210_dat_m$OS,
                             cv.gene=intersect(tcga_risk$module_gene,colnames(gse31210_dat_m)),
                             step=F,cutoff='median')

gse31210_roc_km<-risk_plot(time=gse31210_risk$result$time/365,
                           event=gse31210_risk$result$status,
                           riskscore=gse31210_risk$result$riskscore,
                           group=gse31210_risk$result$risktype,
                           mk=c(1,2,3,4,5),labs=c('High','Low'),
                           palette=ggsci::pal_aaas()(10)[c(2,1)])
gse31210_roc_km
ggsave(plot = gse31210_roc_km,filename = 'gse31210_roc_km.pdf',height = 5,width = 10)
gse72094_risk<-get_riskscore(dat=gse72094_dat_m[,-c(1,2)],
                             time=gse72094_dat_m$OS.time,
                             event=gse72094_dat_m$OS,
                             cv.gene=intersect(tcga_risk$module_gene,colnames(gse72094_dat_m)),
                             step=F,cutoff='median')
gse72094_roc_km<-risk_plot(time=gse72094_risk$result$time/365,
                           event=gse72094_risk$result$status,
                           riskscore=gse72094_risk$result$riskscore,
                           group=gse72094_risk$result$risktype,
                           mk=c(1,2,3,4,5),labs=c('High','Low'),#mk
                           palette=ggsci::pal_aaas()(10)[c(2,1)])
gse72094_roc_km
ggsave(plot = gse72094_roc_km,filename = 'gse72094_roc_km.pdf',height = 5,width = 10)
#Fig 3C
gene.coef=read.delim('coef.txt',sep='\t',header = T)
gene.coef$Type=factor(gene.coef$Type,levels=c('Risk','Protective'))
library(dplyr)
library(ggplot2)
fig1d=gene.coef %>% 
  ggplot(aes(reorder(Gene, Coef), Coef)) +
  geom_col(aes(fill = Type)) +
  coord_flip() +
  scale_fill_manual(values=ggsci::pal_lancet('lanonc',alpha =0.7)(8)[c(7,3)]) +
  coord_flip() +
  labs(x = "") +
  labs(y = "Cox coefficient") +
  theme_bw()+
  theme(axis.text.y = element_text(angle = 0, hjust = 1),legend.position="top")
fig1d
ggsave('fig.pdf',fig1d,height = 6,width = 6)

#Figure 4A-D-model comparison
suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(randomForestSRC)
  library(glmnet)
  library(plsRcox)
  library(superpc)
  library(gbm)
  library(CoxBoost)
  library(survivalsvm)
  library(tibble)
  library(BART)
  library(stringr)
  library(ComplexHeatmap)
})
tcga.data <- read.table("500tpms01A_log2.txt",sep = "\t",row.names = 1,header = T,stringsAsFactors = F,check.names = F)
tcga.cli <- read.table("cli.txt")
meta.data = read.table("GEO.txt",sep = "\t",row.names = 1,header = T,stringsAsFactors = F,check.names = F)
meta.cli = read.table("GEO-cli.txt")
signatures = data.table::fread('signatures.txt')
sig2df = function(sig){
  rs = list()
  sig1 = sig %>% filter(Gene %in% rownames(tcga.data))
  RS1 = colSums(tcga.data[sig1$Gene,] * sig1$Coef) %>% data.frame(RS=.)
  rs$TCGA <- tcga.cli %>% select(OS.time,OS) %>% cbind(RS1)
  
  sig2 = sig %>% filter(Gene %in% rownames(meta.data))
  RS2 = colSums(meta.data[sig2$Gene,] * sig2$Coef) %>% data.frame(RS=.)
  rs$META <- meta.cli %>% select(OS.time,OS) %>% cbind(RS2)
  rs2df <- function(rs,modelName){
    cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])})) %>% 
      rownames_to_column('ID')
    cc$Model  <- modelName
    return(cc)
  }
  cc = rs2df(rs,sig$Author[1])
  return(cc)
}
calcoef <- function(cli,data,module.gene)
module.name='Lasso+StepCox'
cindex.NMBRS <- read.table("cindex.txt")[module.name,] %>% 
  pivot_longer(names_to = "ID",values_to = "Cindex",cols = 1:2) %>% 
  mutate(Model="NMBRS")
cindex.22signatures = split(signatures,signatures$Author) %>% lapply(sig2df) %>% Reduce(rbind,.)
cindex = rbind(cindex.NMBRS,cindex.22signatures)
cindex.select = cindex %>% filter(ID=="TCGA") %>% arrange(Cindex)
cindex.select$Model <- factor(cindex.select$Model,levels = cindex.select$Model)
p <- cindex.select %>% mutate(current=Model=="NMBRS") %>%  ggplot() + theme_classic() +  
  geom_vline(xintercept = 0.6,linewidth=1.5,linetype = "dashed") +
  geom_segment(linewidth=1.5,linetype = "dashed",
               aes(x = 0, y = Model,xend = Cindex, yend = Model,color=current)
  ) + geom_point(aes(x=Cindex,y=Model,color=current),size=6) +
  labs(x="TCGA Cindex",y='') + scale_color_manual(values=c("#4F9ACA", "red")) +
  theme(legend.position = 'none')
ggsave("Fig1.pdf",p,height = 8,width = 4)
cindex.select = cindex %>% filter(ID=="GEO") %>% arrange(Cindex)
cindex.select$Model <- factor(cindex.select$Model,levels = cindex.select$Model)
p <- cindex.select %>% mutate(current=Model=="NMBRS") %>%  ggplot() + theme_classic() +  
  geom_vline(xintercept = 0.6,linewidth=1.5,linetype = "dashed") +
  geom_segment(linewidth=1.5,linetype = "dashed",
               aes(x = 0, y = Model,xend = Cindex, yend = Model,color=current)
  ) + geom_point(aes(x=Cindex,y=Model,color=current),size=6) +
  labs(x="GEO Cindex",y='') + scale_color_manual(values=c("#4F9ACA", "red")) +
  theme(legend.position = 'none')
ggsave("Fig2.pdf",p,height = 8,width = 4)
#Figure 4E-F-Clinicopathological characterization
sig_boxplot<-function(dat,leg,ylab,palette=ggsci::pal_aaas()(10)[3:4]){
  library(ggpubr)
  dat=na.omit(dat)
  colnames(dat)=c('group','gene')
  dat=dat[order(dat$group),]
  all.combn=combn(as.character(unique(dat$group)),2)
  my_comparisons=lapply(seq_len(ncol(all.combn)), function(i) all.combn[,i])
  pp=ggboxplot(dat, 
               x='group', y='gene', fill = 'group',
               palette =  palette,
               short.panel.labs = T,outlier.shape = NA)+
    stat_compare_means(comparisons=my_comparisons,method="wilcox.test",label = "p.signif")+
    ylab(ylab)+xlab('')+labs(fill=leg)
  return(pp)
}
tcga.subtype.cli<-read.delim('clinical01A.txt',sep='\t',header = T)
tcga.risk<-read.delim('risk.txt',sep='\t',header = T)
rownames(tcga.risk)=tcga.risk$Samples
rownames(tcga.subtype.cli)=tcga.subtype.cli$Samples

tcga.subtype.cli=cbind.data.frame(tcga.subtype.cli[tcga.risk$Samples,],
                                  riskscore=tcga.risk$riskscore,
                                  RiskType=tcga.risk$risktype)
head(tcga.subtype.cli)
write.table(tcga.subtype.cli,'tcga.risk.cli.txt',quote = F,row.names = F,sep='\t')
fig9a=list()
fig9a[[1]]=sig_boxplot(tcga.subtype.cli[,c("Gender","riskscore")],
                       leg = 'Gender',ylab = 'NMBRS',
                       palette =ggsci::pal_aaas()(10))
fig9a[[1]]
fig9a[[2]]=sig_boxplot(tcga.subtype.cli[,c("T.Stage","riskscore")],
                       leg = 'T.Stage',ylab = 'NMBRS',
                       palette =ggsci::pal_aaas()(10))
fig9a[[2]]
fig9a[[3]]=sig_boxplot(tcga.subtype.cli[,c("N.Stage","riskscore")],
                       leg = 'N.Stage',ylab = 'NMBRS',
                       palette =ggsci::pal_aaas()(10))
fig9a[[3]]
fig9a[[4]]=sig_boxplot(tcga.subtype.cli[,c("M.Stage","riskscore")],
                       leg = 'M.Stage',ylab = 'NMBRS',
                       palette =ggsci::pal_aaas()(10))
fig9a[[4]]
fig9a[[5]]=sig_boxplot(tcga.subtype.cli[,c("Stage","riskscore")],
                       leg = 'Stage',ylab = 'NMBRS',
                       palette =ggsci::pal_aaas()(10))
fig9a[[5]]
fig9a[[6]]=sig_boxplot(tcga.subtype.cli[,c("Age","riskscore")],
                       leg = 'Age',ylab = 'NMBRS',
                       palette =ggsci::pal_aaas()(10))
fig9a[[6]]
library(ggpubr)
ggsave('fig.pdf',fig9a,height = 9,width = 10)

#Supplementary Figure-Nomogram
tcga.risk<-read.delim('tcga.risk.txt',sep='\t',header = T)
head(tcga.risk)
tcga.cli<-read.delim('clinical01A.txt',sep='\t',header = T)
tcga.risk.cli<-merge(tcga.risk,tcga.cli,by='Samples')
write.table(tcga.risk.cli,'tcga.risk.cli.txt',quote = F,row.names = F,sep='\t')
library(forestplot)
library(survcomp)
tcga_cox_datas=tcga.risk.cli
table(tcga_cox_datas$T.Stage)
tcga_cox_datas$T.Stage[tcga_cox_datas$T.Stage == 'T1' | tcga_cox_datas$T.Stage == 'T2'] <- 'T1+T2'
tcga_cox_datas$T.Stage[tcga_cox_datas$T.Stage == 'T3' | tcga_cox_datas$T.Stage == 'T4'] <- 'T3+T4'

table(tcga_cox_datas$N.Stage)
tcga_cox_datas$N.Stage[tcga_cox_datas$N.Stage == 'N1' | tcga_cox_datas$N.Stage == 'N2' | tcga_cox_datas$N.Stage == 'N3'] <- 'N1+N2+N3'

table(tcga_cox_datas$M.Stage)

table(tcga_cox_datas$Stage)
tcga_cox_datas$Stage[tcga_cox_datas$Stage == 'I' | tcga_cox_datas$Stage == 'II'] <- 'I+II'
tcga_cox_datas$Stage[tcga_cox_datas$Stage == 'III' | tcga_cox_datas$Stage == 'IV'] <- 'III+IV'

table(tcga_cox_datas$Gender)

Age_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age,
                             data=tcga_cox_datas))
Age_sig_cox_dat <- data.frame(Names=rownames(Age_sig_cox[[8]]),
                              HR = round(Age_sig_cox[[7]][,2],3),
                              lower.95 = round(Age_sig_cox[[8]][,3],3),
                              upper.95 = round(Age_sig_cox[[8]][,4],3),
                              p.value=round(Age_sig_cox[[7]][,5],3))
Age_sig_cox_dat
T_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~T.Stage,
                           data=tcga_cox_datas))
T_sig_cox_dat <- data.frame(Names=rownames(T_sig_cox[[8]]),
                            HR = round(T_sig_cox[[7]][,2],3),
                            lower.95 = round(T_sig_cox[[8]][,3],3),
                            upper.95 = round(T_sig_cox[[8]][,4],3),
                            p.value=round(T_sig_cox[[7]][,5],3))
T_sig_cox_dat
N_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~N.Stage,
                           data=tcga_cox_datas))
N_sig_cox_dat <- data.frame(Names=rownames(N_sig_cox[[8]]),
                            HR = round(N_sig_cox[[7]][,2],3),
                            lower.95 = round(N_sig_cox[[8]][,3],3),
                            upper.95 = round(N_sig_cox[[8]][,4],3),
                            p.value=round(N_sig_cox[[7]][,5],3))
N_sig_cox_dat
M_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~M.Stage,
                           data=tcga_cox_datas))
M_sig_cox_dat <- data.frame(Names=rownames(M_sig_cox[[8]]),
                            HR = round(M_sig_cox[[7]][,2],3),
                            lower.95 = round(M_sig_cox[[8]][,3],3),
                            upper.95 = round(M_sig_cox[[8]][,4],3),
                            p.value=round(M_sig_cox[[7]][,5],3))
M_sig_cox_dat
Stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Stage,
                               data=tcga_cox_datas))
Stage_sig_cox_dat <- data.frame(Names=rownames(Stage_sig_cox[[8]]),
                                HR = round(Stage_sig_cox[[7]][,2],3),
                                lower.95 = round(Stage_sig_cox[[8]][,3],3),
                                upper.95 = round(Stage_sig_cox[[8]][,4],3),
                                p.value=round(Stage_sig_cox[[7]][,5],3))
Stage_sig_cox_dat
RiskType_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~riskscore,
                                  data=tcga_cox_datas))
RiskType_sig_cox_dat <- data.frame(Names=rownames(RiskType_sig_cox[[8]]),
                                   HR = round(RiskType_sig_cox[[7]][,2],3),
                                   lower.95 = round(RiskType_sig_cox[[8]][,3],3),
                                   upper.95 = round(RiskType_sig_cox[[8]][,4],3),
                                   p.value=round(RiskType_sig_cox[[7]][,5],3))
RiskType_sig_cox_dat
Gender_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Gender,
                                data=tcga_cox_datas))
Gender_sig_cox_dat <- data.frame(Names=rownames(Gender_sig_cox[[8]]),
                                 HR = round(Gender_sig_cox[[7]][,2],3),
                                 lower.95 = round(Gender_sig_cox[[8]][,3],3),
                                 upper.95 = round(Gender_sig_cox[[8]][,4],3),
                                 p.value=round(Gender_sig_cox[[7]][,5],3))
Gender_sig_cox_dat
sig_cox_dat <- rbind(Age_sig_cox_dat,
                     Gender_sig_cox_dat,
                     T_sig_cox_dat,
                     N_sig_cox_dat,
                     M_sig_cox_dat,
                     Stage_sig_cox_dat,
                     RiskType_sig_cox_dat)
data.sig <- data.frame(Names=sig_cox_dat$Names,
                       p.value=sig_cox_dat$p.value,
                       sig_cox_dat$HR,
                       sig_cox_dat$lower.95,
                       sig_cox_dat$upper.95)
rownames(data.sig) <- c("Age",
                        "Gender",
                        "T.Stage",
                        "N.Stage",
                        "M.Stage",
                        "Stage",
                        "RiskScore")
data.sig$Names <- rownames(data.sig)
data.sig$p.value=ifelse(data.sig$p.value==0.000,'<0.001',data.sig$p.value)
write.table(data.sig,'results/data.sig.txt',sep='\t',row.names = T,quote = F)
vio_fores<-function(dat){
  library(forestplot)
  col=fpColors(box='red',summary='#8B008B',lines = 'grey',zero = 'grey')
  smary=rep(F,length(as.numeric(dat[,ncol(dat)-2])))
  lower=as.numeric(dat[,ncol(dat)-1])
  upper=as.numeric(dat[,ncol(dat)])
  nind=which(is.na(lower)|is.na(upper)|is.na(as.numeric(dat[,ncol(dat)-2])))
  smary[nind]=T
  labeltext=as.matrix(dat[,1:(ncol(dat)-3)])
  adt=paste0(dat[,3],'(',dat[,4],',',dat[,5],')')
  
  labeltext=cbind(labeltext,adt)
  colnames(labeltext)[3]=c('Hazard Ratio(95% CI)')
  
  hz_list=list('2'=gpar(lty=1,col='#8B008B'),
               '3'=gpar(lty=1,col='#8B008B'))
  names(hz_list)=c(2,nrow(labeltext)+2)
  forestplot(labeltext = rbind(colnames(labeltext),labeltext),
             hrzl_lines = hz_list,
             mean = c(NA,as.numeric(dat[,ncol(dat)-2])),
             lower =c(NA,lower),
             upper = c(NA,upper),
             is.summary=c(T,smary),
             zero = 1,
             boxsize = 0.4, 
             lineheight = unit(10,'mm'),
             colgap = unit(8,'mm'),
             lwd.zero = 2,
             lwd.ci = 2,
             col=col,
             xlab='HR',
             lwd.xaxis=2,
             lty.ci = 'solid',
             clip = c(min(lower,na.rm = T),max(upper,na.rm = T)),
             xlog=T,
             graph.pos = 2)
  
  
}
#dat=data.sig
pdf('results/sig_cox.pdf', width = 8, height = 5)
vio_fores(data.sig)
dev.off()
#Results of multifactor analysis
muti_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age+Gender+T.Stage+N.Stage+M.Stage+Stage+riskscore, 
                              data=tcga_cox_datas))
muti_cox_dat <- data.frame(Names=rownames(muti_sig_cox[[8]]),
                           HR = round(muti_sig_cox[[7]][,2],3),
                           lower.95 = round(muti_sig_cox[[8]][,3],3),
                           upper.95 = round(muti_sig_cox[[8]][,4],3),
                           p.value=round(muti_sig_cox[[7]][,5],3))
data.muti <- data.frame(Names=muti_cox_dat$Names,
                        p.value=muti_cox_dat$p.value,
                        muti_cox_dat$HR,
                        muti_cox_dat$lower.95,
                        muti_cox_dat$upper.95)
rownames(data.muti) <- c("Age","Gender","T.Stage","N.Stage","M.Stage","Stage","RiskScore")
data.muti$Names <- rownames(data.muti)
data.muti$p.value=ifelse(data.muti$p.value==0.000,'<0.001',data.muti$p.value)
write.table(data.muti,'results/data.muti.txt',sep='\t',row.names = T,quote = F)
pdf('results/muti_cox.pdf', width = 8, height = 5)
vio_fores(data.muti)
dev.off()
#Line diagrams and DCA diagrams 
crbind2DataFrame=function(dat,full=F){
  print(class(dat))
  if(class(dat)[1]=='table'){
    if(!is.na(ncol(dat))){
      dat=apply(dat,2,function(x){
        return(x)
      })
    }
  }
  if(class(dat)[1]!='data.frame'){
    dat1=as.data.frame(as.matrix(dat))
  }else{
    dat1=dat
  }
  dat1.class=apply(dat1, 2, class)
  #which(dat1.class!='numeric')
  #print(head(dat1))
  for(i in which(dat1.class!='numeric')){
    dat1[,i]=as.character(dat1[,i])
    if(full){
      dat1[,i]=as.numeric(dat1[,i])
    }else{
      dt=dat1[which(gsub(' ','',dat1[,i])!=''&!is.na(dat1[,i])),i]
      dt=dt[which(dt!='Inf'&dt!='NaN'&dt!='NA')]
      #dt[which(is.na(as.numeric(dt)))]
      if(sum(is.na(as.numeric(dt)))<length(dt)*0.1){
        #print(dat1[,i])
        dat1[,i]=as.numeric(dat1[,i])
      }
    }
  }
  return(dat1)  
}
mg_plotDCA=function(status,fmlas,modelNames,data){
  set.seed(123)
  all.mod=list()
  for(i in 1:length(fmlas)){
    fmla <- as.formula(paste0("status~",fmlas[i]))
    model<-rmda::decision_curve(fmla,
                                data=data,
                                bootstraps=500)
    all.mod=c(all.mod,list(model))
  }
  rmda::plot_decision_curve(all.mod,
                            curve.names=modelNames,
                            xlim=c(0,1),legend.position="topright",
                            confidence.intervals=FALSE)
}
nomogram_plot=function(clinical_riskscore,
                       os,
                       status,
                       title='Nomogram',
                       quick=T,
                       mks = c(1,3,5)){
  mg_colors= ggsci::pal_npg("nrc")(10)
  norm.stat.al=data.frame(clinical_riskscore,time=os,status=status)
  norm.stat.al=as.data.frame(norm.stat.al)
  library(rms)
  env <- globalenv()
  env$MG_Grobal_DDSet <- rms::datadist(norm.stat.al) 
  options(datadist='MG_Grobal_DDSet')
  fmla <- as.formula(paste0("Surv(time, status) ~",paste0(colnames(clinical_riskscore),collapse = '+')))
  cox2 <- cph(fmla, data = norm.stat.al,surv = T,x = T,y = T)
  
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  #cindex.orig=1-rcorr.cens(fp,Surv(norm.stat.al$time,norm.stat.al$status))
  cut.time=c()
  if(quantile(os[!is.na(os)])['75%']<12){
    cut.time=mks
  }else if(quantile(os[!is.na(os)])['75%']<365){
    cut.time=c(12*mks[1],12*mks[2],12*mks[3])
  }else{
    cut.time=c(365*mks[1],365*mks[2],365*mks[3])
  }
  cut.time=cut.time[which(cut.time<quantile(os,seq(0,1,0.01))['95%'])]
  print(cut.time)
  surv=Survival(cox2)
  survs=list()
  cal_all=list()
  for(i in 1:length(cut.time)){
    f1<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[i]) 
    cal1<-calibrate(f1, cmethod="KM", method="boot",u=cut.time[i],m=floor(sum(f1$n)/3)) 
    cal_all=c(cal_all,list(cal1))
    #    surv0 <- function(x)surv(cut.time[i],lp=x) 
    #    survs=c(survs,list(surv0))
  }
  #f2<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[2]) 
  #cal3<-calibrate(f2, cmethod="KM", method="boot",u=cut.time[2],m=100,B=200) 
  #f3<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[3]) 
  #cal5<-calibrate(f3, cmethod="KM", method="boot",u=cut.time[3],m=100,B=200) 
  if(length(cut.time)==1){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    survs=list(surv1)
  }else if(length(cut.time)==2){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    survs=list(surv1,surv2)
  }else if(length(cut.time)==3){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    surv3 <- function(x)surv(cut.time[3],lp=x) 
    survs=list(surv1,surv2,surv3)
  }
  nom=nomogram(cox2,fun=survs,lp= F
               ,funlabel=c(paste0(mks[1], '-Year Survival'),
                           paste0(mks[2], '-Year Survival'),
                           paste0(mks[3], '-Year Survival'))[1:length(cut.time)]
               ,maxscale=100
               ,fun.at=seq(0,1,0.2)
  )
  
  if(!quick){
    cal_all=list()
    for(i in 1:length(cut.time)){
      cal1=get_best_calibrate(cox2,cut.time[i])
      cal_all=c(cal_all,list(cal1))
    }
    #cal3=get_best_calibrate(cox2,cut.time[2])
    #cal5=get_best_calibrate(cox2,cut.time[3])
  }
  lay2 <- customLayout::lay_new(matrix(1:2))
  lay1 <- customLayout::lay_new(matrix(1:1))
  cl <- customLayout::lay_bind_col(lay2, lay1, widths = c(1, 1.5)) 
  #customLayout::lay_show(cl)
  customLayout::lay_set(cl) 
  
  plot(cal_all[[1]],lwd = 2,lty = 1,errbar.col = mg_colors[1], bty = "l",xlim = c(0,1),ylim= c(0,1),xlab = "Nomogram-prediced OS (%)"
       ,ylab = "Observed OS (%)",col = mg_colors[1],cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6,subtitles = F)
  #lines(cal_all[[1]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[1], pch = 16)
  mtext("")
  if(length(cal_all)>1){
    for(i in 2:length(cal_all)){
      plot(cal_all[[i]],lwd = 2,lty = 1,errbar.col = mg_colors[i],xlim = c(0,1),ylim= c(0,1),col = mg_colors[i],add = T)
      #lines(cal_all[[i]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[i], pch = 16)
    }
    #plot(cal3,lwd = 2,lty = 0,errbar.col = mg_colors[3],xlim = c(0,1),ylim= c(0,1),col = mg_colors[3],add = T)
    #lines(cal3[,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[3], pch = 16)
  }
  abline(0,1, lwd = 2, lty = 3, col = 'black')
  legend("topleft", legend = c(paste0(mks[1], '-Year'),
                               paste0(mks[2], '-Year'),
                               paste0(mks[3], '-Year'))[1:length(cut.time)],col =mg_colors[1:length(cut.time)],lwd = 2,cex = 1.2,bty = "n")
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  dca_dat=data.frame(Nomogram=fp,status=norm.stat.al$status)
  fp.al=cbind(fp)
  for(i in 1:ncol(clinical_riskscore)){
    fmla1 <- as.formula(paste0("Surv(time, status) ~",colnames(clinical_riskscore)[i]))
    cox1 <- cph(fmla1, data = norm.stat.al,surv = T,x = T,y = T)
    fp1 <- predict(cox1)
    fp.al=cbind(fp.al,fp1)
  }
  colnames(fp.al)=c('Nomogram',colnames(clinical_riskscore))
  fp.al=crbind2DataFrame(fp.al)
  fp.al$status=norm.stat.al$status
  mg_plotDCA(fp.al$status
             ,c('Nomogram',colnames(clinical_riskscore))
             ,c('Nomogram',colnames(clinical_riskscore)),fp.al)
  #plot(cal1,xlim=c(0,1),ylim=c(0,1))
  #plot(cal3,xlim=c(0,1),ylim=c(0,1))
  #plot(cal5,xlim=c(0,1),ylim=c(0,1))
  plot(nom)
  options(datadist=NULL)
  return(list(Mod=cox2,Cindex=cindex,CutTime=cut.time))
}
nomogram_buti=function(cox_model,cut.time,title='Nomogram'){
  library(regplot)
  regplot(cox_model
          ,observation=pbc[2,] 
          ,title=title
          ,failtime = cut.time
          ,prfail = TRUE 
          ,showP = T 
          ,droplines = F
          ,points=TRUE)}
pdf('results/nomogram.pdf', width = 12, height = 10)
nom.plot=nomogram_plot(data.frame(N.Stage=tcga_cox_datas$N.Stage,
                                  T.Stage=tcga_cox_datas$T.Stage,
                                  RiskScore=tcga_cox_datas$riskscore,
                                  Gender=tcga_cox_datas$Gender,
                                  Age=tcga_cox_datas$Age1,
                                  Stage=tcga_cox_datas$Stage),
                       os = tcga_cox_datas$OS.time,
                       status = tcga_cox_datas$OS,
                       mks = c(1,3,5))
dev.off()
nomogram_buti(nom.plot$Mod,cut.time = c(1*365,3*365,5*365))
write.table(tcga_cox_datas,'results/tcga_cox_datas.txt',quote = F,row.names = F,sep='\t')
library("timeROC")
tcga_cox_auc=tcga_cox_datas
tcga_cox_auc$Riskscore=as.numeric(tcga_cox_auc$riskscore)
tcga_cox_auc$T.Stage=as.numeric(as.factor(tcga_cox_auc$T.Stage))
tcga_cox_auc$N.Stage=as.numeric(as.factor(tcga_cox_auc$N.Stage))
tcga_cox_auc$M.Stage=as.numeric(as.factor(tcga_cox_auc$M.Stage))
tcga_cox_auc$Stage=as.numeric(as.factor(tcga_cox_auc$Stage))
tcga_cox_auc$Gender=as.numeric(as.factor(tcga_cox_auc$Gender))
tcga_cox_auc$Age=as.numeric(as.factor(tcga_cox_auc$Age))
ROC.DSST.Age=timeROC(T=tcga_cox_auc$OS.time/365,
                     delta=tcga_cox_auc$OS,
                     marker=tcga_cox_auc$Age,
                     cause=1,weighting="marginal",
                     times=1:5,
                     iid=TRUE)
ROC.DSST.gender=timeROC(T=tcga_cox_auc$OS.time/365,
                        delta=tcga_cox_auc$OS,
                        marker=tcga_cox_auc$Gender,
                        cause=1,weighting="marginal",
                        times=1:5,
                        iid=TRUE)
ROC.DSST.T.Stage=timeROC(T=tcga_cox_auc$OS.time/365,
                         delta=tcga_cox_auc$OS,
                         marker=tcga_cox_auc$T.Stage,
                         cause=1,weighting="marginal",
                         times=1:5,
                         iid=TRUE)
ROC.DSST.N.Stage=timeROC(T=tcga_cox_auc$OS.time/365,
                         delta=tcga_cox_auc$OS,
                         marker=tcga_cox_auc$N.Stage,
                         cause=1,weighting="marginal",
                         times=1:5,
                         iid=TRUE)
ROC.DSST.M.Stage=timeROC(T=tcga_cox_auc$OS.time/365,
                         delta=tcga_cox_auc$OS,
                         marker=tcga_cox_auc$M.Stage,
                         cause=1,weighting="marginal",
                         times=1:5,
                         iid=TRUE)
ROC.DSST.Stage=timeROC(T=tcga_cox_auc$OS.time/365,
                       delta=tcga_cox_auc$OS,
                       marker=tcga_cox_auc$Stage,
                       cause=1,weighting="marginal",
                       times=1:5,
                       iid=TRUE)
ROC.DSST.Risk=timeROC(T=tcga_cox_auc$OS.time/365,
                      delta=tcga_cox_auc$OS,
                      marker=tcga_cox_auc$riskscore,
                      cause=1,weighting="marginal",
                      times=1:5,
                      iid=TRUE)
ROC.DSST.Nomo<-timeROC(T=tcga_cox_auc$OS.time/365,
                       delta=tcga_cox_auc$OS,
                       marker=tcga_cox_auc$Riskscore,
                       other_markers=as.matrix(tcga_cox_auc[,c("T.Stage","N.Stage")]),
                       cause=1,
                       weighting="cox",
                       times=1:5,
                       iid=F)
ROC.DSST.Age$AUC
ROC.DSST.gender$AUC
ROC.DSST.T.Stage$AUC
ROC.DSST.N.Stage$AUC
ROC.DSST.M.Stage$AUC
ROC.DSST.Stage$AUC
ROC.DSST.Risk$AUC
ROC.DSST.Nomo$AUC
mg_colors=ggsci::pal_lancet()(9)
pdf('results/AUC.pdf',height = 7,width = 7)
plotAUCcurve(ROC.DSST.Nomo,conf.int=F,col=mg_colors[1])
plotAUCcurve(ROC.DSST.Risk,conf.int=F,col=mg_colors[2],add=TRUE)
plotAUCcurve(ROC.DSST.Stage,conf.int=F,col=mg_colors[3],add=TRUE)
plotAUCcurve(ROC.DSST.Age,conf.int=F,col=mg_colors[4],add=TRUE)
plotAUCcurve(ROC.DSST.T.Stage,conf.int=F,col=mg_colors[5],add=TRUE)
plotAUCcurve(ROC.DSST.N.Stage,conf.int=F,col=mg_colors[6],add=TRUE)
plotAUCcurve(ROC.DSST.M.Stage,conf.int=F,col=mg_colors[7],add=TRUE)
plotAUCcurve(ROC.DSST.gender,conf.int=F,col=mg_colors[8],add=TRUE)
legend("topright",c("Nomogram","RiskScore",'Stage',"Age",'T.Stage','N.Stage',"M.Stage",'Gender')
       ,col=mg_colors[c(1:8)],lty=1,lwd=2)
dev.off()

#Figure 5A-B-Genomic mutation analysis
library(pacman)
library(tidyverse)
library(maftools)
library(dplyr)
tmp <- read.table("5dfa07b9-3d69-4974-9a29-336b47b0cdab.wxs.aliquot_ensemble_masked.maf.gz",
                  sep = "\t", 
                  header = T )
tmp[1:10,1:8]
dir.path <- "gdc_download_20230816_060435.496462"
all.maf <- list.files(path = dir.path, pattern = ".gz", 
                      full.names = T, recursive = T)
all.maf[1:3]
maf.list <- lapply(all.maf, data.table::fread, 
                   sep = "\t", 
                   header = T,
                   skip = 7 )
maf.merge <- do.call(rbind,maf.list)
write.table(maf.merge,"maf.merge.txt",quote = F,row.names = F,sep='\t')
dim(maf.merge)
maf1 <- read.maf(maf.merge)
maf1@data$Tumor_Sample_Barcode 
mut<-maf1@data
mut$Tumor_Sample_Barcode<-substring(mut$Tumor_Sample_Barcode,1,16)
mut$Tumor_Sample_Barcode[1]
all_surv_expr<-read.delim('Risktype.txt',sep='\t',header = T,row.names = 1,check.names = F)
mut$Tumor_Sample_Barcode<-substring(mut$Tumor_Sample_Barcode,1,16)
mut.A <- mut[(mut$Tumor_Sample_Barcode %in% rownames(all_surv_expr)[all_surv_expr$Risktype=="High"]),]
mut.B <- mut[(mut$Tumor_Sample_Barcode %in% rownames(all_surv_expr)[all_surv_expr$Risktype=="Low"]),]
maf.A <-read.maf(maf = mut.A,isTCGA = T)
maf.B <- read.maf(maf = mut.B,isTCGA = T)
col = RColorBrewer::brewer.pal(n = 10, name = 'Paired')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Ins','In_Frame_Ins', 'Splice_Site', 'In_Frame_Del','Nonstop_Mutation','Translation_Start_Site','Multi_Hit')
racecolors = RColorBrewer::brewer.pal(n = 4,name = 'Spectral')
names(racecolors) = c("ASIAN", "WHITE", "BLACK_OR_AFRICAN_AMERICAN",  "AMERICAN_INDIAN_OR_ALASKA_NATIVE")
library(RColorBrewer)
vc_cols <- brewer.pal(9,"Set1")
names(vc_cols) <- levels(maf.all@data$Variant_Classification)
head(vc_cols)
pdf("onco_plot1.pdf", height = 6, width = 9)
oncoplot(maf = maf.A, top = 20,colors = vc_cols)
dev.off()
library(RColorBrewer)
vc_cols <- brewer.pal(9,"Set1")
names(vc_cols) <- levels(maf.B@data$Variant_Classification)
head(vc_cols)
pdf("onco_plot2.pdf", height = 6, width = 9)
oncoplot(maf = maf.B, top = 20,colors = vc_cols)
dev.off()
#Figure 5C-D、F-G-Comparison of TMB, MATH
data=read.delim('DATA.txt',sep='\t',header = T,check.names = F)
library(ggplot2)
library(ggpubr)
library(reshape2)
ggviolin(data, x="Risktype", y="TMB", 
         color = "Risktype",
         fill="Risktype",
         palette =c("#F9BEBB","#89C9C8"),
         add = "boxplot",
         add.params = list(color="white"),
         xlab = F, 
         legend = "none"
)+
  stat_compare_means(size = 5,label.x =1.5,label.y = 30,symnum.args=list(cutpoints = c(0,0.0001,0.001, 0.01, 0.05, 1), symbols = c("****","***", "**", "*", "ns")),label = "p.signif",method = 'wilcox.test')         
#Figure 5E、H-TMB, MATH survival curves
library("ggplot2")
library("ggpubr")
library("survival")
data=read.delim(file.choose())
# Fit survival curves
fit <- survfit(Surv(futime, fustat) ~ Risktype, data = data) 
# Plot informative survival curves
pdf(file="KM.pdf",width=8,height=7)
ggsurvplot(fit, data = data,
           #title = "TCGA-LUAD",          
           pval = TRUE, pval.method = F,    
           legend.title = "Group",               
           legend.labs = c("H-TMB+L-Risk", "H-TMB+H-Risk","L-TMB+L-Risk","L-TMB+H-Risk"),
           palette = c("#FBEA30","#5698C4","#E46C2D","#5A9A50"),    
           risk.table = T,  
           cumevents = F,cumcensor = F,conf.int = F,  
           tables.height = 0.15, 
           tables.theme = theme_cleantable(),  
           tables.y.text = FALSE,   
           xlab = "Survival(year)",   
           break.time.by = 1,    
           ggtheme = theme_bw(),
           legend = c(0.8,0.85)
)
dev.off()

#Figure 6A-B-GSEA enrichment analysis
library(clusterProfiler)
library(tidyverse)
library(GseaVis)
df = read.table("logFC.txt",header = T)
df.id<-bitr(df$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
head(df.id)
easy.df<-merge(df,df.id,by="SYMBOL",all=F)
sortdf<-arrange(easy.df,desc(log2FoldChange))
gene.expr = sortdf$log2FoldChange
names(gene.expr) <- sortdf$ENTREZID
kk <- gseKEGG(gene.expr, organism = "hsa")
###result
sortkk<-arrange(kk,desc(enrichmentScore))
write.table(sortkk,"gsea_output.txt",sep = "\t",quote = F,col.names = T,row.names = F)
terms<-c('hsa03030',
         'hsa04110',
         'hsa04115',
         'hsa03420',
         'hsa03008',"hsa00020")
lapply(terms,function(x){
  gseaNb(object=kk,
         geneSetID=x,
         addPval=T,
         pvalX=0.75,pvalY=0.75,
         pCol='black',
         pHjust=0)
})->gseaList
pdf("mutilpathway_gsea.pdf",width=16,height=16)
cowplot::plot_grid(plotlist=gseaList,ncol=2,align='hv')	
dev.off()	 
#Figure 6C-D-GSVA enrichment analysis
library(dplyr)
library(GSEABase)
library(GSVA)
library(limma)
library(stringr)
library(ggprism)
library(ggplot2)
dat <- read.table("500tpms01A_log2.txt",sep = "\t",row.names = 1,header = T,stringsAsFactors = F,check.names = F)
geneSets <- getGmt("h.all.v7.5.1.symbols.gmt")  
GSVA_hall <- gsva(expr=as.matrix(dat),
                  gset.idx.list=geneSets,
                  mx.diff=T, 
                  kcdf="Gaussian", 
                  parallel.sz=4)
write.csv(GSVA_hall,file="GSVA_hall.csv")
ann <- read.table("group.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
ann<-arrange(ann,desc(group))
samplename<-match(ann$sample,colnames(GSVA_hall))
GSVA_hall<-GSVA_hall[,samplename]
group<-rep(c("Low","High"),times=c(250,250))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(GSVA_hall)
compare <- makeContrasts(High - Low, levels=design)
fit <- lmFit(GSVA_hall, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1, number=200)
dat_plot <- data.frame(id = row.names(Diff),
                       t = Diff$t)
dat_plot$id <- str_replace(dat_plot$id , "HALLMARK_","")
dat_plot$threshold=factor(ifelse(dat_plot$t>-2,ifelse(dat_plot$t>=2 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
dat_plot <- dat_plot %>% arrange(t)
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+ 
  coord_flip() + 
  scale_fill_manual(values = c('Up'= '#8BC7C3','NoSignifi'='#cccccc','Down'='#E2AB78')) + 
  geom_hline(yintercept = c(-2,2),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('t value of GSVA score, high-risk versus low-risk') + 
  guides(fill=F)+ 
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p
low1 <- dat_plot %>% filter(t < -2) %>% nrow
low0 <- dat_plot %>% filter( t < 0) %>% nrow()
high0 <- dat_plot %>% filter(t < 2) %>% nrow()
high1 <- nrow(dat_plot)
p <- p + geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
                   hjust = 0,color = 'black') + 
  geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'grey') +
  geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'grey') + 
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black') 
ggsave("gsva_bar.pdf",p,width = 7.35,height  = 8.5)
#Figure 6E-H-GO, KEGG enrichment analysis
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
go <- read.table("GO.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
go <- as.data.frame(go)
dotplot(go, split="ONTOLOGY",showCategory = 8,label_format=50)+ 
  facet_grid(ONTOLOGY~.,scale="free")+
  theme(panel.grid = element_blank())+
  theme(axis.title =element_text(size = 12, color = 'black'),
        axis.text.y =element_text(size = 12),
        legend.title=element_text(size=12))+
  scale_color_gradient(high="#F67B2F",low="#A3C7DA")
p <- ggplot(go, aes(Term, GeneRatio)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_size(range = c(2, 6)) +
  scale_color_gradient(high="#6DB9E5",low="#FFABAD") +
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        panel.grid = element_blank(), legend.key = element_blank()) +
  coord_flip() +
  labs(x = '', y = 'Enrichment Score')
p
p1=p + facet_grid(Category~., scale = 'free_y', space = 'free_y')
p1
ggsave(p1,file="re.pdf",width =7,height = 9)
GO2  <-   read.csv("GO-up.csv", header = T, sep = ",")
library(WeightedTreemaps)
library(RColorBrewer)
oct_coord <- list(
  x = sin(seq(0, 2, 2/8)*pi) * 1000 + 1000,
  y = cos(seq(0, 2, 2/8)*pi) * 1000 + 1000
)
GO2<- voronoiTreemap(
  data = GO2,
  levels = c("num", "Id"),
  cell_size = "num",
  shape = oct_coord,
  seed = 123)
pdf("GO-up.PDF",width = 5,height = 5) 
drawTreemap(GO2, label_size = 2.5, label_color = "black",
            color_palette = c(brewer.pal(5,"Set3"),brewer.pal(5,"Pastel1")) ,color_level = 2 )
dev.off() 

#Figure 7A-D-Estimate analysis
rm(list = ls())
options(stringsAsFactors = F)
sig_violin<-function(dat,leg,ylab,palette=ggsci::pal_lancet()(10)[3:4]){
  library(ggpubr)
  dat=na.omit(dat)
  colnames(dat)=c('group','gene')
  dat=dat[order(dat$group),]
  all.combn=combn(as.character(unique(dat$group)),2)
  my_comparisons=lapply(seq_len(ncol(all.combn)), function(i) all.combn[,i])
  p<-ggplot(dat, aes(x=group, y=gene,color=group)) + 
    geom_violin(trim=FALSE,fill="white")+
    geom_boxplot(aes(fill=group),color='black')+
    stat_compare_means(comparisons=my_comparisons,method="wilcox.test",label = "p.signif")+
    ylab(ylab)+xlab('')+labs(color=leg,fill=leg)+scale_color_manual(values = palette)+theme_classic()+
    scale_fill_manual(values = palette)
  return(p)
}
tcga.caf.cli=read.delim('ES.txt',sep='\t',header = T,check.names = F)
fig3a=list()
fig3a[[1]]<-sig_violin(dat = tcga.caf.cli[,c("Risktype","TumorPurity")],
                       leg = 'Risktype',ylab = 'TumorPurity',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3a[[1]]
#Figure 7E-CIBERSORT analysis
library(limma)
library(pheatmap)
library(ggpubr)
library(vioplot)
library(corrplot)
immFile="CIBERSORT-Results.txt"     
riskFile="risk.txt"             
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<0.05,]
data=as.matrix(immune[,1:(ncol(immune)-3)])
rownames(data)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data))
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
lowSample=rownames(risk)[risk[,"risk"]=="Low"]
highSample=rownames(risk)[risk[,"risk"]=="High"]
lowSameSample=intersect(row.names(data), lowSample)
highSameSample=intersect(row.names(data), highSample)
data=t(data[c(lowSameSample,highSameSample),])
conNum=length(lowSameSample)
treatNum=length(highSameSample)
rt=t(data)
pdf("vioplot1.pdf", height=8, width=12)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,63), ylim=c(min(rt), max(rt)+0.02),
     main="", xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")
#Figure 7F-ssGSE analysis
library(limma)
library(GSVA)
library(GSEABase)
library(ggpubr)
library(reshape2)
expFile="500tpms01A_log2.txt"        
gmtFile="immune.gmt"        
riskFile="cluster.txt"     
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=mat[rowMeans(mat)>0,]
geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())
ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
data=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(data), data)
write.table(ssgseaOut, file="immFunScore.txt", sep="\t", quote=F, col.names=F)
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=t(data[,group==0])
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)
cluster=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data),row.names(cluster))
data=data[sameSample,,drop=F]
cluster=cluster[sameSample,"Cluster",drop=F]
rt1=cbind(data, cluster)
data=melt(rt1,id.vars=c("Cluster"))
colnames(data)=c("Cluster","Type","Score")
data$Cluster=factor(data$Cluster, levels=c("Low","High"))
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
theme <- theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size=16), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18), 
        axis.line = element_line(size = 1),
        plot.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key = element_rect(fill = 'transparent'), 
        legend.background = element_rect(fill = 'transparent'), 
        legend.position = "top",
        legend.title = element_blank())
p=ggboxplot(data, x="Type", y="Score", fill = "Cluster",
            notch=T, outlier.shape = NA,
            xlab="",ylab="Score",add = "none", palette=c("#89C9C8","#F9BEBB"))
p=p+rotate_x_text(50)+theme
p1=p+stat_compare_means(aes(group=Cluster),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
pdf(file="immFunction2.pdf", width=14, height=8)
print(p1)
dev.off()
#Figure 7G-H-Comparative expression of immune molecules
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
expFile="500tpms01A_log2.txt"         
RisktypeFile="Risktype.txt"      
geneFile="gene.txt"          
#Read gene expression files and process the data
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#Read the gene list file to obtain the expression of immune checkpoint-related genes.
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(row.names(data), as.vector(gene[,1]))
data=t(data[sameGene,])
data=avereps(data)
Risktype=read.table(RisktypeFile, sep="\t", header=T, check.names=F, row.names=1)
sameSample=intersect(row.names(data),row.names(Risktype))
rt1=cbind(data[sameSample,],Risktype[sameSample,])
rt1=rt1[,c(sameGene,"Risktype")]
sigGene=c()
for(i in colnames(rt1)[1:(ncol(rt1)-1)]){
  if(sd(rt1[,i])<0.001){next}
  wilcoxTest=wilcox.test(rt1[,i] ~ rt1[,"Risktype"])
  pvalue=wilcoxTest$p.value
  if(wilcoxTest$p.value<0.05){
    sigGene=c(sigGene, i)
  }
}
sigGene=c(sigGene, "Risktype")
rt1=rt1[,sigGene]
rt1=melt(rt1,id.vars=c("Risktype"))
colnames(rt1)=c("Risktype","Gene","Expression")
group=levels(factor(rt1$Risktype))
rt1$Risktype=factor(rt1$Risktype, levels=c("Low","High"))
comp=combn(group,2)
my_comparisons=list()
for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
theme <- theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size=16), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18), 
        axis.line = element_line(size = 1),
        plot.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key = element_rect(fill = 'transparent'), 
        legend.background = element_rect(fill = 'transparent'), 
        legend.position = "top",
        legend.title = element_blank())
boxplot=ggboxplot(rt1, x="Gene", y="Expression", fill="Risktype",
                  xlab="",
                  ylab="Gene expression",
                  legend.title="Risktype",
                  width=0.8,
                  #outlier.shape = NA,
                  palette = c("#89C9C8","#F9BEBB") )+theme+
  rotate_x_text(50)+
  stat_compare_means(aes(group=Risktype),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
pdf(file="checkpoint.diff2.pdf", width=10, height=5)
print(boxplot)
dev.off()
#Figure 7I-Comparison of immune function
library(limma)
library(GSVA)
library(GSEABase)
library(pheatmap)
library(reshape2)
expFile="500tpms01A_log2.txt"         
gmtFile="imm.gmt"         
riskFile="Risk.txt"     
# Read the expression input file and process the input file
rt=read.table(expFile, header=T, sep="\t", check.names=F,row.names = 1)
rt=as.matrix(rt)
mat=avereps(rt)
mat=mat[rowMeans(mat)>0,]
geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())
#ssgsea analysis
ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
data=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(data), data)
write.table(ssgseaOut, file="immFunScore.txt", sep="\t", quote=F, col.names=F)
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=t(data[,group==0])
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3-\\4", rownames(data))
data=t(avereps(data))
# Read the risk file and get the high and low risks.
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
lowSample=row.names(risk[risk$risk=="Low",])
highSample=row.names(risk[risk$risk=="High",])
lowData=data[,lowSample]
highData=data[,highSample]
data=cbind(lowData, highData)
conNum=ncol(lowData)        
treatNum=ncol(highData)     
sampleType=c(rep(1,conNum), rep(2,treatNum))
# Variance analysis
sigVec=c()
for(i in row.names(data)){
  test=wilcox.test(data[i,] ~ sampleType)
  pvalue=test$p.value
  Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
  sigVec=c(sigVec, paste0(i, Sig))
}
row.names(data)=sigVec
#Heatmap Visualisation
Type=c(rep("Low risk",conNum), rep("High risk",treatNum))
Type=factor(Type, levels=c("Low risk", "High risk"))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf("heatmap1.pdf", width=8, height=4.6)
pheatmap(data,
         annotation=Type,
         color = colorRampPalette(c(rep("#AAD1E5",4), "white", rep("#F3A683",4)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=7,
         fontsize_row=7,
         fontsize_col=7)
dev.off()
#Figure 7J-Immune Cell Correlation
inputFile="cor.result.txt"       
data = read.table(inputFile, header=T, sep="\t", check.names=F)
#colour
p.col = c('#86C0CB',"#A3D2E2",'#D9E0E6','#F9BEBB','#F3A683')
fcolor = function(x,p.col){
  color = ifelse(x>0.8,p.col[1],ifelse(x>0.6,p.col[2],ifelse(x>0.4,p.col[3],
                                                             ifelse(x>0.2,p.col[4], p.col[5])
  )))
  return(color)}
p.cex = seq(2.5, 5.5, length=5)
fcex = function(x){
  x=abs(x)
  cex = ifelse(x<0.1,p.cex[1],ifelse(x<0.2,p.cex[2],ifelse(x<0.3,p.cex[3],
                                                           ifelse(x<0.4,p.cex[4],p.cex[5]))))
  return(cex)}
points.color = fcolor(x=data$pvalue,p.col=p.col)
data$points.color = points.color
points.cex = fcex(x=data$cor)
data$points.cex = points.cex
data=data[order(data$cor),]
xlim = ceiling(max(abs(data$cor))*10)/10        
pdf(file="1.pdf", width=10, height=7)     
layout(mat=matrix(c(1,1,1,1,1,0,2,0,3,0),nc=2),width=c(8,2.2),heights=c(1,2,1,2,1))
par(bg="white",las=1,mar=c(5,18,2,4),cex.axis=1.5,cex.lab=2)
plot(1,type="n",xlim=c(-xlim,xlim),ylim=c(0.5,nrow(data)+0.5),xlab="Correlation Coefficient",ylab="",yaxt="n",yaxs="i",axes=F)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col="#F5F5F5",border="#F5F5F5")
grid(ny=nrow(data),col="white",lty=1,lwd=2)
segments(x0=data$cor,y0=1:nrow(data),x1=0,y1=1:nrow(data),lwd=4)
points(x=data$cor,y = 1:nrow(data),col = data$points.color,pch=16,cex=data$points.cex)
text(par('usr')[1],1:nrow(data),data$Cell,adj=1,xpd=T,cex=1.5)
#pvalue
pvalue.text=ifelse(data$pvalue<0.001,'<0.001',sprintf("%.03f",data$pvalue))
redcutoff_cor=0
redcutoff_pvalue=0.05
text(par('usr')[2],1:nrow(data),pvalue.text,adj=0,xpd=T,col=ifelse(abs(data$cor)>redcutoff_cor & data$pvalue<redcutoff_pvalue,"red","black"),cex=1.5)
axis(1,tick=F)
par(mar=c(0,4,3,4))
plot(1,type="n",axes=F,xlab="",ylab="")
legend("left",legend=c(0.1,0.2,0.3,0.4,0.5),col="black",pt.cex=p.cex,pch=16,bty="n",cex=2,title="abs(cor)")
par(mar=c(0,6,4,6),cex.axis=1.5,cex.main=2)
barplot(rep(1,5),horiz=T,space=0,border=NA,col=p.col,xaxt="n",yaxt="n",xlab="",ylab="",main="pvalue")
axis(4,at=0:5,c(1,0.8,0.6,0.4,0.2,0),tick=F)
dev.off()
#Figure 7K-Gene-Cell Correlation
library(tidyverse)
library(ggplot2)
library(reshape2)
library(corrplot)
library(xCell)
exp<-read.table("EXP.txt",sep = "\t",row.names = 1,header = T,stringsAsFactors = F,check.names = F)
CIBERSORT<-read.table("CIBERSORT.txt",sep = "\t",row.names = 1,header = T,stringsAsFactors = F,check.names = F)
genelist<-c('GENE')
goal_exp<-filter(exp,rownames(exp) %in%genelist)
combine<-rbind(goal_exp,xcell)
comcor<-cor(t(combine))
comp<-cor.mtest(comcor,conf.level=0.95)
goalcor<-select(as.data.frame(comcor),genelist)%>%rownames_to_column(var="celltype")
goalcor<-filter(goalcor,!(celltype %in% genelist))
goalcor<-melt(goalcor,id.vars="celltype")
colnames(goalcor)<-c("celltype","Gene","correlation")
pval<-select(as.data.frame(pval),genelist)%>%rownames_to_column(var="celltype")
pval<-filter(pval,!(celltype %in% genelist))
pval<-melt(pval,id.vars="celltype")
colnames(pval)<-c("celltype","gene","pvalue")
final<-left_join(goalcor,pval,by=c("celltype"="celltype","Gene"="gene"))
final$sign<-case_when(final$pvalue<0.05 &final$pvalue>0.01 ~"*",
                      final$pvalue<0.01 &final$pvalue>0.001 ~"**",
                      final$pvalue<0.001 ~"***",
                      final$pvalue>0.05 ~"")
ggplot(data=final,aes(x=Gene,y=celltype))+
  geom_tile(aes(fill=correlation),colour="white",size=1)+
  scale_fill_gradient2(low="#2b8cbe",mid="white",high="#e41a1c")+
  geom_text(aes(label=sign),colour="black")+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45,hjust=1,size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  guides(fill=guide_legend(title="* p<0.05\n\n** p<0.01\n\n*** p<0.001\n\ncorrelation"))
ggsave("correlation.pdf",width=8,height=9)

#Figure 8A-D-Comparison of TIDE, IPS
rm(list = ls())
options(stringsAsFactors = F)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
es <- read.table("input.txt",sep = "\t",row.names = 1,header = T,stringsAsFactors = F,check.names = F)
Palette <- c('#FAB598','#7ac7e2')
comparisons = list(c('High','Low'))
p <- ggplot(data=es,aes(Risktype,TIDE,fill=Risktype))+
  geom_violin(alpha=1,width=0.9,
              position=position_dodge(width=0.8),
              size=1.75)+
  geom_boxplot(width=0.2,alpha=0.5,
               size=1.75,outlier.colour = NA)+
  scale_fill_manual(values = Palette)+
  theme_light() +
  theme(legend.position="none",
        panel.grid.major.x = element_line(size = 1),
        panel.grid.major.y = element_line(size = 1),
        panel.grid.minor.x = element_line(size = 1),
        panel.grid.minor.y = element_line(size = 1),
        axis.line = element_line(colour = 'black', size = 0.75)) + 
  theme(text = element_text(size=16)) + 
  ylab("TIDE")+xlab('Risktype')+
  stat_compare_means(comparisons = comparisons,label = "p.signif") ###or'p.format'
p
ggsave(p,filename = 'TIDE.pdf',width = 4.5,height=6)
#Figure 8E-submap heatmap
library(pheatmap)
heatmap.YlGnPe <- c("#3C79B3","#78a3cc","#B3CDE4","#CFF4D2","#F9E181")
cherry    <- "grey"
lightgrey <- "white"
#Import data
tmp <- matrix(c("submap-result"),
              nrow = 4,byrow = T,dimnames = list(c("Low_p","High_p","Low_b","High_b"),c("PR","SD","CR","PD")))
pheatmap(tmp, cellwidth = 30, cellheight = 30,
         cluster_rows = F,cluster_cols = F,
         color = heatmap.YlGnPe[5:1],
         gaps_row = 2,
         annotation_row = data.frame(pvalue=c("Nominal p value","Nominal p value","Bonferroni corrected","Bonferroni corrected"),row.names = rownames(tmp)),
         annotation_colors = list(pvalue=c("Nominal p value"=lightgrey,"Bonferroni corrected"=cherry)),
         filename = "heatmap_submap2.pdf")
#Figure 8F-TIDE analysis
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(survival)
  library(survminer)
  library(ggpubr)
  library(patchwork)
  library(viridis)
  library(ComplexHeatmap)
})
tcga.cli <- read.table("tcga.cli.txt")
tide <- read.table("tide.res.txt",sep = "\t",row.names = 1,header = T,stringsAsFactors = F,check.names = F)
tide <- tide[rownames(tcga.cli),] %>% cbind(data.frame(binRS=tcga.cli$binRS)) %>% rownames_to_column("Samples")
ctable <- table(tide$Responder,tide$binRS)
write.csv(tide,file="tide.csv")
tide <- read.table("Tide.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
p1 <- tide %>% 
  ggbarplot(x="Samples",y="TIDE",fill="Responder",color = "transparent",
            title = sprintf("Immune checkpoint inhibitors\np.fisher = %.3f",p.value)) + 
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.1),
        plot.title = element_text(hjust = 0.5, face = "bold")
  ) + 
  scale_fill_manual(labels=c("Non-responder","Responder"),values=c("#FAB598","#7ac7e2"))
p2 <- ctable %>% t %>% {round(100 * ./rowSums(.),2)} %>% as.data.frame() %>% 
  ggbarplot(x="Var1",y="Freq",fill="Var2",
            label = T,lab.pos="in",lab.col = "white",width = .4)  + 
  theme(legend.position = "none") +  
  scale_fill_manual(labels=c("Non-responder","Responder"),values=c("#FAB598","#7ac7e2")) +
  labs(x="NMBRS",y="Percent weight")
p <- p1 + p2 + plot_layout(ncol = 2, widths = c(2, 1))
ggsave("results/Fig.7E.pdf",p,width = 9,height = 5)
#Figure 8G-N-Immunotherapy cohort analysis
library(dplyr)
library(ggplot2)
library(ggsci)
rt1=read.delim('RE.txt',sep='\t',header = T,check.names = F)
dat = count(rt1,group,Responder)
dat = dat %>% group_by(group) %>% 
  summarise(Responder = Responder,n = n/sum(n))
dat$Response = factor(dat$Responder,levels = c("PD/SD","PR/CR"))
library(ggplot2)
pdf('Fig.pdf',width=5,height=8)
ggplot(dat,aes(group,n,fill= Responder))+
  geom_bar(stat="identity",position = position_stack())+
  scale_fill_manual(values = c("#FAB598","#7ac7e2"),label=c("PD/SD","CR/PR"))+
  scale_y_continuous(labels = scales::percent_format(scale = 100))+ #Percentage y-axis
  labs(x="",y="Percent Weight",
       fill="")+
  geom_text(aes(label = scales::percent(n)),position = position_stack(vjust = 0.5),size=5,color="black")+
  theme_bw()+
  theme(legend.position = "top",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"))
dev.off()
#ggviolion
library(ggplot2)
library(ggpubr)
rt1=read.delim('RE.txt',sep='\t',header = T,check.names = F,row.names = 1)
comparisons <- list(c("PD/SD","CR/PR")) 
ggviolin(rt1, x = "responder", y = "riskscore", fill = "responder", palette = c("npg"),
         add = "boxplot", add.params = list(fill = "white"), order = c("PD/SD","CR/PR")) + 
  stat_compare_means(comparisons = comparisons)
library("ggplot2")
library("ggpubr")
library("survminer")
library("survival")
data=read.delim(file.choose())
# Fit survival curves
fit <- survfit(Surv(OS.time1, OS) ~ group, data = data) 
# Plot informative survival curves
pdf(file="KM.pdf",width=8.5,height=8)
ggsurvplot(fit, data = data,
           title = "IMvigor210",          # Title
           pval = TRUE, pval.method = T,    
           legend.title = "Group",               
           legend.labs = c("High", "Low"),      
           palette = c("#FAB598","#7ac7e2"),    
           risk.table = T,  
           cumevents = F,cumcensor = F,conf.int = T,  
           tables.height = 0.15, 
           tables.theme = theme_cleantable(),  
           tables.y.text = FALSE,   
           xlab = "Survival(year)",  
           break.time.by = 1,    
           ggtheme = theme_bw(),
           legend = c(0.8,0.9)
)
dev.off()

#Figure 9A-Single-cell analysis of genes
load("Macrodata.RData")
library(Seurat)
FeaturePlot(Macrodata,features = c('ALG3',
                                   'SNRPA',
                                   'UCK2',
                                   'ZIC2',
                                   'ZNF490',
                                   'ADSS1','NT5E'),color <- c('lightgrey', '#FB8072','#DC050C'))
#Figure 9B
rm(list=ls())
# devtools::install_github("sajuukLyu/ggunchull", type = "source")
# install.packages("tidydr")
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggunchull)
library(tidydr)
library(ggsci)
library(Cairo)
load("sce.all.RData")
color = c(pal_d3("category20")(20),
          pal_d3("category20b")(20),
          pal_d3("category20c")(20),
          pal_d3("category10")(10))

color2=c("#e5192c","#3a77b7","#3cac4c","#813c93","#f36c24",
         "#37b8c3","#a54922","#6b7627","#28996b",
         "#965b6a","#e9148f","#595b5e",
         "#80d08a","#d29099","#f2e010")

DimPlot(sce.all, reduction = "tsne", group.by = "Cell_type",
        cols = color2,
        pt.size = 1,
        label = T,label.box = T
) 
ggsave(file="tsne.pdf",width = 12,height = 9)
df=sce.all@reductions$tsne@cell.embeddings %>%
  as.data.frame() %>%
  cbind(Cell_type=sce.all@meta.data$Cell_type)

head(df)
#tsne
CairoPDF("tsne_ggunchull.pdf", width=12, height=9) 
p_tsne=ggplot(df, aes(tSNE_1, tSNE_2, color=Cell_type,fill = Cell_type))+
  geom_point(size = 1) + 
  scale_color_manual(values=color2)+
  theme(panel.border = element_blank(), 
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        #panel.grid = element_blank(),
        # panel.grid.major = element_blank(), 
        # panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"),
        legend.title = element_blank(), 
        legend.key=element_rect(fill='white'), 
        legend.text = element_text(size=14), 
        legend.key.size=unit(0.6,'cm') )+
  stat_unchull(
    fill = "white",
    alpha = 0,
    show.legend = FALSE,
    nsm = 30,
    nbin = 150,
    sfac = 1.2)+
  guides(fill= guide_legend(override.aes = list(size = 4)))+ 
  scale_color_manual(values=color2);p_tsne
dev.off()
CairoPDF("tsne_ellipse.pdf", width=12, height=9) 
p_tsne1=ggplot(df, aes(tSNE_1, tSNE_2, color=Cell_type))+
  geom_point(size = 1) + 
  scale_color_manual(values=color2)+
  theme(panel.border = element_blank(), 
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        # panel.grid = element_blank(),
        # panel.grid.major = element_blank(), 
        # panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"),
        legend.title = element_blank(), #legend.title 
        legend.key=element_rect(fill='white'), 
        legend.text = element_text(size=14), #legend.size
        legend.key.size=unit(0.6,'cm') )+
  stat_ellipse(aes(x = tSNE_1, y = tSNE_2, fill = Cell_type),
               geom = "polygon",
               linetype=2, 
               alpha = 0.05,  
               show.legend = FALSE, 
               level = 0.93)+ 
  guides(fill= guide_legend(override.aes = list(size = 4)))+ 
  scale_color_manual(values=color2); p_tsne1 
dev.off()
CairoPDF("tsne_theme_dr.pdf", width=12, height=9) 
p_tsne2<-p_tsne1+theme_dr(xlength = 0.22, ylength = 0.22,
                          arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank());p_tsne2
dev.off()
Cell_type_position <- df %>%
  group_by(Cell_type) %>%
  dplyr::summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2))
library(ggrepel)
CairoPDF("tsne_label.pdf", width=12, height=9) 
p_tsne3<-p_tsne2+#geom_point(aes(color=factor(Cell_type)), size=3)+
  geom_label_repel(data = Cell_type_position,aes(label=Cell_type), 
                   fontface="bold",point.padding=unit(0.1, "lines"))+theme(legend.position = "none");p_tsne3
dev.off()
#Figure 9C-Comparison of riskscores
data=read.delim('input.txt',sep='\t',header = T,check.names = F)
library(ggplot2)
library(ggpubr)
library(reshape2)
p=ggviolin(data, x="group_copykat", y="NMBRS", 
           color = "group_copykat",
           fill="group_copykat",
           palette =c("#9AC9D9","#E17877"),
           add = "boxplot",
           add.params = list(color="white"),
           xlab = F,
           legend = "none"
)+
  stat_compare_means(size = 5,label.x =1.5,label.y = 1.2,symnum.args=list(cutpoints = c(0,0.0001,0.001, 0.01, 0.05, 1), symbols = c("****","***", "**", "*", "ns")),label = "p.signif",method = 'wilcox.test')         
p
ggsave(p,file="output.pdf",width = 5,height = 8)
#Figure 9D-Single-cell GSEA
library(clusterProfiler)
library(tidyverse)
library(GseaVis)
df = read.table("logFC.txt",header = T)
df.id<-bitr(df$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
head(df.id)
easy.df<-merge(df,df.id,by="SYMBOL",all=F)
sortdf<-arrange(easy.df,desc(log2FoldChange))
gene.expr = sortdf$log2FoldChange
names(gene.expr) <- sortdf$ENTREZID
kk <- gseKEGG(gene.expr, organism = "hsa")
sortkk<-arrange(kk,desc(enrichmentScore))
write.table(sortkk,"gsea_output.txt",sep = "\t",quote = F,col.names = T,row.names = F)
terms<-c('hsa00760',
         'hsa00240',
         'hsa05230',
         'hsa01232',
         'hsa04066',"hsa03010","hsa05202","hsa04510","hsa04151")
pdf("gsea_mutilpathway.pdf",width=14,height=8)
gseaNb(object=kk,
       geneSetID=terms,
       subPlot=2,
       curveCol=ggsci::pal_npg()(10),
       addPval=T,
       pvalX=1,pvalY=1)
dev.off()
#Figure 9E-I-Analysis of intercellular communication
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SingleR)
library(cowplot)
library(tidyverse)
library(CellChat)
load(sce.all.RData)
af=sce.all
cellchat = createCellChat(object = af,meta =af@meta.data, group.by = "Cell_type")
levels(cellchat@idents) 
groupSize <- as.numeric(table(cellchat@idents)) 
CellChatDB = CellChatDB.human 
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use = CellChatDB
#cellchat@netP$pathways
#[1] "TGFb"       "NRG"        "PDGF"       "CCL"        "CXCL"       "MIF"        "IL2"        "IL6"        "IL10"       "IL1"        "CSF"       
#[12] "IL16"       "IFN-II"     "LT"         "LIGHT"      "FASLG"      "TRAIL"      "BAFF"       "CD40"       "VISFATIN"   "COMPLEMENT" "PARs"      
#[23] "FLT3"       "ANNEXIN"    "GAS"        "GRN"        "GALECTIN"   "BTLA"       "BAG"       
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB = CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(object=cellchat,raw.use = TRUE)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
par(mfrow = c(1,2), xpd=TRUE)    
par(mfrow = c(1,1), xpd=TRUE)
#Number
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
#weights
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
par(mfrow = c(2,3), xpd=TRUE)        
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
cellchat@netP$pathways
pathways.show <- c("EGF") 
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
par(mfrow=c(1,1))
pdf(file="EGF.pdf", width=10, height=6)
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()
pdf(file="EGF.pdf", width=10, height=6)
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()

#Figure 10A-Differences in gene expression
library(tidyverse)
library(introdataviz)
library(ggpubr)
library(ggsci)
library(ggprism)
library(scales)
library(patchwork)
library(ggplot2)
dd2<- read_tsv("input.txt") %>% pivot_longer(- lncRNA_signature)
colors = c("#e64b35","#4dbbd5")
a <- ggplot(dd2,aes(x = name, y = value,fill= lncRNA_signature))+
  geom_split_violin(trim =F,color = NA,adjust = 1.5)+
  guides(fill=guide_legend(title="group"))+
  scale_fill_manual(values = colors)+
  stat_summary(fun.data = "mean_sd",position=position_dodge(0.15),geom = "errorbar",width = .1) +
  stat_summary(fun = "mean", geom = "point", position=position_dodge(0.15),show.legend = F)+
  stat_compare_means(aes(group = lncRNA_signature), label = "p.signif",label.y=12, method="t.test")+
  # label = p.signif & p.format
  labs(x=NULL,y=NULL)+
  theme(axis.title = element_blank(),
        axis.text.x=element_text(angle =45,hjust =1,vjust =1,color="black",size = 10,margin = margin(b =2)),
        axis.text.y=element_text(color="black",size = 10,margin = margin(r =1)),
        panel.background = element_rect(fill = NA,color = NA),
        panel.grid.minor= element_line(size=0.2,color="#e5e5e5"),
        panel.grid.major = element_line(size=0.2,color="#e5e5e5"),
        panel.border = element_rect(fill=NA,color="black",size=1,linetype="solid"),
        legend.key=element_blank(),
        legend.title = element_text(color="black",size=10),
        legend.text = element_text(color="black",size=8),
        legend.spacing.x=unit(0.1,'cm'),
        legend.key.width=unit(0.5,'cm'),
        legend.key.height=unit(0.5,'cm'),
        legend.box.background=element_rect(colour = "black"), 
        legend.position = c(1,0), legend.justification = c(1,0),
        legend.background=element_blank())
ggsave(a,file="Exp.pdf",height = 6,width = 12)
#Figure 10B、D、E-KM curve
library("ggplot2")
library("ggpubr")
library("survminer")
library("survival")
data=read.delim(file.choose())
# Fit survival curves
fit <- survfit(Surv(OS.time, OS) ~ Risktype, data = data) 
# Plot informative survival curves
pdf(file="KMM.pdf",width=8,height=6)
ggsurvplot(fit, data = data,
           title = "Gene",       
           pval = TRUE, pval.method = T,    
           legend.title = "Group",               
           legend.labs = c("High", "Low"),      
           palette = c("#E64B35","#4DBBD5"),    
           risk.table = F,  
           cumevents = F,cumcensor = F,conf.int = T,  
           tables.height = 0.15, 
           tables.theme = theme_cleantable(),  
           tables.y.text = FALSE,   
           xlab = "Survival(year)",   
           break.time.by = 1,    
           ggtheme = theme_bw(),
           legend = c(0.8,0.9)
)
dev.off()
#Figure 10C-ROC diagnostic curve
library(pROC)
library(limma)
library(tidyverse)
gene="Gene"
rt=read.table(file ="tpms.txt" , header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=as.data.frame(data)
exp_data_T = data%>% dplyr::select(str_which(colnames(.), "-01A")) 
nT = ncol(exp_data_T) 
exp_data_N = data%>% dplyr::select(str_which(colnames(.), "-11A"))
nN = ncol(exp_data_N) 
data= cbind(exp_data_N, exp_data_T)
data=avereps(data)
data=t(data[gene,,drop=F])
data=as.data.frame(data)
group=sapply(strsplit(rownames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
group=as.data.frame(group)  
data$group=group$group
y=data$group
roc1=roc(y, as.numeric(data[,gene]))
ci1=ci.auc(roc1, method="bootstrap")
ciVec=as.numeric(ci1)
pdf(file=paste0("ROC.",gene,".pdf"), width=5, height=5)
plot(roc1, print.auc=TRUE, col="#6DB9E5", legacy.axes=T, main=gene)
text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="#6DB9E5")
dev.off()
