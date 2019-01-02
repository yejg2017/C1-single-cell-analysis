library(Seurat)
library(stringr)
source("tools.R")

data_dir="./a-9_geneCounts.csv"
geneCounts=read.csv(data_dir,header = TRUE)[-c(1:5),]
genes=as.character(geneCounts$gene)

geneCounts=geneCounts[!duplicated(genes),]
rownames(geneCounts)=as.character(geneCounts$gene)           
geneCounts=geneCounts[,-1]

### each geneCounts in variable
geneSum=apply(geneCounts,1,sum)
var.genes=rownames(geneCounts)[geneSum>3000]

### create Seurat Oject
obj=CreateSeuratObject(raw.data = geneCounts[var.genes,],project = "Mul.Cell")
obj=NormalizeData(obj,display.progress =FALSE)
obj=ScaleData(obj,display.progress = FALSE)
obj=FindVariableGenes(obj,do.plot = FALSE,display.progress = FALSE)
obj@meta.data$sample=paste0("a.",1:nrow(obj@meta.data))


### cv2 model for gene selection
#v=get.variable.genes(geneCounts,min.cv2 =2,do.plot = FALSE)
#var.genes<-as.character(rownames(v))[v$p.adj<0.01]


obj=RunPCA(obj,do.print = FALSE,pcs.compute =3)
obj=RunTSNE(obj,reduction.use = "pca",dims.use = 1:3,perplexity=2)


p2=DimPlot(obj,reduction.use = "tsne",do.return = T,plot.title = "tSNE",pt.size = 2,group.by = "sample")
p1=DimPlot(obj,reduction.use = "pca",do.return = T,plot.title = "PCA",pt.size = 2,group.by = "sample")
plot_grid(p1,p2,ncol = 1)


### caculate avg expression
avg.gene=AverageExpression(obj,genes.use =rownames(obj@raw.data))
avg.gene=avg.gene[order(avg.gene$Mul.Cell,decreasing = TRUE),,drop=FALSE]

#markers <- FindAllMarkers(object =obj)
#markers %>% group_by(cluster) %>% top_n(20, avg_logFC) -> top20

markers<-str_to_title(c('TP63','TJP1','KRT19','GJA1','ITGB4','ITGB6','ITGB1','KRT14','KRT12',
           'CDH1','CLDN4','BMP2','TFAP2A','EHF','KLF4',
           'KLF5','KRT3','GJA1','CDH1','CTNNB1'))
DoHeatmap(object = obj,genes.use =rownames(avg.gene)[1:20],col.high = "red",col.low = "lightblue")
