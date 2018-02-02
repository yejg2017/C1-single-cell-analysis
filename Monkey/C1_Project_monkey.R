library(Seurat)
library(data.table)
library(NMF)
library(rsvd)
library(Rtsne)
library(ggplot2)
library(cowplot)
library(sva)
library(igraph)
library(cccd)
library(KernSmooth)
library(beeswarm)
library(stringr)
library(formatR)
source('tools.R')


#   Read data
data_dir<-paste('../C1data/',list.files('../C1data/'),'/',sep='')
monkey<-Merge_GeneCount(data_dir[2])

monkey.sample.names<-colnames(monkey)
monkey.sample.names<-unlist(lapply(monkey.sample.names,function(x)return(str_replace_all(x,'-','_'))))
colnames(monkey)<-monkey.sample.names
sample.groups<-unlist(lapply(monkey.sample.names,function(x)return(str_split(x,"_")[[1]][1])))
sample.groups.len<-unlist(lapply(monkey.sample.names,function(x){
  v<-str_split(x,'_')[[1]]
  v<-str_c(v[1],v[2],sep = '_')
  return(v)
}))  # include length 

cell.size<-unlist(lapply(monkey.sample.names,function(x)return(str_split(x,"_")[[1]][2])))
#  sample len
sample.len<-unlist(lapply(sample.groups.len,function(x)return(str_split(x,'_')[[1]][2])))


#   Create Seurat object
monkey.pbmc<-CreateSeuratObject(monkey,min.cells = 3,min.genes =100,is.expr = 1,
                                names.field = 1,names.delim = '_')  # remove those genes all count is < 3
monkey.pbmc@meta.data[,'group.len']<-sample.groups.len[monkey.sample.names%in%monkey.pbmc@cell.names]
monkey.pbmc@meta.data[,'length']<-sample.len[monkey.sample.names%in%monkey.pbmc@cell.names]

monkey.pbmc<-NormalizeData(monkey.pbmc)
monkey.pbmc<-ScaleData(monkey.pbmc)
monkey.pbmc<-FindVariableGenes(monkey.pbmc)  # Find highly variable gene

monkey.variable.gene100<-rownames(x=head(x=monkey.pbmc@hvg.info,100))

#  Basic exploration of data
#  Violin plot for highest variable genes ,maybe
VlnPlot(monkey.pbmc,monkey.variable.gene100[1:4])

#   Explore
# Plot two cells against each other
# Set do.ident=TRUE to identify outliers by clicking on points (ESC to exit after)
par(mfrow=c(2,2))
CellPlot(monkey.pbmc,cell1 = monkey.pbmc@cell.names[1],cell2 = monkey.pbmc@cell.names[2])
CellPlot(monkey.pbmc,cell1 = monkey.pbmc@cell.names[3],cell2 = monkey.pbmc@cell.names[4])

GenePlot(monkey.pbmc,monkey.variable.gene100[1],monkey.variable.gene100[2])
GenePlot(monkey.pbmc,monkey.variable.gene100[1],monkey.variable.gene100[2],
         cell.ids =WhichCells(object = monkey.pbmc,ident = 'mkc001'))


#  Perform linear dimensional reduction (PCA)

monkey.pbmc<-RunPCA(object = monkey.pbmc,do.print = FALSE)
DimPlot(object = monkey.pbmc,reduction.use = 'pca',pt.size = 4)
DimPlot(object = monkey.pbmc,reduction.use = 'pca',pt.size = 4,group.by = 'length')
PrintPCA(monkey.pbmc,c(1,2))


#   Visualize PCA genes
VizPCA(object =monkey.pbmc,pcs.use = 1:2,do.balanced = TRUE,num.genes = 10,font.size =0.8)
PCElbowPlot(monkey.pbmc)

# Draw a heatmap where both cells and genes are ordered by PCA score
# Options to explore include do.balanced (show equal # of genes with +/- PC scores), 
#and use.full (after projection, see below)

DimHeatmap(object = monkey.pbmc,reduction.type = 'pca',do.balanced = TRUE,cells.use = 20)


#  Determine statistically significant principal components
monkey.pbmc<-JackStraw(monkey.pbmc,do.print=FALSE)

# The jackStraw plot compares the distribution of P-values for each PC with a uniform distribution (dashed line)
# 'Significant' PCs will have a strong enrichment of genes with low p-values (solid curve above dashed line)
# In this case PC1-8 are strongly significant
JackStrawPlot(monkey.pbmc,PCs = 1:9)


#  Optional step, grow gene list through PCA projection

# Previous analysis was performed on < 400 variable genes. To identify a larger gene set that may drive
# biological differences, but did not pass our mean/variability thresholds, we first calculate PCA scores
# for all genes (PCA projection)

monkey.pbmc<-ProjectPCA(monkey.pbmc,do.print = TRUE) #projects this onto the entire dataset (all genes)

# Visualize the full projected PCA, which now includes new genes which 
# were not previously included (use.full=TRUE)
PCHeatmap(object = monkey.pbmc,use.full = TRUE,do.balanced = TRUE)

#  # Choose significant genes for PC1-9, allow each PC to contribute a max of 200 genes
# (to avoid one PC swamping the analysis)

pc.sig.genes<-PCASigGenes(object = monkey.pbmc,pcs.use = 1:9,pval.cut = 1e-5,max.per.pc = 200)


#  Now redo the PCA analysis with the new gene list
monkey.pbmc<-RunPCA(object = monkey.pbmc,pc.genes = pc.sig.genes,do.print = FALSE)
DimPlot(object = monkey.pbmc,reduction.use = 'pca',pt.size = 4)

monkey.pbmc<-JackStraw(monkey.pbmc,num.replicate = 200,do.print = FALSE)
JackStrawPlot(monkey.pbmc,PCs = 1:15)



#   Run Non-linear dimensional reduction (tSNE)
#  Based on the tNSE plot.It seems that we can not split the dataset ???
monkey.pbmc<-RunTSNE(object = monkey.pbmc,reduction.use = 'pca',dims.use = 1:12,
                     genes.use = pc.sig.genes)
DimPlot(object = monkey.pbmc,reduction.use = 'tsne',pt.size = 3)
DimPlot(object = monkey.pbmc,reduction.use = 'tsne',pt.size = 3,group.by = 'length')
DimPlot(object = monkey.pbmc,reduction.use = 'tsne',pt.size = 3,group.by = 'group.len')



#  Cluster the cells
# Perform spectral density clustering on single cells

#  can not clutser(only one cluster)
# monkey.pbmc<-DBClustDimension(monkey.pbmc,reduction.use = 'tsne',G.use = 4,set.ident = TRUE)
# TSNEPlot(monkey.pbmc,pt.size = 4)

#  Find clutser
# monkey.pbmc<-FindClusters(monkey.pbmc,genes.use = pc.sig.genes,dims.use = 1:12,plot.SNN = TRUE,
#                           save.SNN = TRUE)


