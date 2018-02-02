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
human<-Merge_GeneCount(data_dir[1])
human.sample.names<-colnames(human)
human.sample.names<-unlist(lapply(human.sample.names,function(x)return(str_replace_all(x,'-','_'))))

#  Extract the name 
#  consistent name

human.sample.names<-colnames(human)
human.sample.names<-unlist(lapply(human.sample.names,function(x)return(str_replace_all(x,'-','_'))))
human.sample.names<-unlist(lapply(human.sample.names,function(x){
  v<-str_split(x,'_')[[1]]
  if((str_detect(x,'^hc'))&(length(v)==6)){
    y<-str_split(x,'_',2)[[1]]
    z<-str_split(y[1],'')[[1]]
    if(length(z)==9){
      k<-str_c(str_c(z[1:5],collapse = ''),
               str_c(str_c(z[6:7],collapse = ''),'um'),str_c(z[8:9],collapse = ''),sep='_')
      k<-str_c(k,'_',y[2])
      return(k)
    }
    
    if(length(z==8)){
      k<-str_c(str_c(z[1:5],collapse = ''),
               str_c(str_c(z[6],collapse = ''),'um'),str_c(z[7:8],collapse = ''),sep='_')
      k<-str_c(k,'_',y[2])
      return(k)
    }
  }else{
    return(x)
  }
}))


human.sample.names<-unlist(lapply(human.sample.names,function(x){
  v<-str_split(x,'_')[[1]]
  if(length(v)==3){
    v<-str_c(str_c('hc001','20um',sep = '_'),str_c(x,collapse = '_'),sep = '_')
    return(v)
  }else{
    return(x)
  }
}))


human.sample.names<-unlist(lapply(human.sample.names,function(x){
  if(str_detect(x,'^[A-Z]{1}')){
    v<-str_c('hc006','_','20um')
    v<-str_c(v,'_',x)
    return(v)
  }else{
    return(x)
  }
}))

human.sample.names<-unlist(lapply(human.sample.names,function(x){
  if(str_detect(x,'^[0-9]')){
    len<-str_extract(x,'^[0-9]{1,}um')
    v<-str_split(x,'um')[[1]][2]
    v<-str_c(str_c('hc009','_',len),'_',v)
    return(v)
  }else{
    return(x)
  }
}))


human.sample.names<-unlist(lapply(human.sample.names,function(x){
  if(str_detect(x,"shoutiao")){
    v<-str_split(x,'6um')[[1]]
    v<-str_c(v[1],'6um',v[2],sep = '_')
    return(v)
  }else{
    return(x)
  }
}))

#   sample identity
groups<-unlist(lapply(human.sample.names,function(x){
  v<-str_split(x,'_')[[1]]
  return(v[1])
}))

#   each cell size 
cell.size<-unlist(lapply(human.sample.names,function(x){
  v<-str_split(x,'_')[[1]]
  return(v[2])
}))

#  groups len
groups.size<-unlist(lapply(human.sample.names,function(x){
  v<-str_split(x,'_')[[1]]
  v<-str_c(v[1],'_',v[2])
  return(v)
}))



#  Analysis
colnames(human)<-human.sample.names
human.pbmc<-CreateSeuratObject(human,min.cells =100,min.genes = 3 )
human.pbmc@meta.data[,'group.size']<-groups.size[human.sample.names%in%human.pbmc@cell.names]
human.pbmc@meta.data[,'size']<-cell.size[human.sample.names%in%human.pbmc@cell.names]
human.pbmc@meta.data[,'group']<-groups[human.sample.names%in%human.pbmc@cell.names]

#   Prepare data
human.pbmc<-NormalizeData(human.pbmc)
human.pbmc<-ScaleData(human.pbmc)
human.pbmc<-FindVariableGenes(human.pbmc)  # Find highly variable gene

human.variable.gene100<-rownames(x=head(x=human.pbmc@hvg.info,100))


#  Basic exploration of data
#  Violin plot for highest variable genes ,maybe
VlnPlot(human.pbmc,human.variable.gene100[1:4],y.lab.rot = 90)


#   Explore
# Plot two cells against each other
# Set do.ident=TRUE to identify outliers by clicking on points (ESC to exit after)
par(mfrow=c(2,2))
CellPlot(human.pbmc,cell1 = human.pbmc@cell.names[1],cell2 = human.pbmc@cell.names[2])
CellPlot(human.pbmc,cell1 = human.pbmc@cell.names[3],cell2 = human.pbmc@cell.names[4])

GenePlot(human.pbmc,human.variable.gene100[1],human.variable.gene100[2])
GenePlot(human.pbmc,human.variable.gene100[1],human.variable.gene100[2],
         cell.ids =WhichCells(object = human.pbmc,ident = 'hc009'))


#  Perform linear dimensional reduction (PCA)

human.pbmc<-RunPCA(object = human.pbmc,do.print = FALSE)
DimPlot(object = human.pbmc,reduction.use = 'pca',pt.size = 4)
DimPlot(object = human.pbmc,reduction.use = 'pca',pt.size = 4,group.by = 'size')
PrintPCA(human.pbmc,c(1,2))


#   Visualize PCA genes
VizPCA(object =human.pbmc,pcs.use = 1:2,do.balanced = TRUE,num.genes = 10,font.size =0.8)
PCElbowPlot(human.pbmc)

# Draw a heatmap where both cells and genes are ordered by PCA score
# Options to explore include do.balanced (show equal # of genes with +/- PC scores), 
#and use.full (after projection, see below)

DimHeatmap(object = human.pbmc,reduction.type = 'pca',do.balanced = TRUE)



#  Determine statistically significant principal components
human.pbmc<-JackStraw(human.pbmc,do.print=FALSE)

# The jackStraw plot compares the distribution of P-values for each PC with a uniform distribution (dashed line)
# 'Significant' PCs will have a strong enrichment of genes with low p-values (solid curve above dashed line)
# In this case PC1-15 are strongly significant
JackStrawPlot(human.pbmc,PCs = 1:15)


#  Optional step, grow gene list through PCA projection

# Previous analysis was performed on < 400 variable genes. To identify a larger gene set that may drive
# biological differences, but did not pass our mean/variability thresholds, we first calculate PCA scores
# for all genes (PCA projection)

human.pbmc<-ProjectPCA(human.pbmc,do.print = TRUE) #projects this onto the entire dataset (all genes)

# Visualize the full projected PCA, which now includes new genes which 
# were not previously included (use.full=TRUE)
PCHeatmap(object = human.pbmc,use.full = TRUE,do.balanced = TRUE)

#  # Choose significant genes for PC1-15, allow each PC to contribute a max of 200 genes
# (to avoid one PC swamping the analysis)

pc.sig.genes<-PCASigGenes(object = human.pbmc,pcs.use = 1:15,pval.cut = 1e-5,max.per.pc = 200)


#   Run Non-linear dimensional reduction (tSNE)
#  Based on the tNSE plot.It seems that we can not split the dataset ???
human.pbmc<-RunTSNE(object = human.pbmc,reduction.use = 'pca',dims.use = 1:15,
                     genes.use = pc.sig.genes)
DimPlot(object = human.pbmc,reduction.use = 'tsne',pt.size = 3)
DimPlot(object = human.pbmc,reduction.use = 'tsne',pt.size = 3,group.by = 'size')
DimPlot(object = human.pbmc,reduction.use = 'tsne',pt.size = 3,group.by = 'group.size')


#  Feature plot
Features<-pc.sig.genes[1:6]
FeaturePlot(human.pbmc,features.plot = Features,
            no.legend = FALSE,dark.theme = TRUE)

#  Find markers
human.pbmc<-SetIdent(human.pbmc,ident.use ='newcells',cells.use = human.pbmc@cell.names[1:8])
markers<-FindMarkers(human.pbmc,ident.1 = 'newcells',genes.use = pc.sig.genes)



library(DESeq2)
human.condition.matrix<-as.matrix(human[,human.condition.sample])
#human.condition.matrix<-human.condition.matrix[which(unlist(apply(human.condition.matrix,1,sum))>0),]
condition<-factor(condition)
coldata<-data.frame(row.names = human.condition.sample,condition)
dds<-DESeqDataSetFromMatrix(countData = human.condition.matrix,colData =coldata,
                            design = ~condition)

X<-apply(human.condition.matrix,1,function(x){
  if(all(x!=0)){
    return(x)
  }else{
    return(x+1)
  }
})

xdds<-DESeqDataSetFromMatrix(countData = X,colData =coldata,
                            design = ~condition)

xxds<-DESeq(xdds)
plotDispEsts(xdds, main="Dispersion plot")

