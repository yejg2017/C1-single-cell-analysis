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


#   Read data(txt)
Load_data<-function(data_dir,condition.sample_dir=NULL){
  dat<-read.table(data_dir)
  
  if(!is.null(condition.sample_dir)){
    condition.sample.data<-read.table(condition.sample_dir)
    condition.sample<-as.character(condition.sample.data[,'sample'])
    condition<-as.character(condition.sample.data[,'condition'])
    
    return(list(data.frame(dat[,condition.sample]),condition))
  }else{
    return(dat)
  }
}


#  Not control condition

##   Analysis:
## QA
human.only.pro<-Load_data(data_dir = './human.txt')

important.genes<-c('ITGB4','ABCB5','KRT19','ACTB','KRT12','KRT5',
                   'GAPDH','KRT3','PAX6','WNT7A','KRT14','TP63','KRT10')

DESeq_SeuratObj<-function(X,cells.1=NULL,cells.2=NULL,min.cells=0,min.genes=0,condition=NULL,
                          FindVar=TRUE,condition.name=NULL,simple.rm=TRUE,DESq=TRUE){
  
  #  data frame or matrix of raw data
  #  condition : control condtion:Positive,Negative,Tube_1,Tube_2
  #  condition.name: when condition is not NULL,the condition.name is the respective sample names
  #  FindVar: whether is use FindVariableGenes function
  if(simple.rm){
    X<-Simplify_Select_Genes(X)
  }else{
    X<-X
  }
  
  if(!is.null(condition)){
    X<-X[,condition.name]
    con.pbmc<-CreateSeuratObject(X,min.cells = min.cells,min.genes = min.genes)
    con.pbmc<-NormalizeData(con.pbmc)
    con.pbmc<-ScaleData(con.pbmc)
    if(FindVar){
      con.pbmc<-FindVariableGenes(con.pbmc,do.plot = TRUE)
    }
    con.pbmc@meta.data[,'condition']<-condition[colnames(X)%in%rownames(con.pbmc@meta.data)]
    con.pbmc@meta.data[,'group']<-unlist(lapply(rownames(con.pbmc@meta.data),function(x)return(str_split(x,'_')[[1]][1])))
    con.pbmc@meta.data[,'cell.size']<-unlist(lapply(rownames(con.pbmc@meta.data),function(x)return(str_split(x,'_')[[1]][2])))

    #   Change the ident
    con.pbmc<-SetIdent(con.pbmc,cells.use = rownames(con.pbmc@meta.data),
                       ident.use = condition[colnames(X)%in%rownames(con.pbmc@meta.data)])
    
    if(DESq){
      con.DESq<-DESeq2DETest(con.pbmc,cells.1 = WhichCells(con.pbmc,ident = cells.1),
                             cells.2 = WhichCells(con.pbmc,ident = cells.2))
      return(list(con.pbmc,con.DESq))
    }else{
    return(con.pbmc)
    }
  }else{
    pbmc<-CreateSeuratObject(X,min.cells = min.cells,min.genes = min.genes)
    pbmc<-NormalizeData(pbmc)
    pbmc<-ScaleData(pbmc)
    
    if(FindVar){
      pbmc<-FindVariableGenes(pbmc,do.plot = TRUE)
    }
    
    pbmc@meta.data[,'group']<-unlist(lapply(rownames(pbmc@meta.data),function(x)return(str_split(x,'_')[[1]][1])))
    pbmc@meta.data[,'cell.size']<-unlist(lapply(rownames(pbmc@meta.data),function(x)return(str_split(x,'_')[[1]][2])))
    
    if(DESq){
      all.DESq<-DESeq2DETest(pbmc,cells.1 = WhichCells(pbmc,ident = cells.1),
                             cells.2 = WhichCells(pbmc,ident = cells.2))
      return(list(pbmc,all.DESq))
    }else{
    return(pbmc)
    }
  }
}

#  Create Seurate
##  DESeq analysis
#human.all<-DESeq_SeuratObj(X=human.only.pro,DESq = FALSE)

#  Here we try the DESeq(gene differential expression analysis) 

# only select the cells contain 10 genes expressed at least,select the genes must be expressed in two cells at least
human.all.DESeq<-DESeq_SeuratObj(X=human.only.pro,DESq = FALSE,min.cells = 10,min.genes = 2)
# check the important genes whether is still in the data after QA
important.genes[!important.genes%in%rownames(human.all.DESeq@raw.data)]  # KRT12 miss



# human.all.TTest<-DiffTTest(human.all,cells.1 = WhichCells(human.all,ident = 'hc001'),
#                            cells.2 = WhichCells(human.all,ident = 'hc006'))



#  Create sub Seurate object
Sub_Seurat<-function(object,ident.1,ident.2){
  #  ident.1:  the first ident selected in the object
  #  ident.2:  the second ident selected  in the object
  cells1<-WhichCells(object=object,ident = ident.1)
  cells2<-WhichCells(object=object,ident = ident.2)
  
  cells<-c(cells1,cells2)
  Sub.ident<-object@ident[object@cell.names%in%cells]
  Subobj<-SubsetData(object = object,cells.use = cells,ident.use = Sub.ident)
  
  return(Subobj)
}


#  Let’s examine the effective library size 
#(estimated by the average read count per gene) for the different cells.


#   Plot
###  why I can not add colors???
###  Barplot
Group_Bar<-function(count.data,group,Gv.return=FALSE){
  
  if(dim(count.data)[2]!=length(group)){
    stop('Invalid count data')
  }
  nGroups<-length(unique(group))
  nGene<-dim(count.data)[1]
  uniq.group<-unique(group)
  
  Group.per.gene.cov<-lapply(uniq.group,function(x){
    counts<-rowSums(count.data[,x%in%group])/nGene
    counts<-data.frame(counts)
    colnames(counts)<-x
    return(counts)
  })
  
  Groups.converage<-as.data.frame(Group.per.gene.cov)
  barplot(as.matrix(Groups.converage),beside = TRUE,ylab="Counts per gene")
  #axis(side=1,col =as.integer(as.factor(nGroups)),each=nGene))
  
  if(Gv.return){
    return(Groups.converage)
  }
}

#  some genes exactly have high expression in all sample excep control condition
#  but the most gene expressed not sinificant

#  Return the genes after DESqe analysis
DESeq_Genes<-function(p.data,p.val=0.05,num.genes=20,Order=TRUE,bar=TRUE){
  #  p.data: the data from DESeq analysis:genes,p.val
  #  p.val :the p.val cutoff
  #  Order: Wheter return the genes after p.val rank
  
  p.genes<-data.frame(na.omit(p.data))
  
  if(Order){
    p.genes<-p.genes[order(p.genes$p_val),,drop=FALSE]
  }
  
  genes<-rownames(p.genes)[which(p.genes$p_val<p.val)]
  p.genes<-p.genes[genes,,drop=FALSE]
  p.genes$gene<-genes
  
  if(bar){
    print(ggplot(data = p.genes[1:num.genes,,drop=FALSE])+xlab('gene')+ggtitle('Differential expression genes vs adj,p.value')+
            geom_bar(aes(x=reorder(gene,-p_val),y=p_val,fill=gene),stat = 'identity')+
            theme(axis.text.x = element_text(face="bold",size=10,hjust = 1,angle = 90),
                  legend.position = 'none'))
  }
  return(genes)
}

#  The genes select
#  but ITGB4 gene is not in the selected all genes which were selected using DESeq
human.p.genes<-DESeq_Genes(p.data = human.all.DESeq[[2]])




#  ITGB4 gene : plot
GenePlot(human.all.DESeq[[1]],gene1 = 'ITGB4',gene2 = human.p.genes[1])
VlnPlot(human.all.DESeq[[1]],features.plot = 'ITGB4',y.lab.rot = 90)  # Violinn plot of gene ITGB in all sample
VlnPlot(human.all.DESeq[[1]],features.plot = human.p.genes[1:6],y.lab.rot = 90)  # Violinn plot of gene ITGB in all sample


#  Perform linear dimensioal reduction (PCA),tsne

PCA.TSNE<-function(object,pcs.compute=FALSE,#pcs.plot=FALSE,tsne.plot=TRUE,
                   pc.genes=NULL,sig.genes=NULL,num.pcs=18){
  
  # object : Seurate object
  # pcs.compute : whether to compute the significant  pc components
  # pcs.plot: plot pcs(dimPlot),true or false
  # tsne.plot : plot tsne(dimplot),true or false
  
  count.data<<-object@raw.data
  scale.data<-as.data.frame(object@scale.data)
  
  if(pcs.compute){
    #   decide the number of significant pc
    y = sig.pcs.perm(dat=t(count.data), center=T, scale=T,
                     max.pc=100, B=1000, n.cores=20, randomized=T) 
    
    object<-RunPCA(object = object,pc.genes = pc.genes,pcs.compute = y$r)
    object<-ProjectPCA(object = object,do.print = FALSE)  #  Project on all genes in dataset
    
    object<-RunTSNE(object = object,reduction.use = 'pca',genes.use = sig.genes,
                      dims.use = 1:y$r,do.fast = FALSE)
    
    return(list(object,y))
  }else{
    object<-RunPCA(object = object,pc.genes = pc.genes,pcs.compute = num.pcs,do.print = FALSE)
    object<-ProjectPCA(object = object,do.print = FALSE)
    
    object<-RunTSNE(object = object,reduction.use = 'pca',genes.use = sig.genes,
                    dims.use = 1:num.pcs,do.fast = FALSE)
    
    return(object)
  }
  
  # if(pcs.plot){
  #   DimPlot(object = object,reduction.use = 'pca',pt.size = 4)
  #   DimHeatmap(object = object,reduction.type ='pca')
  #   #plot_grid(p1,p2)
  # }
  # 
  # if(tsne.plot){
  #   DimPlot(object =object,reduction.use = 'tsne',pt.size =3)
  #   
  #   #  Feature plot
  #   
  #   if(!'ITGB4'%in%rownames(count.data)){
  #     FeaturePlot(object = object,features.plot =object@var.gene[1:5])
  #   }
  #   
  #   if('ITGB4'%in%rownames(count.data)){
  #     FeaturePlot(object = object,features.plot =unique(c('ITGB',object@var.gene[1:5])))
  #   }
  # }
}

all.pbmc<-PCA.TSNE(object = human.all.DESeq[[1]],pcs.compute = FALSE)

#  Figure analysis

FeaturePlot(object = all.pbmc,features.plot ='ITGB4',pt.size = 4)  # ITGB4 gene in all dataset
FeaturePlot(object = all.pbmc,features.plot =human.p.genes[1:6],pt.size = 4)  # ITGB4 gene in all dataset


DimPlot(all.pbmc,reduction.use = 'tsne',pt.size = 4)  #  grour by sample
DimHeatmap(all.pbmc,reduction.type = 'pca')


FeatureHeatmap(all.pbmc,features.plot = 'ITGB4',pt.size = 3,plot.horiz = TRUE,
               cols.use = c("lightgrey", "blue"))

FeatureHeatmap(all.pbmc,features.plot =human.p.genes[1:6],pt.size = 3,plot.horiz = TRUE,
               cols.use = c("lightgrey", "blue"))


#  Heatmap for select genes
human.heatmap<-Heatmap_fun(genes = human.p.genes[1:20],tpm.data = human.all.DESeq[[1]]@raw.data,
                           condition = unique(as.character(human.all.DESeq[[1]]@ident)),
                           all.condition = as.character(human.all.DESeq[[1]]@ident))

NMF::aheatmap(human.heatmap[[2]], Rowv = NA, Colv = NA, annCol = human.heatmap[[1]], 
              scale = "row")






#   Condition control 
#  eg.  Tube_1 VS   Tube_2

#  Read condition data
human.condition.pro<-Load_data(data_dir = './human.txt',
                               condition.sample_dir = './human.condition.sample.txt')
#  Control condition
control.condition<-human.condition.pro[[2]]
control.sample.names<-colnames(human.condition.pro[[1]])

#  Control condition data frame(Tube_1  Tube_2)
human.part.DESeq<-DESeq_SeuratObj(X=human.condition.pro[[1]],cells.1 = 'Tube_1',cells.2 = 'Tube_2',
                                  condition = control.condition,condition.name = control.sample.names)


#   genes select under the control of Tube_1 Tube_2
human.p.part.genes<-DESeq_Genes(p.data = human.part.DESeq[[2]])

GenePlot(human.part.DESeq[[1]],gene1 = 'ITGB4',gene2 = human.p.part.genes[1])

#  violon plot
VlnPlot(human.part.DESeq[[1]],features.plot = 'ITGB4',y.lab.rot = 90)  # Violinn plot of gene ITGB in all sample
VlnPlot(human.part.DESeq[[1]],features.plot = human.p.part.genes[1:6],y.lab.rot = 90)  # Violinn plot of gene ITGB in all sample


# PCA  TSNE
part.pbmc<-PCA.TSNE(object = human.part.DESeq[[1]],pcs.compute = FALSE,num.pcs = 7)


#  Figure analysis
FeaturePlot(object = part.pbmc,features.plot ='ITGB4',pt.size = 4,no.legend = FALSE)  # ITGB4 gene in part dataset
FeaturePlot(object = part.pbmc,features.plot =human.p.part.genes[1:6],pt.size = 4,no.legend = FALSE)  # ITGB4 gene in part dataset


DimPlot(part.pbmc,reduction.use = 'tsne',pt.size = 4)  #  grour by sample
DimHeatmap(part.pbmc,reduction.type = 'pca')


FeatureHeatmap(part.pbmc,features.plot = 'ITGB4',pt.size = 3,plot.horiz = TRUE,
               cols.use = c("lightgrey", "blue"))
 
human.heatmap<-Heatmap_fun(genes = human.p.part.genes[1:20],tpm.data = human.part.DESeq[[1]]@raw.data,
                           condition = unique(as.character(human.part.DESeq[[1]]@ident)),
                           all.condition = as.character(human.part.DESeq[[1]]@ident))

NMF::aheatmap(human.heatmap[[2]], Rowv = NA, Colv = NA, annCol = human.heatmap[[1]], 
              scale = "row")








########################################################################
# another approach to handle gene count data and use the DESeq  test 
#to detect genes Differential expression.

DESeq_CT<-function(count.data,condition.1,condition.2=NULL,simple=TRUE){
  if(simple){
  count.data<-count.data[which(unlist(apply(count.data,1,sum)>0)),]
  }else{count.data<-count.data}
  
  # each genes count +1:for caculation
  count.data<-t(apply(count.data,1,function(x){
    if(all(x!=0)){
      return(x)
    }else{
      return(x+1)    
    }
  }))
  if(is.null(condition.2)){
     condition.1<-as.factor(condition.1)
     coldata<-data.frame(row.names =colnames(count.data),condition.1)
     xdds<-DESeqDataSetFromMatrix(countData = count.data,colData =coldata,
                               design = ~condition.1)
     xxs<-DESeq(xdds)
  }
  
  if(!is.null(condition.2)){
    condition.2<-as.factor(condition.2)
    coldata<-data.frame(row.names =colnames(count.data),condition.1,condition.2)
    xdds<-DESeqDataSetFromMatrix(countData = count.data,colData =coldata,
                                 design = ~condition.1+condition.2)
    xxs<-DESeq(xdds)
  }
  
  return(xxs)
}

condition.1<-unlist(lapply(colnames(human.only.pro),function(x)return(str_split(x,'_')[[1]][1])))
condition.2<-unlist(lapply(colnames(human.only.pro),function(x)return(str_split(x,'_')[[1]][2])))

xdds<-DESeq_CT(count.data = human.only.pro[1:100,],condition.1 = condition.1)
plotDispEsts(xdds)

DESeq_result<-function(object,condition,p.val=0.05){
  
  uniq.condition<-unique(condition)
  combn.con<-combn(uniq.condition,2)
  
  result<-list()
  for(i in 1:dim(combn.con)[2]){
    res<-results(object = object,contrast = c('condition.1',combn.con[,i][1],combn.con[,i][2]))
    index<-which(res$padj<p.val)
    # pair.g<-paste(combn.con[,i][1],'_vs_',combn.con[,i][2],sep = '')
    # result<-c(result,list(pair.g=index))
    genes<-rownames(res[index,])
    #pair.g<-paste(combn.con[,i][1],'_vs_',combn.con[,i][2],sep = '')
    gene.l<-list(genes)
    names(gene.l)<-paste(combn.con[,i][1],'_vs_',combn.con[,i][2],sep = '')
    result<-c(result,gene.l)
  #   }else{
  #     genes<-NULL
  #     gene.l<-list(genes)
  #     names(gene.l)<-paste(combn.con[,i][1],'_vs_',combn.con[,i][2],sep = '')
  #     result<-c(result,gene.l)
  #   }
   }
  return(result)
}

#  Venn diagram
r<-DESeq_result(xdds,condition = condition.1)
x<-as.vector(r)

grid.draw(venn.diagram(x[1:3],filename = NULL,fill= c("dodgerblue", "goldenrod1", "darkorange1")))

grid.draw(venn.diagram(x[1:4],filename = NULL,
                       fill= c("dodgerblue", "goldenrod1", "darkorange1", "darkorchid1")))

