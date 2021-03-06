library(RColorBrewer)
brewer16 = c(brewer.pal(9, "Set1"), brewer.pal(7, "Set2"))
brewer16[6] = "khaki2"
brewer16[8] = "lightskyblue2"

cubehelix1.16 = c('#000000', '#1B0F00', '#411704', '#681B20',
                  '#85214B', '#932D7E', '#9042AF', '#8160D2', '#6F83E3',
                  '#63A6E2', '#65C5D3', '#78DBC2', '#99E9B9', '#C1F0BF', '#E6F5D8', '#FFFFFF')

### Compute the group-wise mean of a dataset.
group.means <- function(counts, groups, fn=mean, use.data.table=F)
{
  counts <- aggregate(t(counts), by=list(groups), FUN=fn)
  rownames(counts) = counts$Group.1
  counts$Group.1 = NULL
  r = t(counts)
  return(r)
}

# Logging utility function
info <- function(text, ...)
{
  cat(sprintf(paste(Sys.time(),"INFO:", text,"\n")))
}

# Logging utility function
warn <- function(text, ...)
{
  cat(sprintf(paste(Sys.time(),"WARN:", text,"\n")))
}

### Compute TPM expression values from raw UMI counts
tpm <- function(counts, mult=10000,eps=1e-6)
{
  info("Running TPM normalisation")
  total.counts = colSums(counts)
  scaled.counts = t(t(counts) / total.counts)
  scaled.counts * mult
}


### Run ComBat batch correction from the SVA package
batch.normalise.comBat <- function(counts, batch.groups, max.val=6)
{
  batch.groups = factor(batch.groups) ## drop zero levels
  batch.id = 1:length(unique(batch.groups))
  names(batch.id) = unique(batch.groups)
  batch.ids = batch.id[batch.groups]
  correct.data = ComBat(counts,batch.ids, prior.plots=FALSE, par.prior=T)
  correct.data[correct.data > max.val] = max.val
  as.data.frame(correct.data)
}


### Get variable genes. Code adapted from:
### | Brennecke et al, Accounting for technical noise in single-cell RNA-seq experiments
### | Nature Methods 10, 1093–1095 (2013), doi:10.1038/nmeth.2645
###     See: https://images.nature.com/original/nature-assets/nmeth/journal/v10/n11/extref/nmeth.2645-S2.pdf
###     and: http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
get.variable.genes <- function(ed, min.cv2=2, pdf=NULL, width=9, height=8, do.plot=T, p.thresh=0.05)
{
  
  means <- rowMeans(ed)
  vars <- apply(ed,1,var)
  cv2 <- vars/means^2
  minMeanForFit <- unname( quantile( means[ which( cv2 > min.cv2 ) ], .95 ) )
  useForFit <- means >= minMeanForFit # & spikeins
  info(sprintf("Fitting only the %s genes with mean expression > %s", sum(useForFit), minMeanForFit))
  fit <- glm.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ), cv2[useForFit] )
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"])
  if(do.plot){par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(means),log(cv2))}
  xg <- exp(seq( min(log(means[means>0])), max(log(means)), length.out=1000 ))
  vfit <- a1/xg + a0
  if(do.plot){lines( log(xg), log(vfit), col="black", lwd=3 )}
  
  df <- ncol(ed) - 1
  # add confidence interval
  if(do.plot){
    lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black")
    lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black")
  }
  afit <- a1/means+a0
  varFitRatio <- vars/(afit*means^2)
  varorder <- order(varFitRatio, decreasing=T)
  oed <- ed[varorder,]
  pval <- pchisq(varFitRatio*df,df=df,lower.tail=F)
  adj.pval <- p.adjust(pval,"fdr")
  r = data.frame(rownames(ed), varFitRatio, pval, adj.pval)
  colnames(r) = c("Gene", "VarianceFitRatio", "p", "p.adj")
  v = r[!is.na(r$p.adj),]
  n.sig = sum(v$p.adj<p.thresh)
  info(sprintf("Found %s variable genes (p<0.05)", n.sig))
  
  # add top 100 genes
  if(do.plot){
    points(log(means[varorder[1:n.sig]]),log(cv2[varorder[1:n.sig]]),col=2)
  }
  r = r[order(r$VarianceFitRatio, decreasing=T), ]
  r$Rank = 1:nrow(r)
  return(r)
}


# Test for significant PCs adapted from:
#
#       ' Permutation Parallel Analysis
#       '
#       ' Estimate a number of significant principal components from a permutation test
#   B is the number of permutations
#   threshold is p-value for significance
#'
sig.pcs.perm <- function (dat, B = 100, threshold = 0.05,
                          randomized=F,
                          verbose=TRUE, seed = NULL,
                          max.pc=100, n.cores=1,
                          center=T, scale=T) {
  ptm <- proc.time()
  if(B %% n.cores != 0){stop("Permutations must be an integer multiple of n.cores")}
  cat(sprintf("Scaling input matrix [center=%s, scale=%s]\n", center, scale))
  dat = t(dat)
  dat = as.matrix(t(scale(t(dat), center=center, scale=scale)))
  if (!is.null(seed)) set.seed(seed)
  n <- min(max.pc, ncol(dat))
  m <- nrow(dat)
  print(paste0("Considering only the top ", n, " PCs. Supply max.pc if you wish to change"))
  cat(sprintf("Running initial PCA\n"))
  if(randomized){
    library(rsvd)
    uu <- rsvd(as.matrix(dat), k=max.pc)
  }else{
    uu <- corpcor::fast.svd(dat, tol = 0)
  }
  ndf <- n - 1
  dstat <- uu$d[1:ndf]^2/sum(uu$d[1:ndf]^2)
  dstat0 <- matrix(0, nrow = B, ncol = ndf)
  if(verbose==TRUE) message("Estimating number of significant principal components. Permutation: ")
  #permutations
  if(n.cores==1){
    for (i in 1:B) {
      if(verbose==TRUE) cat(paste(i," "))
      dat0 <- t(apply(dat, 1, sample, replace = FALSE))
      if(randomized){
        library(rsvd)
        uu0 <- rsvd(as.matrix(dat0), k=max.pc)
      }else{
        uu0 <- corpcor::fast.svd(dat0, tol = 0)
      }
      dstat0[i, ] <- uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
    }
  }else{
    library(parallel)
    library(foreach)
    library(doParallel)
    cl<-makePSOCKcluster(n.cores, outfile="")
    registerDoParallel(cl, n.cores)
    chunksize = B/n.cores
    vals = split(1:B, ceiling(seq_along(1:B)/chunksize))
    dstat0 = foreach(run.id=1:n.cores, .packages="corpcor", .combine=cbind) %dopar% {
      v = vals[[run.id]]
      #cat(sprintf("Core %s will run perms: %s \n", run.id, paste(v, collapse=",")))
      do.call(rbind, lapply(v, function(i) {
        if(verbose==TRUE) cat(paste(i," "))
        dat0 <- t(apply(dat, 1, sample, replace = FALSE))
        
        if(randomized){
          library(rsvd)
          uu0 <- rsvd(as.matrix(dat0), k=max.pc)
        }else{
          uu0 <- corpcor::fast.svd(dat0, tol = 0)
        }
        uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
      }))
    }
    cat("\nUnregistering parallel backend..")
    stopCluster(cl)
    registerDoSEQ()
    cat(" done\n");
  }
  p <- rep(1, n)
  for (i in 1:ndf) {
    p[i] <- mean(dstat0[, i] >= dstat[i])
  }
  for (i in 2:ndf) {
    p[i] <- max(p[(i - 1)], p[i])
  }
  r <- sum(p <= threshold)
  y = proc.time() - ptm
  cat(sprintf("\n\n PC permutation test completed. \n %s PCS significant (p<%s, %s bootstraps)\n Runtime: %s s\n ", r,  threshold, B,signif(y[["elapsed"]], 3)))
  return(list(r = r, p = p))
}

build_knn_graph <- function(dm, k=200, verbose=F)
{
  if(k==0)
  {
    k = floor(sqrt(nrow(dm))/2)
  }
  if(verbose)
  {
    info(sprintf("Building %s-nearest [%s] neighbor graph..", k, dist.type))
  }
  g <- nng(dx=dm,k=k)
  V(g)$name = rownames(dm)
  if(verbose)
  {
    info(sprintf("%s %s-NN computed. Average degree: %s", dist.type, k, mean(degree(g))))
  }
  return(g)
}






# graph.type can be jaccard, invlogweighted or dice, community detect
# can be louvain, infomap or markov.
cluster_graph <- function(      g,
                                graph.type="knn", # can be threshold (binarise the distance matrix), jaccard or knn.
                                dm=NULL,
                                community.detect="infomap",
                                distance.method="euclidean",
                                k=0)
{
  if(identical(toupper(community.detect), toupper("markov")))
  {
    r = igraph::cluster.markov(g)
    clusters = r$Cluster
  }else{
    if(identical(toupper(community.detect), toupper("louvain")))
    {
      r = igraph::multilevel.community(as.undirected(g))
      clusters = r$membership
    }else{
      if(identical(toupper(community.detect), toupper("infomap")))
      {
        r = igraph::infomap.community(g, modularity=TRUE)
        clusters = r$membership
      }else{
        error(sprintf("Unknown community detection method: %s", community.detect))
        return (FALSE)
      }
    }
  }
  n.clusters =length(unique(clusters))
  f = function(i){as.vector(clusters==i)}
  clist= lapply(1:n.clusters, f)
  m = igraph::modularity(g, clusters)
  return (list("result"=r,
               "clustermethod"=paste(graph.type, "-graph clustering [", community.detect,"]", sep=""),
               "nc"=n.clusters,
               "modularity"=m,
               "clusterlist"=clist,
               "partition"=clusters))
}


merge_clusters <- function(clustering, clusters.to.merge, new.name=NULL)
{
  if(length(clustering) < 2){cat("ERROR: Must provide 2 or more cluster ID's to merge!");return (clustering)}
  i = 1
  if(!is.null(new.name)){
    use.id = new.name
    levels(clustering) = c(levels(clustering), use.id)
    clustering[which(clustering == clusters.to.merge[1])] = use.id
  }else
  {use.id = clusters.to.merge[1]}
  for(id in clusters.to.merge)
  {
    if(i > 1)
    {
      cat(sprintf("Merging cluster %s into %s ..\n", id, use.id))
      clustering[which(clustering == id)] = use.id
    }
    i = i + 1
  }
  factor(clustering)
}


#  Create PCA,TNSE scores
PCA_TSNE.scores<-function(data.tpm,data.umis,var_genes,data_name,is.var.genes=TRUE,sig.pcs=TRUE){
  if(is.var.genes){
    X<-data.tpm[var_genes,]
    Y<-data.umis[var_genes,]
    
  }else{
    X<-data.tpm
    Y<-data.umis
  }
  
  pca_name<-paste(data_name,'_pca_scores.txt',sep='')
  if(file.exists(pca_name)){
    cat(sprintf('%s exists\n',pca_name))
    pca<-read.table(pca_name)
  }else{
    pca = rpca(t(X), center=T, scale=T, retx=T, k=100)$x
    write.table(pca,file = pca_name,quote = F)
  }
  
  
  tsne_name<-paste(data_name,'_tsne_scores.txt',sep='')
  if(file.exists(tsne_name)){
    cat(sprintf('%s exists\n',tsne_name))
    tsne.rot = read.table(tsne_name)
  }else{
    if(sig.pcs){
    y = sig.pcs.perm(dat=t(Y), center=T, scale=T,
                    max.pc=100, B=1000, n.cores=20, randomized=T)  
    barnes_hut_tsne = Rtsne(pca[, 1:y$r], check_duplicates=T,pca=FALSE, #dont run PCA again
                            initial_dims = y$r, perplexity = 20, max_iter = 20000, verbose=T, whiten=F)
    }else{
      barnes_hut_tsne = Rtsne(pca[, 1:20], check_duplicates=T,pca=FALSE, #dont run PCA again
                              initial_dims = 20, perplexity = 20, max_iter = 20000, verbose=T, whiten=F)
    }
    
    tsne.rot = barnes_hut_tsne$Y
    tsne.rot<-as.data.frame(tsne.rot)
    colnames(tsne.rot)<-c('tSNE_1','tSNE_2')
    write.table(tsne.rot,file = tsne_name,quote = F)
  }
  return(tsne.rot)
}

get.variable.genes_cvdiff <- function(ed,do.plot=T){
  
  means <- rowMeans(ed)
  vars <- apply(ed,1,var)
  cv <- sqrt(vars)/means
  #minMeanForFit <- unname( quantile( means[ which( cv > min.cv ) ], .95 ) )
  #useForFit <- means >= minMeanForFit # & spikeins
  #info(sprintf("Fitting only the %s genes with mean expression > %s", sum(useForFit), minMeanForFit))
  fit <- glm.fit( cbind( a0 = 1, a1tilde = 1/means),cv) #[useForFit] ), cv[useForFit] )
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"])
  if(do.plot){par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(means),log(cv))}
  xg <- exp(seq( min(log(means[means>0])), max(log(means)), length.out=1000 ))
  vfit <- a1/xg + a0
  if(do.plot){lines( log(xg), log(vfit), col="black", lwd=3 )}
  
  df <- ncol(ed) - 1
  # add confidence interval
  if(do.plot){
    lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black")
    lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black")
  }
  afit <- a1/means+a0  # fited value
  CVdiff<-cv-afit
  
  CVdiff_mean<-mean(CVdiff)
  CVdiff_std<-sqrt(var(CVdiff))
  delta<-CVdiff_mean+1.67*CVdiff_std
  
  r<-data.frame(rownames(ed),CVdiff)
  colnames(r)<-c('Genes','CVdiff')
  r<-r[r$CVdiff>delta,]
  return(r)
}


load_data<-function(data_name,web='ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92332/suppl/'){
	  data_dir<-paste(web,data_name,sep='')
  if(file.exists(data_name)){
	      print('File exists!')
      data<-read.delim(data_name)
        }else{
		    download.file(web,destfile=basename(data_name))
          data<-read.delim(data_name)
        }
	  info(sprintf("Data dimensions: %s" , paste(dim(data), collapse = "x")))
    return(data)
}

#  Create data for TPM (ggplot)
Genes_mean_tpm<-function(genes,tpm_data,tsne_data,title,fun=mean,doplot=TRUE){
  Log2TPM<-as.numeric(apply(tpm_data[genes,],2,fun,na.rm=TRUE))
  if(doplot){
    title_1<-paste(genes,collapse = ',')
    title_2<-paste(title,'(',title_1,')',sep='')
    print(ggplot(tsne_data, aes(x=tSNE_1, y=tSNE_2))+geom_point(aes(color=Log2TPM))+theme(legend.title = element_text(size=8,color='blue',face='bold'),
                                                                                    legend.position = 'right') +ggtitle(title_2)+
                                                                                     scale_color_gradient2(low='lightblue',mid='green',high='red',name='Log2\nTPM+1'))
  }
  else{
    return(Log2TPM)
  }
}




## melt data for heatmap(ggplot2)
Create_plot_data<-function(genes,fun=scale,origin.data,var_genes,cell_groups,
                           if.use.var.genes=FALSE,cells=NULL){
  
  #   origin.data: [genes,cells]
  #   var_genes:the genes to use from get.variable.genes  function;if use,if.use.var.genes=TRUE
  #   cell_groups:  extract from sample names
  #   cells:  whether need to select cells.when NULL,not select
  genes.names<-rownames(origin.data)
  for(g in genes){
    if(!g%in%genes.names){
      cat(sprintf('%s is not exists\n',g))
    }
  }
  
  if(if.use.var.genes){
    origin.data_1<-origin.data[var_genes,]
  }else{
    origin.data_1<-origin.data
  }
  
  origin.data_2<-as.data.frame(t(origin.data_1))
  
  if(!is.null(cells)){
    heatp_1<-origin.data_2[cell_groups%in%cells,]
    
    heatp_2<-heatp_1[,colnames(heatp_1)%in%genes]
    heatp_3<-as.data.frame(apply(heatp_2,2,fun))
    heatp_3[,'Groups']<-cell_groups[cell_groups%in%cells] 
  }
  if(is.null(cells)){
    heatp_2<-origin.data_2[,colnames(origin.data_2)%in%genes]
    heatp_3<-as.data.frame(apply(heatp_2,2,fun))
    heatp_3[,'Groups']<-cell_groups
  }
  info(sprintf("Data dimensions: %s" , paste(dim(heatp_3), collapse = "x")))
  
  heatp_4<-melt(heatp_3)
  heatp_4<-heatp_4[order(heatp_4$Groups),]
  
  
  return(heatp_4)
  
}


Facet_wrap_fun<-function(gene,tpm.data,tsne.data,
                         condition=c('Control','Salm'),all.condition=Salmonellalnfect.condition){
  cat(sprintf('There ara %d conditions\n',length(condition)))
  tsne<-data.frame()
  for(i in 1:length(condition)){
    tsne<-rbind(tsne,tsne.data[all.condition%in%condition[i],])
  }
  cat(sprintf('Whether creat data accurate %d \n',sum(dim(tsne.data)[1]==dim(tsne)[1])))
  
  ###  create  gene expression TPM data
  gene.mp<-c()
  for(i in 1:length(condition)){
    gene.mp<-c(gene.mp,as.numeric(tpm.data[gene,all.condition%in%condition[i]]))
  }
  
  ### create Condition
  Condition<-c()
  for(i in 1:length(condition)){
    Condition<-c(Condition,rep(condition[i],sum(all.condition%in%condition[i])))
  }
  tsne$Gene.Mp<-gene.mp
  tsne$Condition<-Condition
  tsne$Gene<-rep(gene,dim(tsne)[1])
  
  return(tsne)
}


Heatmap_fun<-function(genes,tpm.data,condition,all.condition){
  cat(sprintf('There ara %d conditions\n',length(condition)))
  tpm<-data.frame()
  for(i in 1:length(condition)){
    tpm<-rbind(tpm,t(tpm.data[genes,all.condition%in%condition[i],drop=FALSE]))
  }
  #cat(sprintf('Whether creat data accurate %d \n',sum(dim(tpm.data)[1]==dim(tpm)[1])))
  tpm<-data.frame(t(tpm))
  cat(sprintf('Whether creat data accurate %d \n',sum(dim(tpm.data)[1]==dim(tpm)[2])))
  ### create Condition
  Condition<-c()
  for(i in 1:length(condition)){
    Condition<-c(Condition,rep(condition[i],sum(all.condition%in%condition[i])))
  }
  
  return(list(Condition,tpm))
}


###   Find the k to make 3 cluster
###  k-nearest method
Find_K<-function(K,pca.data,n=3){
  dm<-as.matrix(dist(pca.data))
  for(k in K){
    knn<-build_knn_graph(dm,k=k)
    clustering<-cluster_graph(knn)$partition
    if(length(unique(clustering))==n){
      cat(sprintf('Find the K:%d\n',k))
      return(k)
      break
    }
  }
}


#  Read geneCount data functiuon
Merge_GeneCount<-function(path){
  
  ##  Create data frame X:X has the rownames and colaname we defind
  WeData<-function(df){
    x<-data.frame(df$V2)
    rownames(x)<-as.character(df$V1)
    return(x)
  }
  
  Sample_names<-list.files(path)
  Data_path<-paste(path,Sample_names,sep = '')
  
  df<-fread(Data_path[1])
  X<-WeData(df)
  colnames(X)<-str_split(Sample_names[1],'\\.')[[1]][1]
  
  for(i in 2:length(Data_path)){
    x<-fread(Data_path[i])
    x<-WeData(x)
    colnames(x)<-str_split(Sample_names[i],'\\.')[[1]][1]
    X<-cbind(X,x)
  }
  
  all.genes<-rownames(X)
  all.genes<-setdiff(all.genes,all.genes[str_detect(all.genes,'^[-_]{1,}')])
  X<-X[all.genes,]
  return(X)
}


##  The function for selecting genes simply before variable select
Simplify_Select_Genes<-function(data,method=NULL,seq.len=NULL,id=2,fun=mean){
  #  method:NULL,quantile
  if(is.null(method)){
    data<-data[which(unlist(apply(data,1,sum)>0)),]
    return(data)
  }
  
  if(method=='quantile'){
    #  fun:min,max,mean,sum
    each.data<-unlist(apply(data,1,fun))
    seq.data<-as.numeric(quantile(each.data,probs = seq(0,1,1/seq.len))[id])
    data<-data[which(each.data>seq.data),]
    return(data)
  }
}

#  Thr function :DESeq analysis
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


#  return the genes from DESeq test
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