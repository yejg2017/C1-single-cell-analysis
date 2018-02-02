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
mouse<-Merge_GeneCount(data_dir[3])

mouse.sample.names<-colnames(mouse)
#mouse.sample.names<-unlist(lapply(mouse.sample.names,function(x)return(str_replace_all(x,'-','_'))))
mouse.groups<-unlist(lapply(mouse.sample.names,function(x)return(str_split(x,"_")[[1]][1])))


#  consistent sample names
mouse.sample.names<-unlist(lapply(mouse.sample.names,function(x){
  x<-str_replace_all(x,'-','_')
  v<-str_split(x,'_')[[1]]
  if(str_detect(v[2],'^[0-9]')){
    return(x)
  }else{
    v<-str_split(x,'_',2)[[1]]
    z<-str_split(v[1],'')[[1]]
    if(length(z)>=9){
      k<-str_c(str_c(str_c(z[1:5],collapse =''),'_',str_c(str_c(z[6:7],collapse = ''),'um',
                                                          collapse = '')),str_c(z[8:length(z)],collapse = ''),v[2],sep = '_')
      return(k)
      }
  }
}))

colnames(mouse)<-mouse.sample.names
mouse.groups<-unlist(lapply(mouse.sample.names,function(x){
  return(str_split(x,'_')[[1]][1])
}))

cell.size<-unlist(lapply(mouse.sample.names,function(x){
  return(str_split(x,'_')[[1]][2])
}))

mouse.groups.size<-unlist(lapply(mouse.sample.names,function(x){
  v<-str_split(x,'_')[[1]]
  return(str_c(v[1],'_',v[2]))
}))

#     Other message: condition

