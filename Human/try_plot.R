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
source('../tools.R')
library(DESeq2)
############
# QC data: how to remove bad data,genes
human.only.pro<-Load_data(data_dir = '../data/human.txt')
important.genes<-c('ITGB4','ABCB5','KRT19','ACTB','KRT12','KRT5',
                   'GAPDH','KRT3','PAX6','WNT7A','KRT14','TP63','KRT10')

table(unlist(lapply(colnames(human.only.pro),
                    function(x)return(str_split(x,'_')[[1]][2]))))
table(unlist(lapply(colnames(human.only.pro),
                    function(x)return(str_split(x,'_')[[1]][1]))))

human.only.pro<-human.only.pro[,colnames(human.only.pro)[unlist(lapply(colnames(human.only.pro),
                                                                       function(x)return(str_split(x,'_')[[1]][2])))%in%c('10um','20um','6um')]]

human.only.pro<-human.only.pro[,colnames(human.only.pro)[!unlist(lapply(colnames(human.only.pro),
                                        function(x)return(str_split(x,'_')[[1]][1])))%in%c('hc001','shoutiao')]]

sample.genes.counts<-unlist(apply(human.only.pro,2,sum))



#############

ITGB4<-as.integer(human.only.pro['ITGB4',])
Positive.idx<-which(ITGB4>0)
Negative.idx<-which(ITGB4==0)
Positive.data<-human.only.pro[,Positive.idx,drop=FALSE]
Negative.data<-human.only.pro[,Negative.idx,drop=FALSE]


Positive.data<-data.frame(t(Positive.data[important.genes,]))
Negative.data<-data.frame(t(Negative.data[important.genes,]))
plot.data<-rbind(Positive.data,Negative.data)
plot.data$condition<-rep(c('Positive','Negative'),times=c(dim(Positive.data)[1],dim(Negative.data)[1]))
cell.size<-c(unlist(lapply(rownames(Positive.data),function(x) return(str_split(x,'_')[[1]][2]))),
             unlist(lapply(rownames(Negative.data),function(x) return(str_split(x,'_')[[1]][2]))))


plot.data$cell.size<-cell.size

library(ggplot2)
library(reshape2)

X<-melt(plot.data)
p<-ggplot(data = X[X$variable=="KRT19",],aes(y=value,x=as.factor(condition),fill=as.factor(cell.size)))
p+geom_violin(scale = "width")+geom_jitter()+guides(fill=guide_legend(title="Cell Size"))
p+geom_boxplot()+guides(fill=guide_legend(title="Cell Size"))

#  一起分层
ggplot(data = X,aes(x=value,fill=variable))+
  geom_density(kernel='gaussian')+scale_x_log10()+facet_wrap(~condition+cell.size)


ggplot(data = X,aes(x=value,fill=variable))+
  geom_density(kernel='gaussian',position = 'stack')+scale_x_log10()+facet_wrap(~condition+cell.size)


ggplot(data = X,aes(x=value,fill=variable))+
   geom_histogram()+scale_x_log10()+facet_wrap(~condition+cell.size)


############
#  各自画density,histogram
Positive.data$cell.size<-unlist(lapply(rownames(Positive.data),function(x) return(str_split(x,'_')[[1]][2])))
Negative.data$cell.size<-unlist(lapply(rownames(Negative.data),function(x) return(str_split(x,'_')[[1]][2])))

Positive.melt<-melt(Positive.data)
Negative.melt<-melt(Negative.data)

#  Positive
Pp<-ggplot(data=Positive.melt,aes(value,fill=variable))
Pp+geom_density(kernel='gaussian',position = 'stack')+scale_x_log10()+facet_wrap(~cell.size)
Pp+geom_density(kernel='gaussian')+scale_x_log10()+facet_wrap(~cell.size)

# Negative
Pn<-ggplot(data=Negative.melt,aes(value,fill=variable))
Pn+geom_density(kernel='gaussian',position = 'stack')+scale_x_log10()+facet_wrap(~cell.size)
Pn+geom_density(kernel='gaussian')+scale_x_log10()+facet_wrap(~cell.size)






#########
library(easyGgplot2)
plot.1<-ggplot2.violinplot(data=X[X$variable=="KRT19",],xName = 'condition',yName = 'value',
                   groupName = 'cell.size',addMean = TRUE)

plot.1<-ggplot2.customize(plot.1,mainTitle="KRT19",xtitle='Condition')
print(plot.1)