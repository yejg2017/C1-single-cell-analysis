library(stringr)
source('tools.R')

genecount=read.csv("a-9_geneCounts.csv",header = TRUE)
genecount=genecount[-(1:5),]

uniq.gene=unique(as.character(genecount$gene))
genecount=genecount[genecount$gene%in%uniq.gene,]
#genecount=genecount[!duplicated(genecount[,-1]),]
#rownames(genecount)=as.character(genecount$gene)


#reads=scale(genecount[,-1])
scale.reads=as.data.frame(scale(genecount[,-1]))
tpm.reads=data.frame(log(1+tpm(genecount[,-1])))

scale.reads<-cbind(data.frame(gene=genecount$gene),scale.reads)
tpm.reads<-cbind(data.frame(gene=genecount$gene),tpm.reads)
write.table(scale.reads,file = "scale.reads.csv",row.names = FALSE,sep=",")
write.table(tpm.reads,file = "tpm.reads.csv",row.names = FALSE,sep=",")
