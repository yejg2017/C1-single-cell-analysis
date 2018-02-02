library(data.table)
library(stringr)
source('../tools.R')

data_dir<-'../../C1data/Macaca/'
all.files<-list.files(data_dir)
Count.log<-c()  #log file
geneCounts.file<-c()  # genecount file
for(f in all.files){
  if(str_detect(f,".geneCounts")){
    geneCounts.file<-c(f,geneCounts.file)
  }else{
    Count.log<-c(f,Count.log)
  }
}

Combind_Raw<-function(files){
  
  ##  Create data frame X:X has the rownames and colaname we defind
  WeData<-function(df){
    x<-data.frame(df$V2)
    rownames(x)<-as.character(df$V1)
    return(x)
  }
  
  
  df<-fread(paste(data_dir,files[1],sep=''))
  X<-WeData(df)
  colnames(X)<-str_split(files[1],'\\.')[[1]][1]
  
  for(i in 2:length(files)){
    x<-fread(paste(data_dir,files[i],sep=''))
    x<-WeData(x)
    colnames(x)<-str_split(files[i],'\\.')[[1]][1]
    X<-cbind(X,x)
  }
  
  all.genes<-rownames(X)
  all.genes<-setdiff(all.genes,all.genes[str_detect(all.genes,'^[-_]{1,}')])
  X<-X[all.genes,]
  return(X)
}

Monkey<-Combind_Raw(geneCounts.file)
colnames(Monkey)<-unlist(lapply(colnames(Monkey),function(x)return(str_replace_all(x,'-','_'))))
rownames(Monkey)<-unlist(lapply(rownames(Monkey),function(x)return(str_to_upper(x))))


#  Check whether the important genes is in the Monkey data
important.genes<-c('ITGB4','KRT19','ACTB','KRT12','KRT5',
                   'GAPDH','KRT3','PAX6','WNT7A','KRT14','TP63','KRT10')


#  ABCB5 not found
important.genes.ano<-c('ENSMFAG00000001695','ENSMFAG00000003391','ENSMFAG00000031981',
                       'ENSMFAG00000003650','ENSMFAG00000038566','ENSMFAG00000001880',
                       'ENSMFAG00000001094','ENSMFAG00000019311','ENSMFAG00000044809',
                       'ENSMFAG00000035842','ENSMFAG00000035702','ENSMFAG00000038850')


all.Monkey.genes<-rownames(Monkey)
for(i in 1:length(important.genes)){
  idx<-which(all.Monkey.genes==important.genes.ano[i])
  all.Monkey.genes[idx]<-important.genes[i]
}
rownames(Monkey)<-all.Monkey.genes
#write.table(Monkey,file='../data/monkey.txt',quote = FALSE)