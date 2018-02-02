library(shiny)
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
library(DESeq2)
source('helper.R')

options(shiny.maxRequestSize=500*1024^2)  # set the max size allowed to upload
shinyServer(
  function(input,output)({
    # read the data
    # file1<-input$DataFile
    # file2<-input$ConditionFile
    # if(input$condition){
    #   if(is.null(file2)){stop('Condition file path can not be NULL')}
    #   data<-reactive({
    #       if(is.null(file1)){return()} 
    #       Load_data(data_dir =file1$datapath,condition.sample_dir =file2$datapath)
    # })
    # }else{
    #   data<-reactive({
    #     if(is.null(file1)){return()} 
    #     Load_data(data_dir =file1$datapath)
    #   })
    # }
    
    
    
    data<-reactive({
      file1<-input$DataFile
      file2<-input$ConditionFile
      if(is.null(file1)){
        return()
      }else{
          Load_data(data_dir =file1$datapath)
      }
    })
    
    object.1<-reactive({
      if(is.null(data())){return()}
      DESeq_SeuratObj(X = data(),min.cells = input$min.cells,min.genes = input$min.genes,
                      DESq = FALSE)
    })
    
    
    object.2<-reactive({
      if(is.null(object.1())){return()}
      if(input$group=='body'){
        object.1()
      }
      if(input$group=='size'){
        SetIdent(object.1(),cells.use =object.1()@cell.names,ident.use =
                   unlist(lapply(object.1()@cell.names,function(x)return(str_split(x,'_')[[1]][2]))))
      }
    })
    
    object.3<-reactive({
      if(is.null(object.2())){return()}
      PCA.TSNE(object = object.2(),pcs.compute = FALSE)
    })
    
    output$plot<-renderPlot({
      if(is.null(object.3())){return()}
      if(input$PlotMethod=='Bar'){
        group<-as.character(object.3()@ident)
        Group_Bar(count.data = object.3()@raw.data,group = group)
      }
      if(input$PlotMethod=='Violin'){
        VlnPlot(object = object.3(),features.plot =input$gene)
      }
      
      if(input$PlotMethod=='Geneplot'){
        GenePlot(object.3(),gene1 = 'ITGB4',gene2 =input$gene)
      }
      
      if(input$PlotMethod=='PCA'){
        DimPlot(object.3(),reduction.use = 'pca',pt.size = input$pt.size)
      }
      
      if(input$PlotMethod=='TSNE'){
        DimPlot(object.3(),reduction.use = 'tsne',pt.size = input$pt.size)
      }
      
      if(input$PlotMethod=='FeaturePlot'){
        FeaturePlot(object.3(),features.plot = input$gene,pt.size = input$pt.size,no.legend = FALSE)
      }
    })
    
    output$down<-downloadHandler(
      filename =function(){paste('C1_single_cell',input$var,sep = '.')},
      content = function(file){
        if(input$var=='png'){
          png(file)
        }else{
          pdf(file)
        }
        
        if(is.null(object.3())){return()}
        if(input$PlotMethod=='Bar'){
          group<-as.character(object.3()@ident)
          Group_Bar(count.data = data(),group = group)
        }
        if(input$PlotMethod=='Violin'){
          VlnPlot(object = object.3(),features.plot =input$gene)
        }
        
        if(input$PlotMethod=='Geneplot'){
          GenePlot(object.3(),gene1 = 'ITGB4',gene2 =input$gene)
        }
        
        if(input$PlotMethod=='PCA'){
          DimPlot(object.3(),reduction.use = 'pca',pt.size = input$pt.size)
        }
        
        if(input$PlotMethod=='TSNE'){
          DimPlot(object.3(),reduction.use = 'tsne',pt.size = input$pt.size)
        }
        
        if(input$PlotMethod=='FeaturePlot'){
          FeaturePlot(object.3(),features.plot =input$gene,pt.size = input$pt.size,no.legend = FALSE)
        }
        dev.off()
      })
      }
)
)