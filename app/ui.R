library(shiny)
source('helper.R')

shinyUI(fluidPage(
  titlePanel("UI page for C1 single cell data"),
  sidebarLayout(
    
    sidebarPanel(
      helpText("Read the data from path",
             strong("Seurat package"), "to analysis "),
      h6('C1 single cell data!'),
      fileInput('DataFile',label = 'Upload the data from the selected path',multiple = TRUE),
      numericInput('min.cells',label = 'min.cells',value = 0,min = 0,max = 1000,step = 1),
      numericInput('min.genes',label = 'min.genes',value = 0,min = 0,max = 1000,step = 1),
      hr(),
      # base on the helper.R function,whether use the condition
      h6('Use the group criter'),
      selectInput('group',label = 'GroupType',choices = c('body','size','control'),selected = 'body'),
      hr(),
      helpText('If the group is ',strong('control'),'will use the below file to extract control message'),
      fileInput('ConditionFile',label = 'Condition data path',multiple = TRUE),
      hr(),
      hr(),
      
      h5('Figure explore'),
      selectInput('PlotMethod',label = 'Figure Explore',
                  choices = c('Bar',
                              'Violin',
                              'Geneplot',
                              'PCA',
                              'TSNE',
                              "FeaturePlot"),
                  selected = 'Violin'),
      
      
      hr(),
      h5('Genes the Figure plot'),
      selectInput('gene',label = 'gene',choices = c("ITGB4",
                                                    "A2M",
                                                    "AAAS",
                                                    "AARS",
                                                    "ABCA5",
                                                    "AASDHPPT",
                                                    "AAGAB",
                                                    "AAMDC"),selected = 'ITGB4'),
      hr(),
      h6('Paramenters size for Dimention Reduction plot'),
      numericInput('pt.size',label = 'pt.size',value = 4,min = 0,max = 8,step = 1),
      
      hr(),
      h6('DownLoad the plot'),
      radioButtons(inputId = "var", label = "Select the file type", choices = list("png","pdf"))
    ),
    mainPanel(
      plotOutput("plot"),
      downloadButton(outputId = "down", label = "Download the plot")
      )
  )
))

# server<-function(input,output){}
# shinyApp(ui,server = server)