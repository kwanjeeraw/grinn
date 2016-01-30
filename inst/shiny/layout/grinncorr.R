mainPanel(width=12,
fluidRow(column(12,
  mainPanel(width=12,
    h3("fetchGrinnCorrNetwork"),
    p("Reconstruct a grinn network using information from the internal graph database, then compute and combine with a weighted correlation network, see ",
      a(href='http://kwanjeeraw.github.io/grinn/fetchgrinncorr.html',target='_blank','here'),' for argument details.'
    )
  )#end mainPanel
)),
wellPanel(
  h4("Input arguments:"), helpText("* required field"),
  h5("Network query:"),
  fluidRow(
    column(4, fileInput(inputId='txtInput', label='txtInput *', accept=c('text/csv','text/comma-separated-values,text/plain','.csv'))),
    column(8, mainPanel(tableOutput('txtExTable')))
  ),
  fluidRow(
    column(6,radioButtons("from", label = "from",inline = TRUE, 
                          choices = list("metabolite"="metabolite","protein"="protein","gene"="gene","pathway"="pathway"),selected = "metabolite")),
    column(6,radioButtons("to", label = "to",inline = TRUE, 
                          choices = list("metabolite"="metabolite","protein"="protein","gene"="gene","pathway"="pathway"),selected = "pathway"))
  ),
  fluidRow(
    column(6,checkboxGroupInput("filterSource", label ="filterSource",inline = TRUE,
                                choices = list("SMPDB"="SMPDB","KEGG"="KEGG","REACTOME"="REACTOME"),selected = c("SMPDB","KEGG","REACTOME")))
  ),
  fluidRow(
    column(12,radioButtons("dbXref", label ="dbXref",inline = TRUE,
                           choices = list("grinn"="grinn","chebi"="chebi","kegg"="kegg","pubchem"="pubchem",
                                          "inchi"="inchi","hmdb"="hmdb","smpdb"="smpdb","reactome"="reactome",
                                          "uniprot"="uniprot","ensembl"="ensembl","entrezgene"="entrezgene"),selected = "grinn"))
  ),
  hr(),
  h5("Correlation analysis:"),
  fluidRow(
    column(8, radioButtons('sep', 'Delimiter',c(Comma=',',Tab='\t',Semicolon=';'),',',inline=TRUE))
  ),
  fluidRow(
    column(4, fileInput(inputId='datXInput', label='datX *', accept=c('text/csv','text/comma-separated-values,text/plain','.csv'))),
    column(8, mainPanel(tableOutput('datXExTable')))
  ),
  fluidRow(
    column(4, fileInput(inputId='datYInput', label='datY', accept=c('text/csv','text/comma-separated-values,text/plain','.csv'))),
    column(8, mainPanel(tableOutput('datYExTable')))
  ),
  fluidRow(  
    column(6,numericInput("corrCoef", label = "corrCoef *",value = 0.5, min = 0, max = 1, step = 0.1)),
    column(6,numericInput("pval", label = "pval *",value = 0.05, min = 0, max = 0.05, step = 0.001))
  ),
  fluidRow(
    column(6,radioButtons("method", label ="method",inline = TRUE,
                          choices = list("spearman"="spearman","pearson"="pearson","kendall"="kendall"),selected = "spearman"))
  ),
  hr(),
  actionButton("submit","Submit")
),
fluidRow(column(12,
  h4("Output network"),
  sidebarLayout(
    sidebarPanel(width=3,
      downloadButton('downloadEdge', 'Download Edge'),br(),
      downloadButton('downloadNode', 'Download Node')
    ),
    mainPanel(width=9,    
      tabsetPanel(type = "tabs", 
        tabPanel("Summary", verbatimTextOutput("summaryCode")),
        tabPanel("Nodes", tableOutput("nodeTable")),
        tabPanel("Edges", tableOutput("edgeTable"))
      )
    )
  )#end sidebarLayout
))
)