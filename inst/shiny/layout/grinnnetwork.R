mainPanel(width=12,
fluidRow(column(12,
  mainPanel(width=12,
    h3("fetchGrinnNetwork"),
    p("Query a network (grinn network) using information from the internal graph database, see ",
      a(href='http://kwanjeeraw.github.io/grinn/fetchgrinn.html',target='_blank','here'),' for argument details.'
    )
  )#end mainPanel
)),
wellPanel(
  h4("Input arguments:"), helpText("* required field"),
  fluidRow(
    column(4, fileInput(inputId='txtInput', label='txtInput *', accept=c('text/csv','text/comma-separated-values,text/plain','.csv'))),
    column(8, mainPanel(tableOutput('txtExTable')))
  ),
  hr(),
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