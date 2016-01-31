mainPanel(width=12,
fluidRow(column(12,
  mainPanel(width=12,
    h3("convertToGrinnID"),
    p("Convert other database IDs to grinn IDs, see ",
      a(href='http://kwanjeeraw.github.io/grinn/convert.html',target='_blank','here'),' for argument details.'
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
    column(6,radioButtons("from", label = "nodetype",inline = TRUE, 
                          choices = list("metabolite"="metabolite","protein"="protein","gene"="gene","pathway"="pathway"),selected = "metabolite"))
  ),
  fluidRow(
    column(12,radioButtons("dbXref", label ="dbXref",inline = TRUE,
                           choices = list("chebi"="chebi","kegg"="kegg","pubchem"="pubchem",
                                          "inchi"="inchi","hmdb"="hmdb","smpdb"="smpdb","reactome"="reactome",
                                          "uniprot"="uniprot","ensembl"="ensembl","entrezgene"="entrezgene"),selected = "kegg"))
  ),
  hr(),
  actionButton("submit","Submit")
),
fluidRow(column(12,
  h4("Output"),
  sidebarLayout(
    sidebarPanel(width=3,
      downloadButton('downloadId', 'Download')
    ),
    mainPanel(width=9,    
      tabsetPanel(type = "tabs", 
        tabPanel("Summary", verbatimTextOutput("summaryCode")),
        tabPanel("Table", tableOutput("idTable"))
      )
    )
  )#end sidebarLayout
))
)