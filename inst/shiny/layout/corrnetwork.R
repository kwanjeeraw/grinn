mainPanel(width=12,
fluidRow(column(12,
  mainPanel(width=12,
    h3("fetchCorrNetwork"),
    p("Compute a weighted correlation network, see ",
      a(href='http://kwanjeeraw.github.io/grinn/fetchcorr.html',target='_blank','here'),' for argument details.'
    )
  )#end mainPanel
)),
wellPanel(
  h4("Input arguments:"), helpText("* required field"),
  fluidRow(
    column(8, radioButtons('sep', 'Delimiter',c(Comma=',',Tab='\t',Semicolon=';'),',',inline=TRUE))
  ),
  fluidRow(
    column(4, fileInput(inputId='datXInput', label='datNormX *', accept=c('text/csv','text/comma-separated-values,text/plain','.csv'))),
    column(8, mainPanel(tableOutput('datXExTable')))
  ),
  fluidRow(
    column(4, fileInput(inputId='datYInput', label='datNormY', accept=c('text/csv','text/comma-separated-values,text/plain','.csv'))),
    column(8, mainPanel(tableOutput('datYExTable')))
  ),
  hr(),
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