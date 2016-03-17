library(shiny)
library(DT)
library(grinn)
library(shinydashboard)
dashboardPage(skin = "black",
              dashboardHeader(title=img(src="logo.png", width = 200)),
              dashboardSidebar(
                sidebarMenu(
                  menuItem("fetchGrinnNetwork", tabName = "fgrinn", icon = icon("circle")),
                  menuItem("fetchCorrGrinnNetwork", tabName = "fcorrgrinn", icon = icon("circle")),
                  menuItem("fetchDiffCorrGrinnNetwork",  tabName = "fdiffgrinn", icon = icon("circle")),
                  menuItem("fetchGrinnCorrNetwork",  tabName = "fgrinncorr", icon = icon("circle")),
                  menuItem("fetchGrinnDiffCorrNetwork",  tabName = "fgrinndiff", icon = icon("circle")),
                  menuItem("convertToGrinnID",  tabName = "convertid", icon = icon("exchange")),
                  menuItem("setGrinnDb", tabName = "setdb", icon = icon("cogs"))
                )
              ),
              dashboardBody(
                fluidRow(
                  column(width = 12,
                         absolutePanel(
                           tags$head(
                             tags$link(rel = "stylesheet", type = "text/css", href = "grinn.css")
                           )
                         )
                  )
                ),
                tabItems(
                  tabItem(tabName = "fgrinn",
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
                            actionButton("query","Submit")
                          ),
                          fluidRow(column(12,
                                          h4("Output network"),
                                          sidebarLayout(
                                            sidebarPanel(width=3,
                                                         downloadButton('downloadEdgeBio', 'Download Edge'),br(),
                                                         downloadButton('downloadNodeBio', 'Download Node')
                                            ),
                                            mainPanel(width=9,    
                                                      tabsetPanel(type = "tabs", 
                                                                  tabPanel("Summary", verbatimTextOutput("consolemsgbio")),
                                                                  tabPanel("Nodes", dataTableOutput("nodeTableBio")),
                                                                  tabPanel("Edges", dataTableOutput("edgeTableBio"))
                                                      )
                                            )
                                          )#end sidebarLayout
                          ))
                  ),
                  tabItem(tabName = "fcorrgrinn",
                          fluidRow(column(12,
                                          mainPanel(width=12,
                                                    h3("fetchCorrGrinnNetwork"),
                                                    p("Compute a weighted correlation network and expand the network with information from the internal graph database, see ",
                                                      a(href='http://kwanjeeraw.github.io/grinn/fetchcorrgrinn.html',target='_blank','here'),' for argument details.'
                                                    )
                                          )#end mainPanel
                          )),
                          wellPanel(
                            h4("Input arguments:"), helpText("* required field"),
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
                            hr(),
                            fluidRow(  
                              column(6,numericInput("corrCoef", label = "corrCoef *",value = 0.5, min = 0, max = 1, step = 0.1)),
                              column(6,numericInput("pval", label = "pval *",value = 0.05, min = 0, max = 0.05, step = 0.001))
                            ),
                            fluidRow(
                              column(6,radioButtons("method", label ="method",inline = TRUE,
                                                    choices = list("spearman"="spearman","pearson"="pearson","kendall"="kendall"),selected = "spearman"))
                            ),
                            fluidRow(
                              column(6,radioButtons("xTo", label = "xTo",inline = TRUE, 
                                                    choices = list("non"="non","metabolite"="metabolite","protein"="protein","gene"="gene","pathway"="pathway"),selected = "metabolite")),
                              column(6,radioButtons("yTo", label = "yTo",inline = TRUE, 
                                                    choices = list("non"="non","metabolite"="metabolite","protein"="protein","gene"="gene","pathway"="pathway"),selected = "non"))
                            ),
                            fluidRow(
                              column(6,checkboxGroupInput("filterSourceCorr", label ="filterSource",inline = TRUE,
                                                          choices = list("SMPDB"="SMPDB","KEGG"="KEGG","REACTOME"="REACTOME"),selected = c("SMPDB","KEGG","REACTOME")))
                            ),
                            hr(),
                            actionButton("computeCorr","Submit")
                          ),
                          fluidRow(column(12,
                                          h4("Output network"),
                                          sidebarLayout(
                                            sidebarPanel(width=3,
                                                         downloadButton('downloadEdgeCorr', 'Download Edge'),br(),
                                                         downloadButton('downloadNodeCorr', 'Download Node')
                                            ),
                                            mainPanel(width=9,    
                                                      tabsetPanel(type = "tabs", 
                                                                  tabPanel("Summary", verbatimTextOutput("consolemsgcorr")),
                                                                  tabPanel("Nodes", dataTableOutput("nodeTableCorr")),
                                                                  tabPanel("Edges", dataTableOutput("edgeTableCorr"))
                                                      )
                                            )
                                          )#end sidebarLayout
                          ))
                  ),
                  tabItem(tabName = "fdiffgrinn",
                          fluidRow(column(12,
                                          mainPanel(width=12,
                                                    h3("fetchDiffCorrGrinnNetwork"),
                                                    p("Compute a differential correlation network and expand the network with information from the internal graph database, see ",
                                                      a(href='http://kwanjeeraw.github.io/grinn/fetchdiffcorrgrinn.html',target='_blank','here'),' for argument details.'
                                                    )
                                          )#end mainPanel
                          )),
                          wellPanel(
                            h4("Input arguments:"), helpText("* required field"),
                            fluidRow(
                              column(8, radioButtons('sepDiff', 'Delimiter',c(Comma=',',Tab='\t',Semicolon=';'),',',inline=TRUE))
                            ),
                            fluidRow(
                              column(4, fileInput(inputId='datXInputDiff', label='datX1 *', accept=c('text/csv','text/comma-separated-values,text/plain','.csv'))),
                              column(8, mainPanel(tableOutput('datXExTableDiff')))
                            ),
                            fluidRow(
                              column(4, fileInput(inputId='datX2Input', label='datX2 *', accept=c('text/csv','text/comma-separated-values,text/plain','.csv'))),
                              column(8, mainPanel(tableOutput('datX2ExTable')))
                            ),
                            hr(),
                            fluidRow(  
                              column(6,numericInput("pvalDiff", label = "pDiff *",value = 0.05, min = 0, max = 0.05, step = 0.001))
                            ),
                            fluidRow(
                              column(6,radioButtons("methodDiff", label ="method",inline = TRUE,
                                                    choices = list("spearman"="spearman","pearson"="pearson","kendall"="kendall"),selected = "spearman"))
                            ),
                            fluidRow(
                              column(6,radioButtons("xToDiff", label = "xTo",inline = TRUE, 
                                                    choices = list("non"="non","metabolite"="metabolite","protein"="protein","gene"="gene","pathway"="pathway"),selected = "metabolite"))
                            ),
                            fluidRow(
                              column(6,checkboxGroupInput("filterSourceDiff", label ="filterSource",inline = TRUE,
                                                          choices = list("SMPDB"="SMPDB","KEGG"="KEGG","REACTOME"="REACTOME"),selected = c("SMPDB","KEGG","REACTOME")))
                            ),
                            hr(),
                            actionButton("computeDiff","Submit")
                          ),
                          fluidRow(column(12,
                                          h4("Output network"),
                                          sidebarLayout(
                                            sidebarPanel(width=3,
                                                         downloadButton('downloadEdgeDiff', 'Download Edge'),br(),
                                                         downloadButton('downloadNodeDiff', 'Download Node')
                                            ),
                                            mainPanel(width=9,    
                                                      tabsetPanel(type = "tabs", 
                                                                  tabPanel("Summary", verbatimTextOutput("consolemsgdiff")),
                                                                  tabPanel("Nodes", dataTableOutput("nodeTableDiff")),
                                                                  tabPanel("Edges", dataTableOutput("edgeTableDiff"))
                                                      )
                                            )
                                          )#end sidebarLayout
                          ))
                  ),
                  tabItem(tabName = "fgrinncorr",
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
                              column(4, fileInput(inputId='txtInputCorr', label='txtInput *', accept=c('text/csv','text/comma-separated-values,text/plain','.csv'))),
                              column(8, mainPanel(tableOutput('txtExTableCorr')))
                            ),
                            fluidRow(
                              column(6,radioButtons("fromCorr", label = "from",inline = TRUE, 
                                                    choices = list("metabolite"="metabolite","protein"="protein","gene"="gene","pathway"="pathway"),selected = "metabolite")),
                              column(6,radioButtons("toCorr", label = "to",inline = TRUE, 
                                                    choices = list("metabolite"="metabolite","protein"="protein","gene"="gene","pathway"="pathway"),selected = "pathway"))
                            ),
                            fluidRow(
                              column(6,checkboxGroupInput("filterSourceCorr2", label ="filterSource",inline = TRUE,
                                                          choices = list("SMPDB"="SMPDB","KEGG"="KEGG","REACTOME"="REACTOME"),selected = c("SMPDB","KEGG","REACTOME")))
                            ),
                            fluidRow(
                              column(12,radioButtons("dbXrefCorr", label ="dbXref",inline = TRUE,
                                                     choices = list("grinn"="grinn","chebi"="chebi","kegg"="kegg","pubchem"="pubchem",
                                                                    "inchi"="inchi","hmdb"="hmdb","smpdb"="smpdb","reactome"="reactome",
                                                                    "uniprot"="uniprot","ensembl"="ensembl","entrezgene"="entrezgene"),selected = "grinn"))
                            ),
                            hr(),
                            h5("Correlation analysis:"),
                            fluidRow(
                              column(8, radioButtons('sepCorr', 'Delimiter',c(Comma=',',Tab='\t',Semicolon=';'),',',inline=TRUE))
                            ),
                            fluidRow(
                              column(4, fileInput(inputId='datXInputCorr', label='datX *', accept=c('text/csv','text/comma-separated-values,text/plain','.csv'))),
                              column(8, mainPanel(tableOutput('datXExTableCorr')))
                            ),
                            fluidRow(
                              column(4, fileInput(inputId='datYInputCorr', label='datY', accept=c('text/csv','text/comma-separated-values,text/plain','.csv'))),
                              column(8, mainPanel(tableOutput('datYExTableCorr')))
                            ),
                            fluidRow(  
                              column(6,numericInput("corrCoefCorr", label = "corrCoef *",value = 0.5, min = 0, max = 1, step = 0.1)),
                              column(6,numericInput("pvalCorr", label = "pval *",value = 0.05, min = 0, max = 0.05, step = 0.001))
                            ),
                            fluidRow(
                              column(6,radioButtons("methodCorr", label ="method",inline = TRUE,
                                                    choices = list("spearman"="spearman","pearson"="pearson","kendall"="kendall"),selected = "spearman"))
                            ),
                            hr(),
                            actionButton("computeCorr2","Submit")
                          ),
                          fluidRow(column(12,
                                          h4("Output network"),
                                          sidebarLayout(
                                            sidebarPanel(width=3,
                                                         downloadButton('downloadEdgeCorr2', 'Download Edge'),br(),
                                                         downloadButton('downloadNodeCorr2', 'Download Node')
                                            ),
                                            mainPanel(width=9,    
                                                      tabsetPanel(type = "tabs", 
                                                                  tabPanel("Summary", verbatimTextOutput("consolemsgcorr2")),
                                                                  tabPanel("Nodes", dataTableOutput("nodeTableCorr2")),
                                                                  tabPanel("Edges", dataTableOutput("edgeTableCorr2"))
                                                      )
                                            )
                                          )#end sidebarLayout
                          ))
                  ),
                  tabItem(tabName = "fgrinndiff",
                          fluidRow(column(12,
                                          mainPanel(width=12,
                                                    h3("fetchGrinnDiffCorrNetwork"),
                                                    p("Reconstruct a grinn network queried from the internal graph database, then compute and combine with a differential correlation network, see ",
                                                      a(href='http://kwanjeeraw.github.io/grinn/fetchgrinndiffcorr.html',target='_blank','here'),' for argument details.'
                                                    )
                                          )#end mainPanel
                          )),
                          wellPanel(
                            h4("Input arguments:"), helpText("* required field"),
                            h5("Network query:"),
                            fluidRow(
                              column(4, fileInput(inputId='txtInputDiff', label='txtInput *', accept=c('text/csv','text/comma-separated-values,text/plain','.csv'))),
                              column(8, mainPanel(tableOutput('txtExTableDiff')))
                            ),
                            fluidRow(
                              column(6,radioButtons("fromDiff", label = "from",inline = TRUE, 
                                                    choices = list("metabolite"="metabolite","protein"="protein","gene"="gene","pathway"="pathway"),selected = "metabolite")),
                              column(6,radioButtons("toDiff", label = "to",inline = TRUE, 
                                                    choices = list("metabolite"="metabolite","protein"="protein","gene"="gene","pathway"="pathway"),selected = "pathway"))
                            ),
                            fluidRow(
                              column(6,checkboxGroupInput("filterSourceDiff2", label ="filterSource",inline = TRUE,
                                                          choices = list("SMPDB"="SMPDB","KEGG"="KEGG","REACTOME"="REACTOME"),selected = c("SMPDB","KEGG","REACTOME")))
                            ),
                            fluidRow(
                              column(12,radioButtons("dbXrefDiff", label ="dbXref",inline = TRUE,
                                                     choices = list("grinn"="grinn","chebi"="chebi","kegg"="kegg","pubchem"="pubchem",
                                                                    "inchi"="inchi","hmdb"="hmdb","smpdb"="smpdb","reactome"="reactome",
                                                                    "uniprot"="uniprot","ensembl"="ensembl","entrezgene"="entrezgene"),selected = "grinn"))
                            ),
                            hr(),
                            h5("Correlation analysis:"),
                            fluidRow(
                              column(8, radioButtons('sepDiff2', 'Delimiter',c(Comma=',',Tab='\t',Semicolon=';'),',',inline=TRUE))
                            ),
                            fluidRow(
                              column(4, fileInput(inputId='datXInputDiff2', label='datX1 *', accept=c('text/csv','text/comma-separated-values,text/plain','.csv'))),
                              column(8, mainPanel(tableOutput('datXExTableDiff2')))
                            ),
                            fluidRow(
                              column(4, fileInput(inputId='datX2InputDiff2', label='datX2 *', accept=c('text/csv','text/comma-separated-values,text/plain','.csv'))),
                              column(8, mainPanel(tableOutput('datX2ExTableDiff2')))
                            ),
                            hr(),
                            fluidRow(  
                              column(6,numericInput("pvalDiff2", label = "pDiff *",value = 0.05, min = 0, max = 0.05, step = 0.001))
                            ),
                            fluidRow(
                              column(6,radioButtons("methodDiff2", label ="method",inline = TRUE,
                                                    choices = list("spearman"="spearman","pearson"="pearson","kendall"="kendall"),selected = "spearman"))
                            ),
                            hr(),
                            actionButton("computeDiff2","Submit")
                          ),
                          fluidRow(column(12,
                                          h4("Output network"),
                                          sidebarLayout(
                                            sidebarPanel(width=3,
                                                         downloadButton('downloadEdgeDiff2', 'Download Edge'),br(),
                                                         downloadButton('downloadNodeDiff2', 'Download Node')
                                            ),
                                            mainPanel(width=9,    
                                                      tabsetPanel(type = "tabs", 
                                                                  tabPanel("Summary", verbatimTextOutput("consolemsgdiff2")),
                                                                  tabPanel("Nodes", dataTableOutput("nodeTableDiff2")),
                                                                  tabPanel("Edges", dataTableOutput("edgeTableDiff2"))
                                                      )
                                            )
                                          )#end sidebarLayout
                          ))        
                  ),
                  tabItem(tabName = "convertid",
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
                              column(4, fileInput(inputId='txtInputConv', label='txtInput *', accept=c('text/csv','text/comma-separated-values,text/plain','.csv'))),
                              column(8, mainPanel(tableOutput('txtExTableConv')))
                            ),
                            hr(),
                            fluidRow(
                              column(6,radioButtons("fromConv", label = "nodetype",inline = TRUE, 
                                                    choices = list("metabolite"="metabolite","protein"="protein","gene"="gene","pathway"="pathway"),selected = "metabolite"))
                            ),
                            fluidRow(
                              column(12,radioButtons("dbXrefConv", label ="dbXref",inline = TRUE,
                                                     choices = list("chebi"="chebi","kegg"="kegg","pubchem"="pubchem",
                                                                    "inchi"="inchi","hmdb"="hmdb","smpdb"="smpdb","reactome"="reactome",
                                                                    "uniprot"="uniprot","ensembl"="ensembl","entrezgene"="entrezgene"),selected = "kegg"))
                            ),
                            hr(),
                            actionButton("convert","Submit")
                          ),
                          fluidRow(column(12,
                                          h4("Output"),
                                          sidebarLayout(
                                            sidebarPanel(width=3,
                                                         downloadButton('downloadId', 'Download')
                                            ),
                                            mainPanel(width=9,    
                                                      tabsetPanel(type = "tabs", 
                                                                  tabPanel("Summary", verbatimTextOutput("consolemsgconv")),
                                                                  tabPanel("Table", dataTableOutput("idTable"))
                                                      )
                                            )
                                          )#end sidebarLayout
                          ))
                  ),
                  tabItem(tabName = "setdb",
                          fluidRow(column(12,
                                          mainPanel(width=12,
                                                    h3("setGrinnDb"),
                                                    p("Set the graph database location for the currently working environment, see ",
                                                      a(href='http://kwanjeeraw.github.io/grinn/setdb.html',target='_blank','here'),' for argument details.'
                                                    )
                                          )#end mainPanel
                          )),
                          wellPanel(
                            h4("Current database location:"),
                            fluidRow(
                              column(12, verbatimTextOutput("currentdb"))
                            ),
                            h4("Input arguments:"), helpText("* required field"),
                            fluidRow(
                              column(6, textInput(inputId='dburl', label='url *', value=""))
                            ),
                            hr(),
                            actionButton("submit","Submit"),  actionButton("resettingdb","Reset")
                          )
                  )
                )# end tabItems
              )# end dashboardBody
)# end dashboardPage