library(shiny)
options(shiny.maxRequestSize = 30*1024^2) #the upload file size limit is 30MB
shinyServer(function(input, output) {
  #functions
  idConverting <- function(datIn, fromXref){
    grinnID = convertToGrinnID(txtInput=colnames(datIn), nodetype=datIn[1,1], dbXref=fromXref) #call grinn function to convert ids
    grinnID = grinnID[!duplicated(grinnID[,1]),] #keep the first mapped id
    colnames(datIn) = lapply(colnames(datIn),function(x) ifelse(length(which(grinnID[,1] == x))>0,as.character(grinnID$GRINNID[which(grinnID[,1] == x)]),x))
    datIn
  }
  
  #reactiveValues
  values <- reactiveValues()
  
  #inputs
  txtInput <- reactive({
    inFile <- input$txtInput
    if (is.null(inFile))
      return(NULL)
    read.csv(inFile$datapath, header=FALSE, stringsAsFactors=FALSE)
  })
  txtInputCorr <- reactive({
    inFile <- input$txtInputCorr
    if (is.null(inFile))
      return(NULL)
    read.csv(inFile$datapath, header=FALSE, stringsAsFactors=FALSE)
  })
  txtInputDiff <- reactive({
    inFile <- input$txtInputDiff
    if (is.null(inFile))
      return(NULL)
    read.csv(inFile$datapath, header=FALSE, stringsAsFactors=FALSE)
  })
  txtInputConv <- reactive({
    inFile <- input$txtInputConv
    if (is.null(inFile))
      return(NULL)
    read.csv(inFile$datapath, header=FALSE, stringsAsFactors=FALSE)
  })
  xto <- reactive({
   if(input$xTo == "non"){
     return(NULL)
   }else{
     return(input$xTo)
   }
  })
  xtoDiff <- reactive({
    if(input$xToDiff == "non"){
      return(NULL)
    }else{
      return(input$xToDiff)
    }
  })
  yto <- reactive({
    if(input$yTo == "non"){
      return(NULL)
    }else{
      return(input$yTo)
    }
  })
  datXInput <- reactive({
    datXFile <- input$datXInput
    if (is.null(datXFile))
      return(NULL)
    read.csv(datXFile$datapath, header=TRUE, row.names=1, stringsAsFactors=FALSE, sep = input$sep)
  })
  datXInputDiff <- reactive({
    datXFile <- input$datXInputDiff
    if (is.null(datXFile))
      return(NULL)
    read.csv(datXFile$datapath, header=TRUE, row.names=1, stringsAsFactors=FALSE, sep = input$sepDiff)
  })
  datXInputCorr <- reactive({
    datXFile <- input$datXInputCorr
    if (is.null(datXFile))
      return(NULL)
    read.csv(datXFile$datapath, header=TRUE, row.names=1, stringsAsFactors=FALSE, sep = input$sepCorr)
  })
  datXInputDiff2 <- reactive({
    datXFile <- input$datXInputDiff2
    if (is.null(datXFile))
      return(NULL)
    read.csv(datXFile$datapath, header=TRUE, row.names=1, stringsAsFactors=FALSE, sep = input$sepDiff2)
  })
  datYInput <- reactive({
    datYFile <- input$datYInput
    if (is.null(datYFile))
      return(NULL)
    read.csv(datYFile$datapath, header=TRUE, row.names=1, stringsAsFactors=FALSE, sep = input$sep)
  })
  datYInputCorr <- reactive({
    datYFile <- input$datYInputCorr
    if (is.null(datYFile))
      return(NULL)
    read.csv(datYFile$datapath, header=TRUE, row.names=1, stringsAsFactors=FALSE, sep = input$sepCorr)
  })
  datX2Input <- reactive({
    datX2File <- input$datX2Input
    if (is.null(datX2File))
      return(NULL)
    read.csv(datX2File$datapath, header=TRUE, row.names=1, stringsAsFactors=FALSE, sep = input$sepDiff)
  })
  datX2InputDiff2 <- reactive({
    datX2File <- input$datX2InputDiff2
    if (is.null(datX2File))
      return(NULL)
    read.csv(datX2File$datapath, header=TRUE, row.names=1, stringsAsFactors=FALSE, sep = input$sepDiff2)
  })
  datPhenoInput <- reactive({
    datPhenoFile <- input$datPhenoInput
    if (is.null(datPhenoFile))
      return(NULL)
    read.csv(datPhenoFile$datapath, header=TRUE, row.names=1, stringsAsFactors=FALSE, sep = input$sep)
  })

  #execute functions
  observeEvent(input$query, {
    isolate({
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Querying network", value = 0)
      for (i in 1:10) {
        progress$inc(1/10, detail = "...")
        Sys.sleep(0.2)
      }
    })
    values$bioNetwork = fetchGrinnNetwork(txtInput=unlist(txtInput()), from=input$from, to=input$to, filterSource=input$filterSource, dbXref=input$dbXref, returnAs="tab")
  })
  observeEvent(input$computeCorr, {
    isolate({
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Computing network", value = 0)
      for (i in 1:10) {
        progress$inc(1/10, detail = "...")
        Sys.sleep(0.2)
      }
    })
    values$corrNetwork = fetchCorrGrinnNetwork(datX=datXInput(), datY=datYInput(), corrCoef=input$corrCoef, pval=input$pval, method=input$method, 
                                              xTo=xto(), yTo=yto(), filterSource=input$filterSourceCorr, returnAs="tab")
  })
  observeEvent(input$computeDiff, {
    isolate({
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Computing network", value = 0)
      for (i in 1:10) {
        progress$inc(1/10, detail = "...")
        Sys.sleep(0.2)
      }
    })
    values$diffNetwork = fetchDiffCorrGrinnNetwork(datX1=datXInputDiff(), datX2=datX2Input(), pDiff=input$pvalDiff, method=input$methodDiff, 
                                                   xTo=xtoDiff(), filterSource=input$filterSourceDiff, returnAs="tab")
  })
  observeEvent(input$computeCorr2, {
    isolate({
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Computing network", value = 0)
      for (i in 1:10) {
        progress$inc(1/10, detail = "...")
        Sys.sleep(0.2)
      }
    })
    values$corr2Network = fetchGrinnCorrNetwork(txtInput=unlist(txtInputCorr()), from=input$fromCorr, to=input$toCorr, filterSource=input$filterSourceCorr2, dbXref=input$dbXrefCorr,
                                               datX=datXInputCorr(), datY=datYInputCorr(), corrCoef=input$corrCoefCorr, pval=input$pvalCorr, method=input$methodCorr, returnAs="tab")
  })
  observeEvent(input$computeDiff2, {
    isolate({
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Computing network", value = 0)
      for (i in 1:10) {
        progress$inc(1/10, detail = "...")
        Sys.sleep(0.2)
      }
    })
    values$diff2Network = fetchGrinnDiffCorrNetwork(txtInput=unlist(txtInputDiff()), from=input$fromDiff, to=input$toDiff, filterSource=input$filterSourceDiff2, dbXref=input$dbXrefDiff,
                                                    datX1=datXInputDiff2(), datX2=datX2InputDiff2(), pDiff=input$pvalDiff2, method=input$methodDiff2, returnAs="tab")
  })
  observeEvent(input$convert, {
    isolate({
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Converting ids", value = 0)
      for (i in 1:10) {
        progress$inc(1/10, detail = "...")
        Sys.sleep(0.2)
      }
    })
    values$conv = convertToGrinnID(txtInput=unlist(txtInputConv()), nodetype=input$fromConv, dbXref=input$dbXrefConv)
  })
  observeEvent(input$submit, {
    isolate({
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Setting database", value = 0)
      for (i in 1:10) {
        progress$inc(1/10, detail = "...")
        Sys.sleep(0.1)
      }
    })
    values$currdb = setGrinnDb(url=input$dburl)
  })

  #outputs
  output$txtExTable <- renderTable({
    head(txtInput(), n = 5)
  })
  output$txtExTableCorr <- renderTable({
    head(txtInputCorr(), n = 5)
  })
  output$txtExTableDiff <- renderTable({
    head(txtInputDiff(), n = 5)
  })
  output$txtExTableConv <- renderTable({
    head(txtInputConv(), n = 5)
  })
  output$datXExTable <- renderTable({
    head(datXInput()[,1:5], n = 5)
  })
  output$datXExTableDiff <- renderTable({
    head(datXInputDiff()[,1:5], n = 5)
  })
  output$datXExTableCorr <- renderTable({
    head(datXInputCorr()[,1:5], n = 5)
  })
  output$datXExTableDiff2 <- renderTable({
    head(datXInputDiff2()[,1:5], n = 5)
  })
  output$datYExTable <- renderTable({
    head(datYInput()[,1:5], n = 5)
  })
  output$datYExTableCorr <- renderTable({
    head(datYInputCorr()[,1:5], n = 5)
  })
  output$datX2ExTable <- renderTable({
    head(datX2Input()[,1:5], n = 5)
  })
  output$datX2ExTableDiff2 <- renderTable({
    head(datX2InputDiff2()[,1:5], n = 5)
  })
  output$datPhenoExTable <- renderTable({
    head(datPhenoInput(), n = 5)
  })
  output$consolemsgbio <- renderPrint({
    if (input$query == 0) return()
    summary(values$bioNetwork)
  })
  output$nodeTableBio <- renderDataTable(datatable(
    values$bioNetwork$nodes[,-4], rownames = FALSE, extensions = 'Scroller', options = list(pageLength = 25, scrollY = 600, scrollX = TRUE, scrollCollapse = TRUE)
  ))
  output$edgeTableBio <- renderDataTable(datatable(
    values$bioNetwork$edges, rownames = FALSE, extensions = 'Scroller', options = list(pageLength = 25, scrollY = 600, scrollX = TRUE, scrollCollapse = TRUE)
  ))
  output$downloadEdgeBio <- downloadHandler(
    filename = function() { 
      paste('grinnOutEdge.txt') 
    },
    content = function(file) {
      write.table(as.matrix(values$bioNetwork$edges), file, sep='\t', row.names = F, quote = FALSE)
    }
  )
  output$downloadNodeBio <- downloadHandler(
    filename = function() { 
      paste('grinnOutNode.txt') 
    },
    content = function(file) {
      write.table(as.matrix(values$bioNetwork$nodes), file, sep='\t', row.names = F, quote = FALSE)
    }
  )
  output$consolemsgcorr <- renderPrint({
    if (input$computeCorr == 0) return()
    summary(values$corrNetwork)
  })
  output$nodeTableCorr <- renderDataTable(datatable(
    values$corrNetwork$nodes[,-4], rownames = FALSE, extensions = 'Scroller', options = list(pageLength = 25, scrollY = 600, scrollX = TRUE, scrollCollapse = TRUE)
  ))
  output$edgeTableCorr <- renderDataTable(datatable(
    values$corrNetwork$edges, rownames = FALSE, extensions = 'Scroller', options = list(pageLength = 25, scrollY = 600, scrollX = TRUE, scrollCollapse = TRUE)
  ))
  output$downloadEdgeCorr <- downloadHandler(
    filename = function() { 
      paste('corrOutEdge.txt') 
    },
    content = function(file) {
      write.table(as.matrix(values$corrNetwork$edges), file, sep='\t', row.names = F, quote = FALSE)
    }
  )
  output$downloadNodeCorr <- downloadHandler(
    filename = function() { 
      paste('corrOutNode.txt') 
    },
    content = function(file) {
      write.table(as.matrix(values$corrNetwork$nodes), file, sep='\t', row.names = F, quote = FALSE)
    }
  )
  output$consolemsgdiff <- renderPrint({
    if (input$computeDiff == 0) return()
    summary(values$diffNetwork)
  })
  output$nodeTableDiff <- renderDataTable(datatable(
    values$diffNetwork$nodes[,-4], rownames = FALSE, extensions = 'Scroller', options = list(pageLength = 25, scrollY = 600, scrollX = TRUE, scrollCollapse = TRUE)
  ))
  output$edgeTableDiff <- renderDataTable(datatable(
    values$diffNetwork$edges, rownames = FALSE, extensions = 'Scroller', options = list(pageLength = 25, scrollY = 600, scrollX = TRUE, scrollCollapse = TRUE)
  ))
  output$downloadEdgeDiff <- downloadHandler(
    filename = function() { 
      paste('diffOutEdge.txt') 
    },
    content = function(file) {
      write.table(as.matrix(values$diffNetwork$edges), file, sep='\t', row.names = F, quote = FALSE)
    }
  )
  output$downloadNodeDiff <- downloadHandler(
    filename = function() { 
      paste('diffOutNode.txt') 
    },
    content = function(file) {
      write.table(as.matrix(values$diffNetwork$nodes), file, sep='\t', row.names = F, quote = FALSE)
    }
  )
  output$consolemsgcorr2 <- renderPrint({
    if (input$computeCorr2 == 0) return()
    summary(values$corr2Network)
  })
  output$nodeTableCorr2 <- renderDataTable(datatable(
    values$corr2Network$nodes[,-4], rownames = FALSE, extensions = 'Scroller', options = list(pageLength = 25, scrollY = 600, scrollX = TRUE, scrollCollapse = TRUE)
  ))
  output$edgeTableCorr2 <- renderDataTable(datatable(
    values$corr2Network$edges, rownames = FALSE, extensions = 'Scroller', options = list(pageLength = 25, scrollY = 600, scrollX = TRUE, scrollCollapse = TRUE)
  ))
  output$downloadEdgeCorr2 <- downloadHandler(
    filename = function() { 
      paste('corrOutEdge.txt') 
    },
    content = function(file) {
      write.table(as.matrix(values$corr2Network$edges), file, sep='\t', row.names = F, quote = FALSE)
    }
  )
  output$downloadNodeCorr2 <- downloadHandler(
    filename = function() { 
      paste('corrOutNode.txt') 
    },
    content = function(file) {
      write.table(as.matrix(values$corr2Network$nodes), file, sep='\t', row.names = F, quote = FALSE)
    }
  )
  output$consolemsgdiff2 <- renderPrint({
    if (input$computeDiff2 == 0) return()
    summary(values$diff2Network)
  })
  output$nodeTableDiff2 <- renderDataTable(datatable(
    values$diff2Network$nodes[,-4], rownames = FALSE, extensions = 'Scroller', options = list(pageLength = 25, scrollY = 600, scrollX = TRUE, scrollCollapse = TRUE)
  ))
  output$edgeTableDiff2 <- renderDataTable(datatable(
    values$diff2Network$edges, rownames = FALSE, extensions = 'Scroller', options = list(pageLength = 25, scrollY = 600, scrollX = TRUE, scrollCollapse = TRUE)
  ))
  output$downloadEdgeDiff2 <- downloadHandler(
    filename = function() { 
      paste('diffOutEdge.txt') 
    },
    content = function(file) {
      write.table(as.matrix(values$diff2Network$edges), file, sep='\t', row.names = F, quote = FALSE)
    }
  )
  output$downloadNodeDiff2 <- downloadHandler(
    filename = function() { 
      paste('diffOutNode.txt') 
    },
    content = function(file) {
      write.table(as.matrix(values$diff2Network$nodes), file, sep='\t', row.names = F, quote = FALSE)
    }
  )
  output$consolemsgconv <- renderPrint({
    if (input$convert == 0) return()
    summary(values$conv)
  })
  output$idTable <- renderDataTable(datatable(
    values$conv, rownames = FALSE, extensions = 'Scroller', options = list(pageLength = 25, scrollY = 600, scrollX = TRUE, scrollCollapse = TRUE)
  ))
  output$downloadId <- downloadHandler(
    filename = function() { 
      paste('grinnOutId.txt') 
    },
    content = function(file) {
      write.table(as.matrix(values$conv), file, sep='\t', row.names = F, quote = FALSE)
    }
  )
  output$currentdb <- renderPrint({
    if (input$submit == 0 || input$resettingdb == 0) return(nld)
    values$currdb
  })
  #reset
  observeEvent(input$resettingdb, {
    #reset db url to default
    isolate({
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Resetting to default database", value = 0)
      for (i in 1:10) {
        progress$inc(1/10, detail = "...")
        Sys.sleep(0.1)
      }
    })
    values$currdb = setGraphDb(url="http://grinn.genomecenter.ucdavis.edu:7474/db/data/")
  })
})#end shinyServer