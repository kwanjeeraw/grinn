library(shiny)
options(shiny.maxRequestSize = 30*1024^2) #the upload file size limit is 30MB
shinyServer(function(input, output) {
  output$ui <- renderUI({
    if (is.null(input$fnCall))
      return()
    
    # Depending on input$input_type, we'll generate a different
    # UI component and send it to the client.
    switch(input$fnCall,
      "fetchGrinnNetwork" = source("layout/grinnnetwork.R", local=TRUE),
      "fetchCorrGrinnNetwork" = source("layout/corrgrinn.R", local=TRUE), 
      "fetchDiffCorrGrinnNetwork" = source("layout/diffgrinn.R", local=TRUE),
      "fetchGrinnCorrNetwork" = source("layout/grinncorr.R", local=TRUE), 
      "fetchGrinnDiffCorrNetwork" = source("layout/grinndiff.R", local=TRUE),
      "convertToGrinnID" = source("layout/convertid.R", local=TRUE),
      "setGrinnDb" = source("layout/setdb.R", local=TRUE)
    )
  })
  
  #functions
  idConverting <- function(datIn, fromXref){
    grinnID = convertToGrinnID(txtInput=colnames(datIn), nodetype=datIn[1,1], dbXref=fromXref) #call grinn function to convert ids
    grinnID = grinnID[!duplicated(grinnID[,1]),] #keep the first mapped id
    colnames(datIn) = lapply(colnames(datIn),function(x) ifelse(length(which(grinnID[,1] == x))>0,as.character(grinnID$GRINNID[which(grinnID[,1] == x)]),x))
    datIn
  }

  #inputs
  txtInput <- reactive({
    inFile <- input$txtInput
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
  datYInput <- reactive({
    datYFile <- input$datYInput
    if (is.null(datYFile))
      return(NULL)
    read.csv(datYFile$datapath, header=TRUE, row.names=1, stringsAsFactors=FALSE, sep = input$sep)
  })
  datX2Input <- reactive({
    datX2File <- input$datX2Input
    if (is.null(datX2File))
      return(NULL)
    read.csv(datX2File$datapath, header=TRUE, row.names=1, stringsAsFactors=FALSE, sep = input$sep)
  })
  datY2Input <- reactive({
    datY2File <- input$datY2Input
    if (is.null(datY2File))
      return(NULL)
    read.csv(datY2File$datapath, header=TRUE, row.names=1, stringsAsFactors=FALSE, sep = input$sep)
  })
  datPhenoInput <- reactive({
    datPhenoFile <- input$datPhenoInput
    if (is.null(datPhenoFile))
      return(NULL)
    read.csv(datPhenoFile$datapath, header=TRUE, row.names=1, stringsAsFactors=FALSE, sep = input$sep)
  })

  #execute functions
  values <- reactiveValues()
  exeGrinn <- function(){
    if(input$fnCall == "fetchGrinnNetwork"){
      values$myv<-NULL #reset
      values$myv = fetchGrinnNetwork(txtInput=unlist(txtInput()), from=input$from, to=input$to, filterSource=input$filterSource, dbXref=input$dbXref, returnAs="tab")
    }
    if(input$fnCall == "fetchCorrGrinnNetwork"){
      values$myv<-NULL #reset
      values$myv = fetchCorrGrinnNetwork(datX=datXInput(), datY=datYInput(), corrCoef=input$corrCoef, pval=input$pval, method=input$method, 
                                     xTo=xto(), yTo=yto(), filterSource=input$filterSource, returnAs="tab")
    }
    if(input$fnCall == "fetchDiffCorrGrinnNetwork"){
      values$myv<-NULL #reset
      values$myv = fetchDiffCorrGrinnNetwork(datX1=datXInput(), datX2=datX2Input(), datY1=datYInput(), datY2=datY2Input(), pDiff=input$pval, method=input$method, 
                                             xTo=xto(), yTo=yto(), filterSource=input$filterSource, returnAs="tab")
    }
    if(input$fnCall == "fetchGrinnCorrNetwork"){
      values$myv<-NULL #reset
      values$myv = fetchGrinnCorrNetwork(txtInput=unlist(txtInput()), from=input$from, to=input$to, filterSource=input$filterSource, dbXref=input$dbXref,
                                     datX=datXInput(), datY=datYInput(), corrCoef=input$corrCoef, pval=input$pval, method=input$method, returnAs="tab")
    }
    if(input$fnCall == "fetchGrinnDiffCorrNetwork"){
      values$myv<-NULL #reset
      values$myv = fetchGrinnDiffCorrNetwork(txtInput=unlist(txtInput()), from=input$from, to=input$to, filterSource=input$filterSource, dbXref=input$dbXref,
                                         datX1=datXInput(), datX2=datX2Input(), datY1=datYInput(), datY2=datY2Input(), pDiff=input$pval, method=input$method, returnAs="tab")
    }
    if(input$fnCall == "convertToGrinnID"){
      values$myv<-NULL #reset
      values$myv = convertToGrinnID(txtInput=unlist(txtInput()), nodetype=input$from, dbXref=input$dbXref)
    }
    if(input$fnCall == "setGrinnDb"){
      values$myv<-NULL #reset
      values$myv = setGrinnDb(url=input$dburl)
    }
  }
  returnResult <- reactive({
    exeGrinn()
    values$myv
  })
  #outputs
  output$txtExTable <- renderTable({
    head(txtInput(), n = 5)
  })
  output$datXExTable <- renderTable({
    head(datXInput()[,1:5], n = 5)
  })
  output$datYExTable <- renderTable({
    head(datYInput()[,1:5], n = 5)
  })
  output$datX2ExTable <- renderTable({
    head(datX2Input()[,1:5], n = 5)
  })
  output$datY2ExTable <- renderTable({
    head(datY2Input()[,1:5], n = 5)
  })
  output$datPhenoExTable <- renderTable({
    head(datPhenoInput(), n = 5)
  })
  output$summaryCode <- renderPrint({
    if (input$submit == 0) return()
    isolate({
      # Create a Progress object
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      progress$set(message = "Computing network", value = 0)
      for (i in 1:10) {
        # Increment the progress bar, and update the detail text.
        progress$inc(0.1, detail = "...")
        # Pause for 0.1 seconds to simulate a long computation.
        #Sys.sleep(0.1)
      }
      summary(returnResult())
    })
  })
  output$nodeTable <- renderTable({
    if (input$submit == 0) return()
    isolate({
      head(returnResult()$nodes[,c(1,2,4)], n=10)      
    })
  })
  output$edgeTable <- renderTable({
    if (input$submit == 0) return()
    isolate({
      head(returnResult()$edges, n=10)
    })
  })
  output$idTable <- renderTable({
    if (input$submit == 0) return()
    isolate({
      head(returnResult(), n=10)
    })
  })
  output$downloadEdge <- downloadHandler(
    filename = function() { 
      paste('grinnOutEdge.txt') 
    },
    content = function(file) {
      write.table(as.matrix(returnResult()$edges), file, sep='\t', row.names = F, quote = FALSE)
    }
  )
  output$downloadNode <- downloadHandler(
    filename = function() { 
      paste('grinnOutNode.txt') 
    },
    content = function(file) {
      write.table(as.matrix(returnResult()$nodes), file, sep='\t', row.names = F, quote = FALSE)
    }
  )
  output$downloadId <- downloadHandler(
    filename = function() { 
      paste('grinnOutId.txt') 
    },
    content = function(file) {
      write.table(as.matrix(returnResult()), file, sep='\t', row.names = F, quote = FALSE)
    }
  )
  output$currentdb <- renderPrint({
    if (input$submit == 0) return(nld)
    isolate({
      # Create a Progress object
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      progress$set(message = "Computing network", value = 0)
      for (i in 1:10) {
        # Increment the progress bar, and update the detail text.
        progress$inc(0.1, detail = "...")
        # Pause for 0.1 seconds to simulate a long computation.
        #Sys.sleep(0.1)
      }
      returnResult()
    })
  })
})#end shinyServer