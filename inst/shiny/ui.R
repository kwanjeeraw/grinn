library(shiny)
library(DT)
library(grinn)

shinyUI(fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
    tags$link(rel = "shortcut icon", type="image/x-icon", href="favicon.ico")
  ),
  titlePanel(title=img(src="logo.png", width = 250), windowTitle="Grinn"),
  fluidRow(
    column(3, wellPanel(
      selectInput("fnCall", "Select function",
        c("fetchGrinnNetwork", "fetchCorrGrinnNetwork", "fetchDiffCorrGrinnNetwork", 
          "fetchGrinnCorrNetwork", "fetchGrinnDiffCorrNetwork",
          "convertToGrinnID","setGrinnDb")
      ),
      p('Click to visit',a(href='http://kwanjeeraw.github.io/grinn/',target='_blank','homepage')),
      br()
    )),
    column(9,
      # This outputs the dynamic UI component
      uiOutput("ui")
    )
  )#end fluidRow
))#end shinyUI, fluidPage