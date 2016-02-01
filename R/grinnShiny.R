#'Start GUI
#'@description start GUI as shiny app.
#'@usage grinnShiny()
#'@return shiny app
#'@import shiny DT
#'@export
grinnShiny<-function(){
  shiny::runApp(system.file("shiny", package = "grinn"),
                display.mode = "normal",
                launch.browser = TRUE)
}
