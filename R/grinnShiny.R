#' @title grinnShiny
#' @return start shiny app
#' @import shiny DT
#' @export
grinnShiny<-function(){
  shiny::runApp(system.file("shiny", package = "grinn"),
                display.mode = "normal",
                launch.browser = TRUE)
}
