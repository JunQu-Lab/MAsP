#' MAsP Web Application
#'
#' This function initiates the R shiny app 'MAsP'.
#'
#'
#' @export

MAsPShiny <- function(){
  options(shiny.maxRequestSize=100000*1024^2)
  appDir <- system.file("shiny", "MAsP",package = "MAsP")
  if (appDir == "") {
    stop("Could not find UI directory. Try re-installing `MAsP`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}

