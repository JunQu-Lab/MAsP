#' Location Formatting Function
#'
#' This function converts the location information to an coding-tractable format
#'
#' @param location the location information of each micro-slot
#' @return a new format that the codes are easy to deal with
#' @export
#'
#'
#'
#'



loc_convert <- function(location){
  loc_vector <- base::unname(unlist(location))
  loc_vector <- loc_vector[!(loc_vector == '')]
  output <- data.frame(Name = loc_vector, row = NA, column = NA)
  for (i in 1:length(loc_vector)) {
    nrow_location <- nrow(location)
    tmp <- which(location == loc_vector[i])
    tmp_row <- tmp %% nrow_location
    tmp_col <- tmp %/% nrow_location + 1
    if(tmp_row == 0){
      tmp_row = nrow_location
      tmp_col = tmp_col-1
    }
    output$row[i] <- tmp_row
    output$column[i] <- tmp_col
  }
  return(output)
}
