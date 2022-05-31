#' Location Preview Plot
#'
#' Make a location preview plot for user to directly see the map of sample slots
#'
#' @param location the formatted location information
#' @return a location preview plot
#' @export


MyPreview <- function(location){
  new_loc <- loc_convert(location)
  new_loc$row <- factor(new_loc$row, levels = base::sort(unique(new_loc$row), decreasing = TRUE))
  new_loc$column <- factor(new_loc$column, levels = base::sort(unique(new_loc$column), decreasing = FALSE))
  output <- ggplot(new_loc, aes(column, row)) +
    geom_text(aes(label = Name), size = 6) +
    scale_x_discrete(breaks = 1:(ncol(location)+1)) +
    scale_y_discrete(breaks = 1:(nrow(location)+1)) +
    coord_fixed() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="right",panel.grid.minor=element_blank(),
          panel.background=element_rect(fill=NA, colour=NA),
          panel.border=element_rect(size=1,fill=NA,colour="black")) +
    theme(axis.text=element_text(size=14, family="Helvetica", colour="black"),
          axis.title=element_text(size=14, family="Helvetica",colour="black"),
          legend.text=element_text(size=14, family="Helvetica",colour="black"),
          legend.title=element_text(size=14, family="Helvetica", colour="black"))
  return(output)
}

