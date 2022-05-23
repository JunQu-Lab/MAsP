#' Distribution plot
#'
#' Distribution plot and its batch mode
#'
#' @param df the quantitative protein result
#' @param location the original location information
#' @param target_name the protein in the plot
#' @param limit a logical value, if setting the value limits of protein abundance/ratio in the distribution plot
#' @param scale a logical value, if converting data to z-score
#' @param upperL the upperbound of the protein abundance/ratio the plot can show
#' @param lowerL the lowerbound of the protein abundance/ratio the plot can show
#' @param color1 the lower color, color name and HTML code are both acceptable
#' @param color2 the middle color, color name and HTML code are both acceptable
#' @param color3 the upper color, color name and HTML code are both acceptable
#' @param transparent a logical value, if the plot is with a transparent background while downloading
#' @param mask a logical value, if adding a brain-shape cover on the plot
#' @param Legend a logical value, if showing the legend
#' @param g the brain-shape cover
#' @return a protein distribution plot
#' @export
#'
#'
#'



MyDistriPlot <- function(df, location, target_name, limit, scale = TRUE, upperL = 2, lowerL = -2,
                         color1, color2, color3, transparent = FALSE, mask = FALSE, Legend = TRUE, g = g){
  #outputs <- list()
  new_loc <- loc_convert(location)
  rownames(df) <- toupper(rownames(df))
  target_name <- toupper(target_name)
  target <- df[grep(target_name, rownames(df)),]
  target <- target[1,]
  legend_label <- 'Protein Relative Abundance'
  if(scale){
    tmp_value <- as.numeric(target)
    tmp_value <- tmp_value-mean(tmp_value, na.rm = TRUE)/sd(tmp_value, na.rm = TRUE)
    target[1,] <- tmp_value
    legend_label = "Protein Z-score"
  }
  new_loc$value <- NA
  for (i in 1:nrow(new_loc)) {
    tmp <- which(new_loc$Name[i] == names(target))
    if(length(tmp) != 0)
      new_loc$value[i] <- target[tmp]
  }
  new_loc$value <- as.numeric(new_loc$value)
  if(limit == 'none'){
    range_tmp <- range(as.numeric(new_loc$value), na.rm = TRUE)
    UpperLimit <- range_tmp[2]
    LowerLimit <- range_tmp[1]
  }
  if(limit == 'outlierRM'){
    target <- target[!is.na(target)]
    outliers <- boxplot.stats(as.numeric(target))$out
    dataRange <- target[!as.numeric(target) %in% outliers]
    UpperLimit <- range(dataRange)[2]
    LowerLimit <- range(dataRange)[1]
  }
  if(limit == 'customized'){
    UpperLimit <- as.numeric(upperL)
    LowerLimit <- as.numeric(lowerL)
  }

  for (i in 1:nrow(new_loc)) {
    tmp <-new_loc$value[i]
    if(!is.na(tmp)){
      if(tmp > UpperLimit){
        new_loc$value[i] <- UpperLimit
      }
      if(tmp < LowerLimit){
        new_loc$value[i] <- LowerLimit
      }
    }
  }

  if(transparent) transparent2 <- "transparent" else transparent2 <- NA
  if(Legend) legend2 <- 'right' else legend2 <- 'none'


  # read png mask

  PlotTitle <- paste('Heatmap for Protein', rownames(target))
  new_loc$row <- factor(new_loc$row, levels = base::sort(unique(new_loc$row), decreasing = TRUE))
  new_loc$column <- factor(new_loc$column, levels = base::sort(unique(new_loc$column), decreasing = FALSE))

  if(mask){
    output <- ggplot(new_loc, aes(column, row, fill= value)) +
      geom_tile() + xlab('Columns') + ylab('Rows') + labs(fill = legend_label) +
      scale_fill_gradientn(colors = c(color1,color2,color3),na.value = '#A9A9A9') +
      scale_x_discrete(breaks = 1:(ncol(location)+1)) +
      scale_y_discrete(breaks = 1:(nrow(location)+1)) +
      annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      coord_fixed() +
      ggtitle(PlotTitle) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(legend.position=legend2,panel.grid.minor=element_blank(),
            panel.border=element_rect(size=1,fill=NA,colour="black"),
            panel.background = element_rect(fill = "#A9A9A9",colour = NA),
            plot.background = element_rect(fill = transparent2,colour = NA),
            panel.grid.major = element_blank()) +
      theme(axis.text=element_text(size=14, family="Helvetica", colour="black"),
            axis.title=element_text(size=14, family="Helvetica",colour="black"),
            legend.text=element_text(size=14, family="Helvetica",colour="black"),
            legend.title=element_text(size=14, family="Helvetica", colour="black"))
  }else{
    output <- ggplot(new_loc, aes(column, row, fill= value)) +
      geom_tile() + xlab('Columns') + ylab('Rows') + labs(fill = legend_label) +
      scale_fill_gradientn(colors = c(color1,color2,color3),na.value = '#A9A9A9') +
      scale_x_discrete(breaks = 1:(ncol(location)+1)) +
      scale_y_discrete(breaks = 1:(nrow(location)+1)) +
      coord_fixed() +
      ggtitle(PlotTitle) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(legend.position=legend2,panel.grid.minor=element_blank(),
            panel.border=element_rect(size=1,fill=NA,colour="black"),
            panel.background = element_rect(fill = "#A9A9A9",colour = NA),
            plot.background = element_rect(fill = transparent2,colour = NA),
            panel.grid.major = element_blank()) +
      theme(axis.text=element_text(size=14, family="Helvetica", colour="black"),
            axis.title=element_text(size=14, family="Helvetica",colour="black"),
            legend.text=element_text(size=14, family="Helvetica",colour="black"),
            legend.title=element_text(size=14, family="Helvetica", colour="black"))
  }

  return(output)
}


MyDistriPlot_batch <- function(df, location, proteins_csv, directory, limit, scale = TRUE, upperL = 2, lowerL = -2, mask = FALSE,
                               color1, color2, color3, transparent = FALSE, Legend = TRUE, width, height, g = g){
  withProgress(message = 'Batch Plots', value = 0, {
    Proteins <- unname(unlist(proteins_csv))
    current_dir <- getwd()
    setwd(directory)
    for (i in 1:length(Proteins)) {
      incProgress(i/length(Proteins), detail = 'generating...')
      output <- MyDistriPlot(df = df,location = location, target_name = Proteins[i], limit = limit, scale = scale, upperL = upperL, lowerL = lowerL,
                             mask = mask, color1 = color1, color2 = color2, color3 = color3,
                             transparent = transparent, Legend = Legend, g = g)
      ggsave(filename = paste('Distribution_of_', Proteins[i],'.png',sep = ''), plot = output, device='png', width = width, height = height, dpi = 300)
    }
  })
}
