#' Correlation plots
#'
#'
#'
#' @param df_filtered the filtered protein dataset
#' @return a PCA plot colored by label from spectral clustering, and a dataset containing cluster label.
#' @export
#'
#'
#'
#'
#'

correlation_plot <- function(df, target_name, type = 'pearson', location, limit, upperL = 2, lowerL = -2,
                             color1, color2, color3, transparent = FALSE, mask = FALSE, g = g){
  withProgress(message = 'Correlation ', value = 0, {
    incProgress(0.5, detail = 'calculating correlations...')
    rownames(df) <- toupper(rownames(df))
    target_name <- toupper(target_name)
    output <- data.frame(protein_name = rownames(df), value = NA)
    target <- grep(target_name, output$protein_name)
    target <- target[1]
    for (i in 1:nrow(df)) {
      tmp <- rbind(df[target,], df[i,])
      tmp <- tmp[,colSums(is.na(tmp)) == 0] # remove NA values
      if(type == 'pearson') output$value[i] <- cor(as.numeric(tmp[1,]),as.numeric(tmp[2,]))
      if(type == 'cosine') output$value[i] <- lsa::cosine(as.numeric(tmp[1,]),as.numeric(tmp[2,]))
    }

    output <- output[-target,]
    output$abs_value <- abs(output$value)
    output <- output[order(output$abs_value, decreasing = TRUE),]
    output$Index <- 1:nrow(output)

    if(type == 'pearson') tabname <- 'Absolute Value of Pearson Correlation'
    if(type == 'cosine') tabname <- 'Absolute Value of Cosine Similarity'

    corr_plot <- ggplot(output, aes(x = Index, y = abs_value)) +
      geom_point() + xlab('Protein Rank') + ylab(tabname) +
      theme(legend.position="right",panel.grid.minor=element_blank(),
            panel.background=element_rect(fill=NA, colour=NA),
            panel.border=element_rect(size=1,fill=NA,colour="black")) +
      theme(axis.text=element_text(size=12, family="Helvetica", colour="black"),
            axis.title=element_text(size=12, family="Helvetica",colour="black"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.text=element_text(size=12, family="Helvetica",colour="black"),
            legend.title=element_text(size=12, family="Helvetica", colour="black"))

    TOP3 <- output$protein_name[1:3]
    proteins <- c(target_name, TOP3)
    Legend = FALSE
    plots <- list()
    for (i in 1:4) {
      plots[[i]] <- MyDistriPlot(df = df,location = location, target_name = proteins[i], limit = limit, upperL = upperL, lowerL = lowerL,
                                 mask = mask, color1 = color1, color2 = color2, color3 = color3,
                                 transparent = transparent, Legend = Legend, g = g)
    }

    outputs <- list()
    outputs$corr_plot <- corr_plot
    outputs$corr_data <- output
    outputs$total_plots <- plots
    return(outputs)
  })
}
