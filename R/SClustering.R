#' Spectral Clustering on filtered data
#'
#'
#'
#' @param df_filtered the filtered protein dataset
#' @return a PCA plot colored by label from spectral clustering, and a dataset containing cluster label.
#' @export



SClustering <- function(df_filtered){
  withProgress(message = 'Spectral Clustering', value = 0, {
    incProgress(0.5, detail = 'Perform Spectral Clustering...')
    outputs <- list()
    df_filtered <- na.omit(df_filtered)
    test1 <- Spectrum(t(df_filtered))
    SClabel <- test1$assignments
    PCA_result <- prcomp(df_filtered,center=TRUE,scale.=TRUE)
    PC1_VE <- format(round(PCA_result$sdev[1]/sum(PCA_result$sdev) * 100, 2), nsmall = 2)
    PC1_VE <- paste('(',PC1_VE,'%',')',sep = '')
    PC2_VE <- format(round(PCA_result$sdev[2]/sum(PCA_result$sdev) * 100, 2), nsmall = 2)
    PC2_VE <- paste('(',PC2_VE,'%',')',sep = '')

    plotdata <- as.data.frame(cbind(PCA_result$x[,1:2],SClabel))
    plotdata$SClabel <- as.factor(plotdata$SClabel)
    PCA_plot <- ggplot(plotdata, aes(x=PC1, y=PC2, color = SClabel)) +
      geom_point() +
      ggtitle('PCA plot for Spectrum Clustering result') +
      xlab(paste('PC1',PC1_VE, sep = ' ')) +
      ylab(paste('PC2',PC2_VE, sep = ' ')) +
      theme(panel.grid.minor=element_blank(),
            #panel.grid.major=element_line(colour="gray80",linetype="dashed"),
            panel.background=element_rect(fill=NA, colour=NA),
            panel.border=element_rect(size=1,fill=NA,colour="black")) +
      theme(axis.text=element_text(size=14, family="Helvetica", colour="black"),
            axis.title=element_text(size=14, family="Helvetica",colour="black"),
            legend.text=element_text(size=14, family="Helvetica",colour="black"),
            legend.title=element_text(size=14, family="Helvetica", colour="black"))

    output <- cbind(df_filtered, SClabel)
    outputs$df_label <- output
    outputs$PCA_plot <- PCA_plot
    return(outputs)
  })
}
