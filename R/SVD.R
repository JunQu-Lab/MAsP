#' Pattern selection by SVD
#'
#' Make a location preview plot for user to directly see the distribution of sample slots
#'
#' @param location the formatted location information
#' @return a location preview plot
#' @export
imputation <- function(target, location_corr){
  target <- data.frame(Name = names(target), ratio = as.numeric(target))
  location_corr$value <- NA
  for (i in 1:nrow(location_corr)) {
    tmp <- which(target$Name == location_corr$Name[i])
    if(length(tmp) != 0) location_corr$value[i] <- target$ratio[tmp]
  }
  imputation_obs <- which(is.na(location_corr$value))
  for (i in imputation_obs) {
    Row <- location_corr$row[i]
    Column <- location_corr$column[i]
    tmp <- location_corr[location_corr$row %in% (Row-1):(Row+1) & location_corr$column %in% (Column-1):(Column+1),]
    location_corr$value[i] <- mean(as.numeric(tmp$value), na.rm = TRUE)
  }
  target_impu <- as.vector(as.numeric(location_corr$value))
  return(target_impu)
}

#' @export
imputation_df <- function(df, location_corr){
  output <- data.frame(matrix(NA, ncol = nrow(location_corr), nrow = nrow(df)))
  names(output) <- location_corr$Name
  for (i in 1:nrow(df)) {
    output[i,] <- imputation(df[i,], location = location_corr)
  }
  return(output)
}

#' @export
SVD_plot <- function(df, location){
  outputs <- list()
  withProgress(message = 'VE plot', value = 0, {

    df <- as.matrix(df)
    for (i in 1:nrow(df)) {
      df[i,] <- (df[i,]-mean(df[i,],na.rm=TRUE))/sd(df[i,],na.rm = TRUE)
    }
    df <- as.data.frame(df)
    incProgress(0.2, detail = 'Missing value imputation...')
    location_corr <- loc_convert(location)
    df_impu <- imputation_df(df, location_corr = location_corr)
    rownames(df_impu) <- rownames(df)
    outputs$df_impu <- df_impu

    incProgress(0.5, detail = 'Perform SVD...')
    df_VE <- data.frame(protein_name = rownames(df_impu), VE = NA)
    location_tmp <- as.matrix(location)
    for (i in 1:nrow(df_VE)) {
      tmp_protein <- df_impu[i,]
      for (j in 1:length(tmp_protein)) {
        location_tmp[location_corr$row[j], location_corr$column[j]] <- as.numeric(tmp_protein[j])
      }
      location_tmp <- apply(location_tmp,1, as.numeric)
      location_tmp[is.na(location_tmp)] <- min(tmp_protein, na.rm = TRUE)
      out <- svd(location_tmp)
      df_VE$VE[i] <- out$d[1]/sum(out$d)
      location_tmp <- as.matrix(location)
    }

    df_VE <- df_VE[order(df_VE$VE, decreasing = TRUE),]
    df_VE$Index <- 1:nrow(df_VE)
    VE_med <-median(df_VE$VE, na.rm=TRUE)

    VE_plot <- ggplot(df_VE, aes(x = Index, y = VE)) +
      geom_point() + xlab('Protein Rank') + ylab('Variance Explained (%)') +
      geom_hline(yintercept=VE_med, linetype="dashed", color = "red", size=2) +
      theme(legend.position="right",panel.grid.minor=element_blank(),
            panel.background=element_rect(fill=NA, colour=NA),
            panel.border=element_rect(size=1,fill=NA,colour="black")) +
      theme(axis.text=element_text(size=12, family="Helvetica", colour="black"),
            axis.title=element_text(size=12, family="Helvetica",colour="black"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.text=element_text(size=12, family="Helvetica",colour="black"),
            legend.title=element_text(size=12, family="Helvetica", colour="black"))

    outputs$VE_results <- df_VE
    outputs$VE_plot <- VE_plot
    return(outputs)
  })
  return(outputs)
}

#' @export
VE_filtering <- function(df, setting, value = 0.2){
  if(setting == 'med') protein_filtered <- df$protein_name[df$VE > median(df$VE)]
  if(setting == 'customized') protein_filtered <- df$protein_name[df$VE > value]
  return(protein_filtered)
}

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
