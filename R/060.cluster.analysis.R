#' Cluster analysis with normalised backbone measurements and the complete dataset
#' 
#' This function has been designed to perform cluster analysis for the normalised backbone measurements and the complete dataset which includes the normalised backbone measurements and the imputed well-specific markers.
#' 
#' @param paths a vector of characters of paths to store intput, intermediary results, outputs...
#' @param bkb.v a vector of the names of the backbone markers (MUST match to the names in the FCS file).
#' @param yvar the name of the well-specific marker in the FCS files (e.g., "Legend").
#' @param control.wells the well label of the control wells, including the autofluorescence and the isotype controls (format: plate_well, e.g., P1_A01)
#' @param plots logical; if TRUE (default), produce an UMAP embedding plot from the normalised backbone markers and the imputed infinity markers to visualise the structure of the biological clusters.
#' 
#' @author Hsiao-Chi Liao
#' 
#' @import ggplot2
#' @importFrom uwot umap
#' @importFrom iCellR Rphenograph
#' @importFrom igraph membership
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#' 
#' @importFrom grDevices dev.off jpeg
#' @importFrom graphics abline
#' @importFrom methods is
#' @importFrom stats as.formula contr.sum contrasts<- dexp dnorm median model.matrix optim pexp pnorm quantile sd setNames
#' @importFrom utils head read.csv
#' 
#' @return Metadata for cells with group labels from the cluster analysis
#' 
#' @details
#' Updating the metadata for cells in the fcs_metadata_df.rds file, adding the information of the biological clusters from the clean and complete dataset, and visualising the result with the scatter plots in the output directory.
#' 
cluster.analysis <-
function(
    paths,
    bkb.v,
    yvar="Legend",
    control.wells,
    plots=TRUE
    ){
    
    #binding the variable locally to the function
    UMAP1 <- UMAP2 <- Cluster <- NULL
    
    preds <- readRDS(file = file.path(paths["downstream"], "predictions.Rds"))
    metadata.cell <- readRDS(file = file.path(paths["intermediary"], "fcs_metadata_df.rds"))
    names(bkb.v) <- make.names(bkb.v)
    
    for(i in seq_along(preds)){
    impu.mt <- preds[[i]]
    
    if(i == 1){ ## bkb only need to do once
    message("Clustering with normalised backbones")
    bkb.dat <- impu.mt[, match(names(bkb.v), colnames(impu.mt)),drop=FALSE]
    message("Running UMAP...")
    umap.bkb <- umap(
    bkb.dat,
    n_neighbors = 15, min_dist = 0.2, metric = "euclidean", n_epochs = 2000)
    saveRDS(umap.bkb, file = file.path(paths["downstream"], paste0("ClusterAnalysis_", "umap_",length(bkb.v),"bkb.rds")))
    
    message("Running Phenograph...")
    phenog.bkb <- Rphenograph(bkb.dat, k = 50)  #knn_fun = "hnsw", 
    saveRDS(phenog.bkb, file = file.path(paths["intermediary"], paste0("ClusterAnalysis_", "phenog_",length(bkb.v),"bkb.rds")))
    }
    
    ###
    message("Clustering with normalised backbones + imputed infinity markers (", names(preds)[i], ")") #now: use them all - may just pick good ones in the future
    complete.dat <- impu.mt[, -match(c(yvar, control.wells), colnames(impu.mt)),drop=FALSE]
    message("Running UMAP...")
    umap.bkb.impuInf <- umap(
    complete.dat,
    n_neighbors = 15, min_dist = 0.2, metric = "euclidean", n_epochs = 2000)
    saveRDS(umap.bkb.impuInf, 
    file = file.path(paths["downstream"], paste0("ClusterAnalysis_", names(preds)[i], "_umap_", 
    length(bkb.v),"bkb.", (ncol(complete.dat)-length(bkb.v)),"impuInf.rds")))
    
    message("Running Phenograph...")
    phenog.bkb.impuInf <- Rphenograph(complete.dat, k = 50)  #knn_fun = "hnsw", 
    saveRDS(phenog.bkb.impuInf, 
    file = file.path(paths["intermediary"], paste0("ClusterAnalysis_", names(preds)[i], "_phenog_",
    length(bkb.v),"bkb.", (ncol(complete.dat)-length(bkb.v)),"impuInf.rds")))
    
    if(i == 1){
    metadata.cell[which(metadata.cell$train_set == 0),paste0("GP.denoised.bkb")] <- as.factor(membership(phenog.bkb[[2]]))
    }
    metadata.cell[which(metadata.cell$train_set == 0),paste0("GP.denoised.bkb.impuInf_", names(preds)[i])] <- as.factor(membership(phenog.bkb.impuInf[[2]]))
    
    if(plots == TRUE){
    message("\tVisualising clusters...")
    {
    ## bkb
    {
    if(i == 1){ ## would not be affected by imputation method
    graph.dat <- data.frame(umap.bkb, metadata.cell[which(metadata.cell$train_set == 0), paste0("GP.denoised.bkb")])
    colnames(graph.dat) <- c("UMAP1", "UMAP2", "Cluster")
    
    n <- length(unique(graph.dat[,3]))
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',] #up to 74
    col.coeff <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[seq_len(n)]
    
    FN <- paste0("ClusterStructure_UMAP_", length(bkb.v), "bkb_colPhenog_600x750.jpeg")
    jpeg(filename = file.path(paths["graph"], FN), height = 600, width = 750, res = 80)
    p <- ggplot(
    data=graph.dat, mapping = aes(x=UMAP1, y=UMAP2, colour=Cluster)) + geom_point(size = 0.05, alpha = 0.25) + 
    labs(titles = paste0("Normalised ", length(bkb.v), "bkb. \n(", nrow(graph.dat), " cells)")) +
    guides(colour = guide_legend(override.aes = list(size=8))) + theme_bw() + 
    theme(text = element_text(size=22)) + scale_color_manual(values=col.coeff[seq_len(n)])
    plot(p)
    dev.off() 
    }
    }
    ## bkb+impuInf
    {
    graph.dat <- data.frame(umap.bkb.impuInf, metadata.cell[which(metadata.cell$train_set == 0), paste0("GP.denoised.bkb.impuInf_",names(preds)[i])])
    colnames(graph.dat) <- c("UMAP1", "UMAP2", "Cluster")
    
    n <- length(unique(graph.dat[,3]))
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
    col.coeff <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[seq_len(n)]
    
    FN <- paste0("ClusterStructure_UMAP_", names(preds)[i], length(bkb.v), "bkb",
    (ncol(complete.dat)-length(bkb.v)),"impuInf_colPhenog_600x750.jpeg")
    jpeg(filename = file.path(paths["graph"], FN), height = 600, width = 750, res = 80)
    p <- ggplot(
    data=graph.dat, mapping = aes(x=UMAP1, y=UMAP2, colour=Cluster)) + geom_point(size = 0.05, alpha = 0.25) + 
    labs(titles = paste0(names(preds)[i], " ", length(bkb.v), "bkb. ", (ncol(complete.dat)-length(bkb.v)),"impuInf \n(", nrow(graph.dat), " cells)")) +
    guides(colour = guide_legend(override.aes = list(size=8))) + theme_bw() + 
    theme(text = element_text(size=22)) + scale_color_manual(values=col.coeff[seq_len(n)])
    plot(p)
    dev.off()
    }
    }
    }
    }
    ###
    saveRDS(metadata.cell, file = file.path(paths["downstream"], "fcs_metadata_df.rds"))
    
    message("\tCompleted!")
    
    return(metadata.cell)
    }
    