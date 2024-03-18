#' Cluster analysis with normalised backbone measurements
#' 
#' This function has been designed to perform cluster analysis for the normalised backbone measurements.
#' 
#' @param paths a vector of characters of paths to store intput, intermediary results, outputs...
#' @param bkb.v a vector of the names of the backbone markers (MUST match to the names in the FCS file).
#' @param plots logical; if TRUE (default), produce an UMAP embedding plot from the normalised backbone markers to visualise the structure of the biological clusters.
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
#' @return Updating the metadata for cells in the fcs_metadata_df.rds file, adding the information of the biological clusters from the clean and complete dataset. Visualising the result with the scatter plots.
#' 
cluster.analysis.bkbOnly <-
function(
    paths, 
    bkb.v, #bkbs involved in the plotting
    plots=TRUE
    ){
    
    #binding the variable locally to the function
    UMAP1 <- UMAP2 <- Cluster <- NULL
    
    normalised.bkb <- readRDS(file = file.path(paths["downstream"], "bkc.adj.bkb_logScale_mt.rds")) 
    metadata.cell <- readRDS(file = file.path(paths["intermediary"], "fcs_metadata_df.rds"))
    names(bkb.v) <- make.names(bkb.v)
    
    message("Cluster analysis for normalised backbone measurements...")
    message("Clustering with normalised backbones")
    bkb.dat <- normalised.bkb[,match(bkb.v,colnames(normalised.bkb))]
    message("Running UMAP...")
    a <- Sys.time()
    umap.bkb <- umap(bkb.dat, n_neighbors = 15, min_dist = 0.2, metric = "euclidean", n_epochs = 2000)
    b <- Sys.time()
    b-a #12.46998 mins
    saveRDS(umap.bkb, file = file.path(paths["downstream"], paste0("ClusterAnalysis_", "umap_",length(bkb.v),"bkb.rds")))
    
    message("Running Phenograph...")
    a <- Sys.time()
    phenog.bkb <- Rphenograph(bkb.dat, k = 50)  #knn_fun = "hnsw", 
    b <- Sys.time()
    b-a #5.802254 mins
    saveRDS(phenog.bkb, file = file.path(paths["intermediary"], paste0("ClusterAnalysis_", "phenog_",length(bkb.v),"bkb.rds")))
    
    ##
    metadata.cell[,paste0("GP.denoised.bkb")] <- as.factor(membership(phenog.bkb[[2]]))
    
    if(plots == TRUE){
    message("\tVisualising clusters...")
    {
    ## bkb
    {
    graph.dat <- data.frame(umap.bkb, metadata.cell[paste0("GP.denoised.bkb")])
    colnames(graph.dat) <- c("UMAP1", "UMAP2", "Cluster")
    
    n <- length(unique(graph.dat[,3]))
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',] #up to 74
    col.coeff <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[seq_len(n)]
    
    jpeg(filename = file.path(paths["graph"], paste0("/ClusterStructure_UMAP_", length(bkb.v), "bkb_colPhenog_600x750.jpeg")), height = 600, width = 750, res = 80)
    p <- ggplot(
    data=graph.dat, mapping = aes(x=UMAP1, y=UMAP2, colour=Cluster)) + geom_point(size = 0.05, alpha = 0.25) + 
    labs(titles = paste0("Normalised ", length(bkb.v), "bkb. \n(", nrow(graph.dat), " cells)")) +
    guides(colour = guide_legend(override.aes = list(size=8))) + theme_bw() + 
    theme(text = element_text(size=22)) + scale_color_manual(values=col.coeff[seq_len(n)])
    plot(p)
    dev.off() 
    }
    }
    }
    
    ###
    head(metadata.cell)
    saveRDS(metadata.cell, file = file.path(paths["downstream"], "fcs_metadata_df.rds"))
    
    message("\tCompleted!")
    }
