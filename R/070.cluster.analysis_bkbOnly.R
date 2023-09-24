#' Cluster analysis with normalised backbone measurements and the complete dataset
#' 
#' This function has been designed to perform cluster analysis for the normalised backbone measurements and the complete dataset which includes the normalised backbone measurements and the imputed well-specific markers.
#' 
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param bkb.v A vector of the names of the backbone markers (MUST match to the names in the FCS file).
#' 
#' @author Hsiao-Chi Liao
#' 
#' @import ggplot2
#' @importFrom uwot umap
#' @importFrom Rphenograph Rphenograph
#' @importFrom igraph membership
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#' 
#' @return Updating the metadata for cells in the fcs_metadata_df.rds file, adding the information of the biological clusters from the clean and complete dataset. Visualising the result with the scatter plots.
#' 
cluster.analysis.bkbOnly <-
function(
    paths, 
    bkb.v
    ){
      # RDSpath="/Users/chelsea/Desktop/PhD/ProteinExpression/InFlow/Programming/infinity/R_CODE/72.mk_pkg/data/rds/"
      # Graphpath="/Users/chelsea/Desktop/PhD/ProteinExpression/InFlow/Programming/infinity/R_CODE/72.mk_pkg/data/graph/ClusterAnalysis/"
      # bkb.v <- c("FSC-H", "FSC-W", "SSC-H", "SSC-W",
      #            "CD69-CD301b", "MHCII", "CD4", "CD44", "CD8",
      #            "CD11c", "CD11b", "F480", "Ly6C", "Lineage", "CD45a488",
      #            "CD24", "CD103")
      # yvar="Legend"
      # control.wells <- c("P1_A01", "P2_A01", "P3_A01",
      #                    "P3_F04", "P3_F05", "P3_F06", "P3_F07", "P3_F08", "P3_F09",
      #                    "P3_F10", "P3_F11", "P3_F12", "P3_G01", "P3_G02")
      
      normalised.bkb <- readRDS(file = file.path(paths["downstream"], "bkc.adj.bkb_logScale_mt.rds")) 
      metadata.cell <- readRDS(file = file.path(paths["intermediary"], "fcs_metadata_df.rds"))
      names(bkb.v) <- make.names(bkb.v)
      
      message(paste0("Cluster analysis for normalised backbone measurements..."))
      message("Clustering with normalised backbones")
      bkb.dat <- normalised.bkb
      message("Running UMAP...")
      a <- Sys.time()
      umap.bkb <- umap(bkb.dat,
                       n_neighbors = 15, min_dist = 0.2, metric = "euclidean", n_epochs = 2000)
      b <- Sys.time()
      b-a #12.46998 mins
      saveRDS(umap.bkb, file = file.path(paths["downstream"], paste0("ClusterAnalysis_", "umap_",length(bkb.v),"bkb.rds")))
      
      message("Running Phenograph...")
      a <- Sys.time()
      phenog.bkb <- Rphenograph(bkb.dat, 
                                k = 50)  #knn_fun = "hnsw", 
      b <- Sys.time()
      b-a #5.802254 mins
      saveRDS(phenog.bkb, file = file.path(paths["intermediary"], paste0("ClusterAnalysis_", "phenog_",length(bkb.v),"bkb.rds")))
      
      ##
      metadata.cell[,paste0("GP.denoised.bkb")] <- as.factor(membership(phenog.bkb[[2]]))
      
      message("\tVisualising clusters...")
      {
        # pheno color
        # col.v <- c("#6E9A8A", "#6DE5DC", "#9237E9", "#B0E7BE", "#8761DB",
        #            "#71EA4A", "#9ACD49", "#E8813B", "#DEEB45", "#CED0EA",
        #            "salmon4", "#D8868B", "#EAE6C5", "green", "#E2AC7C",
        #            "#E050D7", "#887384", "#73AABE", "#DFD48C", "blue",
        #            "deeppink3", "mediumorchid4", "#ED574E", "#E5EB77", "#CEABDF",
        #            "#CAE7E6", "#6CAF76", "#5FE9B8", "#E0B341", "#B5E893",
        #            "#777FD5", "#64D4EB", "yellow")    ##33
        ## bkb
        {
          graph.dat <- data.frame(umap.bkb, metadata.cell[paste0("GP.denoised.bkb")])
          colnames(graph.dat) <- c("UMAP1", "UMAP2", "Cluster")
          
          n <- length(unique(graph.dat[,3]))
          qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',] #up to 74
          col.coeff <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:n]
          
          jpeg(file = file.path(paths["graph"], paste0("/ClusterStructure_UMAP_", length(bkb.v), "bkb_colPhenog_600x750.jpeg")), height = 600, width = 750, res = 80)
          p <- ggplot(graph.dat, aes(UMAP1, UMAP2, colour=Cluster)) + geom_point(size = 0.05, alpha = 0.25) + 
            labs(titles = paste0("Normalised ", length(bkb.v), "bkb. \n(", nrow(graph.dat), " cells)")) +
            guides(colour = guide_legend(override.aes = list(size=8))) + theme_bw() + 
            theme(text = element_text(size=22)) + scale_color_manual(values=col.coeff[1:n])
          print(p)
          dev.off() 
        }
      }
      ###
      head(metadata.cell)
      saveRDS(metadata.cell, file = file.path(paths["downstream"], "fcs_metadata_df.rds"))
      
      message("\tCompleted!")
    }
