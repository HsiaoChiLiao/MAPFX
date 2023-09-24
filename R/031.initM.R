#' Initial biological clusters
#' 
#' Generating the M matrix for removing well effects.
#' 
#' This function has been designed to find initial biological clusters with centred transformed data.
#' 
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param bkb.v A vector of the names of the backbone markers (MUST match to the names in the FCS file).
#' @param trans.dat Method for transforming raw intensities (default: Logicle).
#' 
#' @author Hsiao-Chi Liao
#' 
#' @import ggplot2
#' @importFrom flowCore logicleTransform
#' @importFrom uwot umap
#' @importFrom Rphenograph Rphenograph
#' @importFrom igraph membership
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#' 
#' @return Updating the metadata for cells in the fcs_metadata_df.rds file, adding the information of the initial biological clusters. Visualising the result with the scatter plots.
#' 
initM <-
function(
    paths,
    bkb.v,
    trans.dat="cent.lgc" #may give "cent.log.bkc" in the future
    ){
      # RDSpath="/Users/chelsea/Desktop/PhD/ProteinExpression/InFlow/Programming/infinity/R_CODE/72.mk_pkg/data/rds/"
      # Graphpath="/Users/chelsea/Desktop/PhD/ProteinExpression/InFlow/Programming/infinity/R_CODE/72.mk_pkg/data/graph/rmWellEffect/"
      # trans.dat <- "cent.lgc"
      # 
      # bkb.v <- c("FSC-H", "FSC-W", "SSC-H", "SSC-W",
      #            "CD69-CD301b", "MHCII", "CD4", "CD44", "CD8",
      #            "CD11c", "CD11b", "F480", "Ly6C", "Lineage", "CD45a488",
      #            "CD24", "CD103")
      
      rawInten <- readRDS(file = file.path(paths["intermediary"], "fcs_rawInten_mt.rds"))
      metadata.cell <- readRDS(file = file.path(paths["intermediary"], "fcs_metadata_df.rds"))
      sel.bkb.raw <- rawInten[,match(bkb.v, colnames(rawInten))]
      
      if(trans.dat == "cent.lgc"){
        ## lgc trans - borrowing inflow people
        {
          ##bkb only
          chans <- bkb.v
          xp <- sel.bkb.raw
          message("\tForming logicle functions...")
          {
            transforms_chan <- setNames(
              lapply(
                chans,
                function(x){
                  data <- xp[,x] ## all cells (together do transformation)
                  
                  t <- max(data)
                  m <- 4.5
                  
                  q <- 0.05
                  r <- .Machine$double.eps + quantile(data, q)
                  w <- max((m-log10(t/abs(r)))/2,0.1)
                  w <- min(w,m/2)
                  
                  a <- 0
                  logicleTransform(w=w,t=t,m=m,a=a) ##Just use summary() to retrive the parameters
                }
              ),
              chans
            )
            saveRDS(transforms_chan, file=file.path(paths["intermediary"],"transform.funcs_bkb.Rds")) 
          }
          
          message("\tLogicle transforming raw intensity...")
          {
            sel.bkb.lgc <- sel.bkb.raw
            for(chan in chans){
              sel.bkb.lgc[,chan] <- transforms_chan[[chan]](xp[,chan])
            }
            head(sel.bkb.lgc)
          }
          
          message("\tCentring logicle transformed intensities...")
          {
            ## shifted by well-wise mean
            well.v <- unique(metadata.cell$Well.lab)
            shifted.ls <- list()
            for(i in 1:length(well.v)){
              dat.here <- sel.bkb.lgc[which(metadata.cell$Well.lab == well.v[i]),]
              dat.meanSF <- apply(dat.here, 2, function(x) (x-mean(x)))
              shifted.ls[[i]] <- dat.meanSF
            }
            
            sel.bkb.cent.lgc <- do.call(rbind, shifted.ls)
            saveRDS(sel.bkb.cent.lgc, file = file.path(paths["intermediary"], "sel.bkb.cent.lgc.rds"))
          }
        }
        message("\tCentred logicle backbone data... Obtained!")
      }
      
      message("\tUMAP with backbones...")
      a <- Sys.time()
      umap.bkb.cent.lgc <- umap(sel.bkb.cent.lgc, 
                                n_neighbors = 15, min_dist = 0.2, metric = "euclidean", n_epochs = 2000) #, n_neighbors = 15, min_dist = 0.2, metric = "euclidean", n_epochs = 2000
      b <- Sys.time()
      print(b-a) #Time difference of 22.56658 mins
      saveRDS(umap.bkb.cent.lgc, file = file.path(paths["intermediary"], "initM_umap.bkb.cent.lgc.rds"))
      
      message("\tDeriving initial clusters with PhenoGraph (forming the M matrix)...")
      a <- Sys.time()
      phenog.bkb.cent.lgc <- Rphenograph(sel.bkb.cent.lgc, 
                                         k = 50) #k = 50, knn_fun = "hnsw" (knn_fun: doesn't work here...); RcppHNSW_0.4.1        Rphenograph_0.99.1
      b <- Sys.time()
      print(b-a) #Time difference of 16.25608 mins
      saveRDS(phenog.bkb.cent.lgc, file = file.path(paths["intermediary"], "initM_phenog.bkb.cent.lgc.rds"))
      
      metadata.cell$init.M <- as.factor(membership(phenog.bkb.cent.lgc[[2]]))
      saveRDS(metadata.cell, file = file.path(paths["intermediary"], "fcs_metadata_df.rds"))
      
      message("\tVisualising clusters...")
      {
        # pheno color
        col.v <- c("#6E9A8A", "#6DE5DC", "#9237E9", "#B0E7BE", "#8761DB",
                   "#71EA4A", "#9ACD49", "#E8813B", "#DEEB45", "#CED0EA",
                   "salmon4", "#D8868B", "#EAE6C5", "green", "#E2AC7C",
                   "#E050D7", "#887384", "#73AABE", "#DFD48C", "blue",
                   "deeppink3", "mediumorchid4", "#ED574E", "#E5EB77", "#CEABDF",
                   "#CAE7E6", "#6CAF76", "#5FE9B8", "#E0B341", "#B5E893",
                   "#777FD5", "#64D4EB", "yellow")    ##33
        graph.dat <- data.frame(umap.bkb.cent.lgc, metadata.cell$init.M)
        colnames(graph.dat) <- c("UMAP1", "UMAP2", "Pheno.gp")
        
        n <- length(unique(graph.dat$Pheno.gp))
        qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
        col.coeff <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:n]
        
        jpeg(file = file.path(paths["graph"], "/", paste0(trans.dat, "_UMAP_colPhenog_600x750.jpeg")), height = 600, width = 750, res = 80)
        p <- ggplot(graph.dat, aes(UMAP1, UMAP2, colour=Pheno.gp)) + geom_point(size = 0.05, alpha = 0.25) + 
          labs(titles = paste0("init. bio. (M matrix) from ", trans.dat," bkb. \n(", nrow(graph.dat)," cells)")) +
          guides(colour = guide_legend(override.aes = list(size=8))) + theme_bw() + 
          theme(text = element_text(size=22)) + scale_color_manual(values=col.coeff[1:n])
        print(p)
        dev.off()
      }
      message("\tCompleted!")
    }
