#' Removing well effect from the data
#' 
#' This function has been designed to remove the unwanted effects (well effects) from the background corrected measurements.
#' 
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param visualisation Whether to visualise the results by using heatmap (default=TRUE).
#' 
#' @author Hsiao-Chi Liao
#' 
#' @import ggplot2
#' @import ggrepel
#' @import ComplexHeatmap
#' @importFrom Rfast lmfit colsums
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#' @importFrom circlize colorRamp2
#' 
#' @return Generating the calibrated measurements and save to bkc.adj.bkb_logScale_mt.rds (on log scale) and bkc.adj.bkb_linearScale_mt.rds (on linear scale). Visualising the result with the heatmaps.
#' 
rmWellEffect <-
function(
    paths, visualisation = TRUE
    ){
      # RDSpath="/Users/chelsea/Desktop/PhD/ProteinExpression/InFlow/Programming/infinity/R_CODE/72.mk_pkg/data/rds/"
      # visualisation = TRUE
      # Graphpath="/Users/chelsea/Desktop/PhD/ProteinExpression/InFlow/Programming/infinity/R_CODE/72.mk_pkg/data/graph/adj/"
      bkc.bkb <- readRDS(file = file.path(paths["intermediary"], "medpara_bkc.bkb_no.bkcPhy_mt.rds")) #no. of selected bkb
      metadata.cell <- readRDS(file = file.path(paths["intermediary"], "fcs_metadata_df.rds"))
      
      ### Rfast to estimate coef for later adj. ###
      {
        message("\tEstimating coefficients for removing well effect (Rfast - pre.adj)...")
        a <- Sys.time()
        {
          metadata.cell[,c("Plate","Column","Row","init.M")] <- lapply(metadata.cell[,c("Plate","Column","Row","init.M")], factor)
          contrasts(metadata.cell[,"Plate"]) <- contr.sum(3, contrasts=TRUE)
          contrasts(metadata.cell[,"Column"]) <- contr.sum(12, contrasts=TRUE)
          contrasts(metadata.cell[,"Row"]) <- contr.sum(8, contrasts=TRUE)
          contrasts(metadata.cell[,"init.M"]) <- contr.sum(length(unique(metadata.cell$init.M)), contrasts=TRUE)
          
          xx <- as.matrix(model.matrix(~ Plate + Column + Row + init.M, data = metadata.cell))
          
          lm.results.ls <- list()
          ln.sig.v <- c()
          res.df <- bkc.bkb
          for(i in 1:ncol(bkc.bkb)){
            yy.val <- as.matrix(bkc.bkb[,i])
            
            lm.result <- tryCatch(lmfit(x = xx, y = log(yy.val+0.000001)), error=function(err) NA)
            lm.results.ls[[i]] <- lm.result
            res <- lm.result$residuals
            ln.sig.v[i] <- sqrt(sum( (res - mean(res) )^2 )/(length(res)-length(lm.results.ls[[1]]$be))) #coef used
            # ln.sig.v[i] <- sd(lm.result$residuals)
            
            res.df[,i] <- res
            
            message(paste0("Processing bkb: ", i))
          }
          
          names(lm.results.ls) <- colnames(bkc.bkb)
          
          quant.bio.pcr.mt <- do.call(cbind, lapply(lm.results.ls, '[[', 1))
          colnames(quant.bio.pcr.mt) <- names(lm.results.ls)
          quant.bio.pcr.mt <- rbind(quant.bio.pcr.mt, ln.sig = ln.sig.v)
          save(lm.results.ls, quant.bio.pcr.mt, file = file.path(paths["intermediary"], "rmWellEffect_log.normal.reg_coef.mdl_pre.RData"))
          save(res.df, file = file.path(paths["intermediary"], "rmWellEffect_log.normal.reg_residuals.RData"))
        }
        b <- Sys.time()
        b-a
        message("\tEstimation completed!")
        #Time difference of 13.76776 secs
      }
      ### adjusting data - subtracting unwanted effect ###
      {
        message("\tRemoving well effect for backbone markers...")
        adj.data = log.adj.data <- bkc.bkb
        a <- Sys.time()
        for(i in 1:ncol(bkc.bkb)){
          unwanted.eff <- (xx[,(2:21)] %*% quant.bio.pcr.mt[(2:21), i])
          
          log.adj.data[,i] <- (log(bkc.bkb[,i]+0.000001) - unwanted.eff) #as the coef was estimated from data on log scale
          
          adj.data[,i] <- (exp(log.adj.data[,i])-0.000001)
        }
        b <- Sys.time()
        b-a
        saveRDS(log.adj.data, file = file.path(paths["downstream"], "bkc.adj.bkb_logScale_mt.rds"))
        saveRDS(adj.data, file = file.path(paths["downstream"], "bkc.adj.bkb_linearScale_mt.rds"))
        message("\tAdjustment completed!")
        #Time difference of 1.667447 secs
      }
      ### Rfast to estimate coef for later adj. ###
      {
        message("\tExamining the existence of well effect in the adjusted data (Rfast - post.adj)...")
        a <- Sys.time()
        {
          metadata.cell[,c("Plate","Column","Row","init.M")] <- lapply(metadata.cell[,c("Plate","Column","Row","init.M")], factor)
          contrasts(metadata.cell[,"Plate"]) <- contr.sum(3, contrasts=TRUE)
          contrasts(metadata.cell[,"Column"]) <- contr.sum(12, contrasts=TRUE)
          contrasts(metadata.cell[,"Row"]) <- contr.sum(8, contrasts=TRUE)
          contrasts(metadata.cell[,"init.M"]) <- contr.sum(length(unique(metadata.cell$init.M)), contrasts=TRUE)
          
          xx <- as.matrix(model.matrix(~ Plate + Column + Row + init.M, data = metadata.cell))
          
          adj.lm.results.ls <- list()
          # glm.results.ls <- list()
          ln.sig.v <- c()
          for(i in 1:ncol(log.adj.data)){
            yy.val <- as.matrix(log.adj.data[,i])
            
            lm.result <- tryCatch(lmfit(x = xx, y = yy.val), error=function(err) NA)
            adj.lm.results.ls[[i]] <- lm.result
            res <- lm.result$residuals
            ln.sig.v[i] <- sqrt(sum( (res - mean(res) )^2 )/(length(res)-length(lm.results.ls[[1]]$be))) #coef used
            # ln.sig.v[i] <- sd(lm.result$residuals)
            
            message(paste0("Processing bkb: ", i))
          }
          
          names(adj.lm.results.ls) <- colnames(log.adj.data)
          
          adj.quant.bio.pcr.mt <- do.call(cbind, lapply(adj.lm.results.ls, '[[', 1))
          colnames(adj.quant.bio.pcr.mt) <- names(adj.lm.results.ls)
          adj.quant.bio.pcr.mt <- rbind(adj.quant.bio.pcr.mt, ln.sig = ln.sig.v)
          save(lm.results.ls, quant.bio.pcr.mt, file = file.path(paths["intermediary"], "rmWellEffect_log.normal.reg_coef.mdl_post.RData"))
        }
        b <- Sys.time()
        b-a
        #Time difference of 10.95444 secs
      }
      
      ### rmWellEffect_bkb_visualisation (heatmap) - bio and pcr effects ###
      if(visualisation == TRUE){
        {
          #adjust y the num. one
          #no intercept and sigma
          in.dat <- list(quant.bio.pcr.mt[-c(1,nrow(quant.bio.pcr.mt)),], adj.quant.bio.pcr.mt[-c(1,nrow(adj.quant.bio.pcr.mt)),])
          in.dat <- list(rbind(in.dat[[1]][1:2,], Plate3 = -colsums(in.dat[[1]][1:2,]),
                               in.dat[[1]][3:13,], Column12 = -colsums(in.dat[[1]][3:13,]),
                               in.dat[[1]][14:20,], Row8 = -colsums(in.dat[[1]][14:20,]),
                               in.dat[[1]][21:nrow(in.dat[[1]]),], init.Mlast = -colsums(in.dat[[1]][21:nrow(in.dat[[1]]),])),
                         rbind(in.dat[[2]][1:2,], Plate3 = -colsums(in.dat[[2]][1:2,]),
                               in.dat[[2]][3:13,], Column12 = -colsums(in.dat[[2]][3:13,]),
                               in.dat[[2]][14:20,], Row8 = -colsums(in.dat[[2]][14:20,]),
                               in.dat[[2]][21:nrow(in.dat[[2]]),], init.Mlast = -colsums(in.dat[[2]][21:nrow(in.dat[[2]]),]))
          )
          rownames(in.dat[[1]])[length(rownames(in.dat[[1]]))] <- paste0("init.M", (length(24:nrow(in.dat[[1]]))))
          rownames(in.dat[[2]])[length(rownames(in.dat[[2]]))] <- paste0("init.M", (length(24:nrow(in.dat[[2]]))))
          
          ## consistent colour key
          bio.lim <- c(min(in.dat[[1]][-(1:23),], in.dat[[2]][-(1:23),]), max(in.dat[[1]][-(1:23),], in.dat[[2]][-(1:23),]))
          well.lim <- c(min(in.dat[[1]][(1:23),], in.dat[[2]][(1:23),]), max(in.dat[[1]][(1:23),], in.dat[[2]][(1:23),]))
          
          pre.post <- c("pre-", "post-")
          bio.sys <- c("bio.", "well.")
          for(i in 1:2){ #pre & post
            for(j in 1:2){ #bio & batch
              ### heatmap
              dat0 <- in.dat[[i]]
              
              if(j == 1){ #biology
                bio.dat <- dat0[-(1:23),]
                # metadata
                meta_mat = data.frame(coeff=rownames(bio.dat))
                meta_mat$Coef = factor(paste0('init.M', 1:nrow(meta_mat)), 
                                       levels = paste0('init.M', 1:nrow(meta_mat)))
                rownames(meta_mat) <- meta_mat$Coef
                
                n <- nrow(meta_mat)
                qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
                col.coeff = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:n]
                names(col.coeff) <- meta_mat$Coef
                row_ha = rowAnnotation(Effect=meta_mat[,'Coef'],
                                       col=list(Effect=col.coeff), show_legend = F)
                #elements in `col` should be named vectors.
                
                # color
                coeff_col_func <- colorRamp2(c(bio.lim[1], 0, bio.lim[2]), c("blue", "gray0", "darkorange1"))
                
                #not cluster rows/ cols
                jpeg(file.path(paths["graph"], paste0("Heatmap_coeffs_", pre.post[i], bio.sys[j], "_2000x1500.jpeg")), height = 2000, width = 1500, res = 250)
                draw(Heatmap(bio.dat, name=paste0("MLE"), #Rfast
                             column_title = paste0(pre.post[i], "well.adj."),
                             row_title = bio.sys[j],
                             right_annotation = row_ha,
                             col = coeff_col_func,
                             column_dend_height = unit(1,"cm"),  # tree size
                             cluster_rows = F,
                             cluster_columns = F,
                             clustering_method_rows = "average",
                             clustering_method_columns = "average",
                             show_column_names = T,
                             show_row_names = T)
                )
                dev.off()
              }else{
                well.dat <- dat0[(1:23),]
                # metadata
                meta_mat = data.frame(coeff=rownames(well.dat))
                meta_mat$Coef = factor(c(rep('logn.mu.Plate', 3), rep('logn.mu.Column', 12), rep('logn.mu.Row', 8)), 
                                       levels = c('logn.mu.Plate', 'logn.mu.Column', 'logn.mu.Row'))
                
                # col. anno for heatmap (be careful of the position of each variable)
                #display.brewer.pal(n = 8, name = 'Dark2')
                col.coeff <- c(brewer.pal(9,"Set1")[c(1,3,2)], brewer.pal(8,"Dark2")[c(4,1,3)])
                row_ha = rowAnnotation(Effect=meta_mat[,'Coef'],
                                       col=list(Effect=c(
                                         'logn.mu.Plate' = col.coeff[4],
                                         'logn.mu.Column' = col.coeff[5], 
                                         'logn.mu.Row' = col.coeff[6])))
                # color
                coeff_col_func <- colorRamp2(c(well.lim[1], 0, well.lim[2]), c("blue", "gray0", "darkorange1"))
                
                #not cluster rows/ cols
                jpeg(file.path(paths["graph"], paste0("Heatmap_coeffs_", pre.post[i], bio.sys[j], "_2000x1500.jpeg")), height = 2000, width = 1500, res = 250)
                draw(Heatmap(well.dat, name=paste0("MLE"), #Rfast
                             column_title = paste0(pre.post[i], "well.adj."),
                             row_title = bio.sys[j],
                             right_annotation = row_ha,
                             col = coeff_col_func,
                             column_dend_height = unit(1,"cm"),  # tree size
                             cluster_rows = F,
                             cluster_columns = F,
                             clustering_method_rows = "average",
                             clustering_method_columns = "average",
                             show_column_names = T,
                             show_row_names = T)
                )
                dev.off()
              }
            }
          }
        }
      }
    }
