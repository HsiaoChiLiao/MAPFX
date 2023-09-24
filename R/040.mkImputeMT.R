#' Making the input for imputation
#' 
#' This function has been designed to combine the normalised backbone measurements and the normalised PE markers for later imputation.
#' 
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' 
#' @author Hsiao-Chi Liao
#' 
#' @return Generating the combined data and save to impu.input_log.mt.rds (on log scale).
#' 
mkImputeMT <-
function(
    paths
    ){
      # RDSpath="/Users/chelsea/Desktop/PhD/ProteinExpression/InFlow/Programming/infinity/R_CODE/72.mk_pkg/data/rds/"
      bkc.adj.bkb <- readRDS(file = file.path(paths["downstream"], "bkc.adj.bkb_logScale_mt.rds")) #no. of selected bkb
      bkc.pe <- readRDS(file = file.path(paths["intermediary"], "bkc.pe_mt.rds"))
      metadata.cell <- readRDS(file = file.path(paths["intermediary"], "fcs_metadata_df.rds"))
      
      impu.input_log.mt <- cbind(bkc.adj.bkb, log(bkc.pe))
      rownames(impu.input_log.mt) <- metadata.cell$NO.in.all
      
      saveRDS(impu.input_log.mt, file = file.path(paths["intermediary"], "impu.input_log.mt.rds"))
    }
