#' Making the input for imputation
#' 
#' This function has been designed to combine the normalised backbone measurements and the normalised PE markers for later imputation.
#' 
#' @param paths a vector of characters of paths to store intput, intermediary results, outputs...
#' 
#' @author Hsiao-Chi Liao
#' 
#' @importFrom grDevices dev.off jpeg
#' @importFrom graphics abline
#' @importFrom methods is
#' @importFrom stats as.formula contr.sum contrasts<- dexp dnorm median model.matrix optim pexp pnorm quantile sd setNames
#' @importFrom utils head read.csv
#' 
#' @return Combined normalised backbone and infinity markers for imputation
#' 
#' @details
#' Generating the combined data and saving to impu.input_log.mt.rds (on log scale) in the output directory.
#' 
mkImputeMT <-
function(
    paths
    ){
    bkc.adj.bkb <- readRDS(file = file.path(paths["downstream"], "bkc.adj.bkb_logScale_mt.rds")) #no. of selected bkb
    bkc.pe <- readRDS(file = file.path(paths["intermediary"], "bkc.pe_mt.rds"))
    metadata.cell <- readRDS(file = file.path(paths["intermediary"], "fcs_metadata_df.rds"))
    
    impu.input_log.mt <- cbind(bkc.adj.bkb, log(bkc.pe))
    rownames(impu.input_log.mt) <- metadata.cell$NO.in.all
    
    saveRDS(impu.input_log.mt, file = file.path(paths["intermediary"], "impu.input_log.mt.rds"))
    }
    