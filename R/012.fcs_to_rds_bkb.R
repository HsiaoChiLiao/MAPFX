#' Converting FCS files to RDS files (for the case without exploratory markers)
#'
#' This function has been designed to convert the raw FCS files to data matrix and export to RDS files.
#'
#' @param paths a vector of characters of paths to store intput, intermediary results, outputs...
#' @param file_meta if the file names of the FCS files are in the specified format, set file_meta="auto"; otherwise set file_meta="usr" and provide "filename_meta.csv" in FCSpath.
#' @param MPC if the data is from MPC experiments, set MPC = TRUE. Setting FALSE represents data from the fluorescence flow cytometry (FFC) assay.
#'
#' @author Hsiao-Chi Liao
#'
#' @importFrom flowCore read.flowSet keyword parameters
#' @importFrom Biobase pData exprs
#' @importFrom stringr str_extract
#' 
#' @importFrom grDevices dev.off jpeg
#' @importFrom graphics abline
#' @importFrom methods is
#' @importFrom stats as.formula contr.sum contrasts<- dexp dnorm median model.matrix optim pexp pnorm quantile sd setNames
#' @importFrom utils head read.csv
#'
#' @return Raw protein intensities and the corresponding metadata from FFC experiments
#'
#' @details
#' Generating fcs_metadata_df.rds and fcs_rawInten_mt.rds files in the output directory.
#'
fcs_to_rds_bkb <-
function(
    paths, file_meta, MPC
    ){
    ## re-make large mt and consider SSC and FSC
    # import fcs.impu.raw files ###
    fs <- read.flowSet(path = file.path(paths["input"], "fcs"))
    name.desc.fs <- pData(parameters(fs[[1]]))
    name.desc.fs$desc[which(is.na(name.desc.fs$desc))] <- name.desc.fs$name[which(is.na(name.desc.fs$desc))]
    
    ####

    #### data - raw ####
    fcs.raw.mt.list <- fcs.raw.meta.list <- list()
    
    if((MPC == FALSE) & (file_meta == "usr")){
    meta.inf <- read.csv(file = file.path(paths["input"], "meta", "filename_meta.csv"))
    for(i in seq_along(fs)){
    fcs.raw.mt.list[[i]] <- exprs(fs[[i]])
    filenam <- unlist(keyword(fs[[i]], "GUID"))
    fcs.raw.meta.list[[i]] <- cbind(
    NO.in.batch=(seq_len(dim(fcs.raw.mt.list[[i]])[1])),
    Filenam=filenam,
    Batch=meta.inf[which(filenam == meta.inf$Filenam), "Batch"])
    }
    fcs.raw.mt <- do.call(rbind, fcs.raw.mt.list)
    colnames(fcs.raw.mt) <- name.desc.fs$desc[match(colnames(fcs.raw.mt), name.desc.fs$name)]
    ####
    
    # making metadata ###
    fcs.raw.meta.df <- as.data.frame(do.call(rbind, fcs.raw.meta.list))
    fcs.raw.meta.df$NO.in.batch <- as.integer(fcs.raw.meta.df$NO.in.batch)
    #use lapply for multiple changing class: lapply(mydf[,2:3], as.factor)
    
    #ordering before giving NO.in.all
    ord.fcs.raw.meta.df <- fcs.raw.meta.df[order(fcs.raw.meta.df$Batch),,drop=FALSE]
    ord.fcs.raw.mt <- fcs.raw.mt[order(fcs.raw.meta.df$Batch),,drop=FALSE]
    
    #global ID for each cell
    ord.fcs.raw.meta.df$NO.in.all <- seq_len(nrow(ord.fcs.raw.mt))
    
    #out rds
    ord.fcs.raw.meta.df.out <- ord.fcs.raw.meta.df[,c(4,1,3),drop=FALSE]
    }else{
    stop("Please provide a \"filename_meta.csv\" file in `FCSpath/meta/`")
    }
    
    saveRDS(ord.fcs.raw.meta.df.out, file = file.path(paths["intermediary"], "fcs_metadata_df.rds"))
    saveRDS(ord.fcs.raw.mt, file = file.path(paths["intermediary"], "fcs_rawInten_mt.rds"))
    }
