#' Normalising data from the Fluorescence Flow Cytometry (FFC) Experiments with mapfx.norm
#'
#' @description This function is used to normalise, including background correction and removal of batch effects, protein intensity data from FFC assays. The input data is in FCS format. The functions include data normalisation and cluster analysis.
#'
#' @param FCSpath path to the input directory where `filename_meta.csv` and FCS files are stored. `filename_meta.csv` should be saved under `FCSpath/FCS/meta/` and FCS files should be saved under `FCSpath/FCS/fcs/`. See Vignette for details.
#' @param Outpath path to the output directory where intermediate results and final results will be stored.
#' @param file_meta if the file names of the FCS files are in the specified format, set `file_meta = "auto"`; otherwise set `file_meta = "usr"` and provide a `filename_meta.csv` file in `FCSpath/FCS/meta/`.
#' @param protein.v a vector of the names of the protein markers (MUST be the same as the names in the FCS files). For example, `protein.v = c("FSC-H","FSC-W","SSC-H","SSC-W","CD3","CD4","CD8","CD45")`.
#' @param protein.upper.quantile the cut-off (default = 0.9) for selecting cells used for estimating the parameter of signal for protein markers.
#' @param protein.lower.quantile the cut-off (default = 0.1) for selecting cells used for estimating the parameters of noise for protein markers.
#' @param protein.min.quantile the cut-off (default = 0.01) for omitting the cells with the smallest values to minimise the impact of outliers during estimation.
#' @param plots.bkc.protein logical; if TRUE (default), produce scatter plots for pre- and post- background adjusted protein markers (calibrated values on y-axis and raw values on x-axis).
#' @param plots.initM logical; if TRUE (default), produce an UMAP embedding plot to visualise the structure of the biological clusters used to form the initial M matrix for removal of batch effects.
#' @param plots.rmBatchEffect logical; if TRUE (default), produce heatmaps to visualise the unwanted (batch) effects and biological effects in the pre- and post- adjusted datasets.
#' @param cluster.analysis.protein  logical; if TRUE (default), perform cluster analysis using normalised protein markers.
#' @param plots.cluster.analysis.protein logical; if TRUE (default), produce an UMAP embedding plot from the normalised protein markers to visualise the structure of the biological clusters.
#'
#' @author Hsiao-Chi Liao, Agus Salim
#'
#' @export
#' @return Normalised protein and infinity measurements for the MPC data; normalised protein measurements for the FFC data. Cluster analysis for normalised proteins. Graphs will be provided if specified.
#' @examples
#' MapfxFFC()
#' 
#' @usage
#' MapfxFFC(FCSpath,
#'          Outpath,
#'          file_meta="auto",
#'          protein.v,
#'          protein.upper.quantile = 0.9, 
#'          protein.lower.quantile = 0.1, 
#'          protein.min.quantile = 0.01,
#'          plots.bkc.protein = TRUE,
#'          plots.initM = TRUE,
#'          plots.rmBatchEffect = TRUE,
#'          cluster.analysis.protein = TRUE, plots.cluster.analysis.protein = TRUE)
#'
MapfxFFC <-
function(
    FCSpath = NULL,
    Outpath = NULL,
    file_meta = "auto",
    protein.v = NULL,
    protein.upper.quantile = 0.9, 
    protein.lower.quantile = 0.1, 
    protein.min.quantile = 0.01,
    plots.bkc.protein = TRUE,
    plots.initM = TRUE,
    plots.rmBatchEffect = TRUE,
    cluster.analysis.protein = TRUE, plots.cluster.analysis.protein = TRUE){
    
    if(is.null(FCSpath) | is.null(Outpath) | is.null(protein.v)){
    message("\nPlease make sure you have specified\n(1) path to the FCS files, and metadata (if file_meta=\"usr\") \n(2) output path \n(3) a list of protein markers")
    }
    
    if(!is.null(FCSpath) & !is.null(Outpath) & !is.null(protein.v)){
    
    settings <- initialize(
    path_to_fcs=FCSpath,
    path_to_output=Outpath,
    verbose=TRUE)
    
    paths <- settings$paths
    
    ###
    
    ##1.fcs to rds
    fcs_to_rds_bkb(
    paths=paths, file_meta=file_meta, MPC = FALSE)
    
    
    ##2-1.bkc protein
    bkc_bkb_ffc(
    paths=paths, bkb.v=protein.v,
    bkb.upper.quantile=protein.upper.quantile, #cells used for estimating parameter of signal
    bkb.lower.quantile=protein.lower.quantile, #cells used for estimating parameters of noise
    bkb.min.quantile=protein.min.quantile,
    plots=plots.bkc.protein) #the lowest 1% of values will not be used to minimise the impact of outliers on sig
    
    ##3-0.initM
    initM(
    paths=paths, assay="FFC", bkb.v=protein.v,
    plots=plots.initM) #may give "cent.log.bkc" in the next version
    
    ##3-1.rmBatchEffect
    rmBatchEffect(
    paths=paths,
    plots = plots.rmBatchEffect)
    
    ##7.cluster analysis proteinOnly
    if(cluster.analysis.protein == TRUE){
    cluster.analysis.bkbOnly(
    paths=paths,
    bkb.v=protein.v,
    plots = plots.cluster.analysis.protein)
    }
    
    message("\tCompleted!")  
    }
}