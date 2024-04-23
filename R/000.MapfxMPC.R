#' MAssively Parallel Flow cytometry Xplorer (MAPFX)
#'
#' @description This function is an end-to-end toolbox for analysing single-cell protein intensity data from the Massively-Parallel Cytometry (MPC) Experiments in FCS format. The functions include data normalisation, imputation (using backbone markers), and cluster analysis.
#'
#' @param runVignette logical; if FALSE (default), specify a path to `FCSpath` argument; TRUE for running Vignette using built-in data. 
#' @param runVignette_meta the argument for the built-in metadata when running Vignette; NULL (default).
#' @param runVignette_rawInten the argument for the built-in raw intensities when running Vignette; NULL (default).
#' @param FCSpath path to the input directory where `filename_meta.csv` and FCS files are stored (one file per well). `filename_meta.csv` should be saved under `FCSpath/FCS/meta/` and FCS files should be saved under `FCSpath/FCS/fcs/` (See Vignette for details.)
#' @param Outpath path to the output directory where intermediate results and final results will be stored.
#' @param file_meta if the file names of the FCS files are in the specified format, set `file_meta = "auto"`; otherwise set `file_meta = "usr"` and provide a `filename_meta.csv` file in `FCSpath/FCS/meta/`.
#' @param bkb.v a vector of the names of the backbone markers (MUST be the same as the names in the FCS files). For example, `bkb.v = c("FSC-H", "FSC-W", "SSC-H", "SSC-W", "CD69-CD301b", "MHCII", "CD4", "CD44", "CD8", "CD11c", "CD11b", "F480", "Ly6C", "Lineage", "CD45a488", "CD24", "CD103")`.
#' @param yvar the name of the well-specific exploratory marker in the FCS files (e.g., "Legend").
#' @param control.wells the well label of the control wells, including the autofluorescence and the isotype controls (format: plate_well, e.g., P1_A01). Users need to provide this information when `cluster.analysis.all = TRUE`. For example, `control.wells = c("P1_A01", "P2_A01", "P3_A01", "P3_F04", "P3_F05", "P3_F06", "P3_F07", "P3_F08", "P3_F09", "P3_F10", "P3_F11", "P3_F12", "P3_G01", "P3_G02")`.
#' @param bkb.upper.quantile the cut-off (default = 0.9) for selecting cells used for estimating the parameter of signal for backbone markers.
#' @param bkb.lower.quantile the cut-off (default = 0.1) for selecting cells used for estimating the parameters of noise for backbone markers.
#' @param bkb.min.quantile the cut-off (default = 0.01) for omitting the cells with the smallest values to minimise the impact of outliers during estimation (backbone).
#' @param inf.lower.quantile the cut-off (default = 0.1) for selecting cells used for estimating the parameters of noise for infinity markers.
#' @param inf.min.quantile the cut-off (default = 0.01) for omitting the cells with the smallest values to minimise the impact of outliers during estimation (infinity).
#' @param plots.bkc.bkb logical; if TRUE (default), produce scatter plots for pre- and post- background adjusted backbone markers (calibrated values on y-axis and raw values on x-axis).
#' @param plots.bkc.inf logical; if TRUE (default), produce scatter plots for pre- and post- background adjusted infinity markers (calibrated values on y-axis and raw values on x-axis).
#' @param plots.initM logical; if TRUE (default), produce an UMAP embedding plot to visualise the structure of the biological clusters used to form the initial M matrix for removal of well effects.
#' @param plots.rmWellEffect logical; if TRUE (default), produce heatmaps to visualise the unwanted (well) effects and biological effects in the pre- and post- adjusted datasets.
#' @param impute logical; if TRUE (default), impute the missing infinity markers.
#' @param models.use a vector of the names of the models used for imputation. For example, `models.use = c("LM", "LASSO3", "SVM", "XGBoost")`.
#' @param extra_args_regression_params a list of the lists of the parameters for running regression models. The order should be the same as the models specified in `models.use`. For example, `extra_args_regression_params = list(list(degree = 1), list(nfolds = 10, degree = 3), list(type = "nu-regression", cost = 8, nu = 0.5, kernel = "radial"), list(nrounds = 1500, eta = 0.03))`.
#' @param prediction_events_downsampling integer (default = NULL); the number of samples used for the downsampling for the prediction.
#' @param impu.training logical; if FALSE (default), not impute the training set (the dataset used to train the imputation models).
#' @param plots.imputation logical; if TRUE (default), visualise the distribution of R-sq values of infinity markers.
#' @param cluster.analysis.bkb  logical; if TRUE (default), perform cluster analysis using normalised backbone markers for all cells.
#' @param plots.cluster.analysis.bkb logical; if TRUE (default), produce an UMAP embedding plot from the normalised backbone markers to visualise the structure of the biological clusters for all cells.
#' @param cluster.analysis.all logical; must set `FALSE` if `impute = FALSE`; if TRUE (default), perform cluster analysis using normalised backbone markers and imputed infinity markers for cells in testing set.
#' @param plots.cluster.analysis.all logical; must set `FALSE` if `impute = FALSE`; if TRUE (default), produce an UMAP embedding plot from the normalised backbone markers and the imputed infinity markers to visualise the structure of the biological clusters for cells in testing set.
#' @param cores the number of cores used to perform parallel computation during the imputation process (default = 4L).
#'
#' @author Hsiao-Chi Liao, Agus Salim, and InfinityFlow (Becht et. al, 2021)
#'
#' @importFrom utils data
#'
#' @export
#' 
#' @return Normalised backbone markers on log scale, background noise corrected infinity markers, imputations, and metadata for cells
#' 
#' @details
#' In the output directory, this function produces the normalised backbone measurements, the background corrected infinity measurements, and imputed infinity markers (if set impute = TRUE), cell group labels from the cluster analysis using both normalised backbones and the completed dataset (if impute = TRUE), and graphs will be provided if specified.
#' 
#' @examples
#' # import built-in data
#' data(ord.fcs.raw.meta.df.out_mpc)
#' data(ord.fcs.raw.mt_mpc)
#' 
#' # create an Output directory for the MapfxMPC function
#' dir.create(file.path(tempdir(), "MPC_impu_Output"))
#' 
#' # When `impute = TRUE`, randomly selecting 50% of the cells in each well for model training
#' set.seed(123) 
#' MapfxMPC_impu_obj <- MapfxMPC(
#'   runVignette = TRUE, #set FALSE if not running this example
#'   runVignette_meta = ord.fcs.raw.meta.df.out_mpc, #set NULL if not running this example
#'   runVignette_rawInten = ord.fcs.raw.mt_mpc, #set NULL if not running this example
#'   FCSpath = NULL, # users specify their own input path if not running this example
#'   Outpath = file.path(tempdir(), "MPC_impu_Output"),
#'   file_meta = "auto",
#'   bkb.v = c(
#'     "FSC-H", "FSC-W", "SSC-H", "SSC-W", "CD69-CD301b", "MHCII", 
#'     "CD4", "CD44", "CD8", "CD11c", "CD11b", "F480", 
#'     "Ly6C", "Lineage", "CD45a488", "CD24", "CD103"),
#'   yvar = "Legend", 
#'   control.wells = c(
#'     "P1_A01", "P2_A01", "P3_A01",
#'     "P3_F04", "P3_F05", "P3_F06", "P3_F07", "P3_F08", 
#'     "P3_F09", "P3_F10", "P3_F11", "P3_F12",
#'     "P3_G01", "P3_G02"),
#'   bkb.upper.quantile = 0.9, 
#'   bkb.lower.quantile = 0.1, 
#'   bkb.min.quantile = 0.01,
#'   inf.lower.quantile = 0.1, 
#'   inf.min.quantile = 0.01, 
#'   plots.bkc.bkb = TRUE, plots.bkc.inf = TRUE, 
#'   plots.initM = TRUE,
#'   plots.rmWellEffect = TRUE,
#'   impute = TRUE,
#'   models.use = c("XGBoost"),
#'   extra_args_regression_params = list(list(nrounds = 1500, eta = 0.03)),
#'   prediction_events_downsampling = NULL,
#'   impu.training = FALSE,
#'   plots.imputation = TRUE,
#'   cluster.analysis.bkb = TRUE, plots.cluster.analysis.bkb = TRUE,
#'   cluster.analysis.all = TRUE, plots.cluster.analysis.all = TRUE,
#'   cores = 2L)
#' 
#' 
MapfxMPC <-
function(
    runVignette = FALSE,
    runVignette_meta = NULL,
    runVignette_rawInten = NULL,
    FCSpath = NULL,
    Outpath = NULL,
    file_meta = "auto",
    bkb.v = NULL,
    yvar = "Legend",
    control.wells = NULL,
    bkb.upper.quantile = 0.9, 
    bkb.lower.quantile = 0.1, 
    bkb.min.quantile = 0.01,
    inf.lower.quantile = 0.1, 
    inf.min.quantile = 0.01, 
    plots.bkc.bkb = TRUE, plots.bkc.inf = TRUE, 
    plots.initM = TRUE,
    plots.rmWellEffect = TRUE,
    impute = TRUE,
    models.use = c("XGBoost"),
    extra_args_regression_params = list(list(nrounds = 1500, eta = 0.03)),
    prediction_events_downsampling = NULL,
    impu.training = FALSE,
    plots.imputation = TRUE,
    cluster.analysis.bkb = TRUE, plots.cluster.analysis.bkb = TRUE,
    cluster.analysis.all = TRUE, plots.cluster.analysis.all = TRUE,
    cores = 4L){
    
    chans <- bkb.v
    
    #run Vignette
    if(runVignette == TRUE){
    if(is.null(Outpath)){
    message("\nPlease make sure the output path has been specified!")
    }
    
    message("\n\n\nCreating directories for output...")
    settings <- initialize(
    path_to_fcs=FCSpath,
    path_to_output=Outpath,
    verbose=TRUE)
    paths <- settings$paths
    
    ###
    
    saveRDS(runVignette_meta, file = file.path(paths["intermediary"], "fcs_metadata_df.rds"))
    saveRDS(runVignette_rawInten, file = file.path(paths["intermediary"], "fcs_rawInten_mt.rds"))
    }
    
    #not running Vignette
    if(runVignette == FALSE){
    
    if(is.null(FCSpath) | is.null(Outpath) | is.null(bkb.v)){
    stop("\nPlease make sure you have specified\n(1) path to the input files (details in the folder diagram in the vignette) \n(2) output path \n(3) a list of backbone markers")
    }
    
    if(is.null(file_meta)){
    stop("\nThe file_meta argument only accepts either \"auto\" or \"usr\".\nAlso, if file_meta=\"usr\", make sure there is a filename_meta.csv file in FCSpath/meta (details in the folder diagram in the vignette).")
    }else if(file_meta != "auto" & file_meta != "usr"){
    stop("\nThe file_meta argument only accepts either \"auto\" or \"usr\".\nAlso, if file_meta=\"usr\", make sure there is a filename_meta.csv file in FCSpath/meta (details in the folder diagram in the vignette).")
    }
    
    if(!is.null(FCSpath) & !is.null(Outpath) & !is.null(bkb.v)){
    
    if(cluster.analysis.all == TRUE & is.null(control.wells)){
    stop("Please specify your control wells because they will be removed for cluster analysis...")
    }
    
    if(length(extra_args_regression_params) != length(models.use)){
    stop("'extra_args_regression_params' and 'models.use' should be lists of the same lengths")
    }
    
    ##/!\ Potentially add a check here to make sure parameters are consistent with FCS files
    message("\n\n\nCreating directories for output...")
    
    settings <- initialize(
    path_to_fcs=FCSpath,
    path_to_output=Outpath,
    verbose=TRUE)
    paths <- settings$paths
    
    ###
    
    ##1.fcs to rds
    message("\n\n\nTransforming FCS to RDS files...")
    fcs_to_rds(paths=paths, file_meta=file_meta, yvar = yvar)
    }
    }
    
    ##2-1.bkc bkb
    message("\n\n\nBackground correcting backbone markers...")
    bkc_bkb(
    paths=paths, bkb.v=bkb.v,
    MPC=TRUE,
    bkb.upper.quantile=bkb.upper.quantile, #cells used for estimating parameter of signal
    bkb.lower.quantile=bkb.lower.quantile, #cells used for estimating parameters of noise
    bkb.min.quantile=bkb.min.quantile,
    plots=plots.bkc.bkb) #the lowest 1% of values will not be used to minimise the impact of outliers on sig
    
    ##2-2.bkc inf
    message("\n\n\nBackground correcting infinity markers...")
    bkc_inf <- bkc_pe(
    paths=paths,
    pe.lower.quantile=inf.lower.quantile, #cells used for estimating parameters of noise
    pe.min.quantile=inf.min.quantile, #the lowest 1% of values will not be used to minimize the impact of outliers on sig
    plots=plots.bkc.inf)
    
    ##3-0.initM
    message("\n\n\nForming a matrix of biology (M) for removal of well effect...")
    metadata <- initM(
    paths=paths, assay="MPC", bkb.v=bkb.v,
    plots=plots.initM) #may give "cent.log.bkc" in the next version
    
    ##3-1.rmWellEffect
    message("\n\n\nRemoval of well effect...")
    norm.bkb <- rmWellEffect(
    paths=paths,
    plots = plots.rmWellEffect)
    
    if(impute == TRUE){
    message("\n\n\nImputation got started...")
    ##4.prep. for imputation - combine normalised.bkb.inf
    mkImputeMT(paths=paths)
    
    ##5.impu. using bkb as predictors
    impu <- imputation_bkb.predictors(
    paths=paths,
    chans=chans,
    yvar="Legend", #have named this in my fcs_to_rds function
    cores=cores,
    models.use = models.use,
    extra_args_regression_params = extra_args_regression_params,
    prediction_events_downsampling = prediction_events_downsampling,
    impu.training = impu.training,
    plots = plots.imputation)
    
    imputation <- impu
    
    ##
    ##6.cluster analysis for completed dataset (bkb, bkb+impu.inf) using cells from testing set
    message("\n\n\nCluster analysis with adjusted backbone markers and completed dataset for cells in the testing set...")
    if(cluster.analysis.all == TRUE){
    metadata <- cluster.analysis(
    paths=paths,
    bkb.v=bkb.v,
    yvar="Legend",
    control.wells=control.wells,
    plots = plots.cluster.analysis.all)
    message("\nCell group labels are saved in \"GP.denoised.bkb\" and \"GP.denoised.bkb.impuInf*\" columns...")
    }
    
    }else{
    imputation <- NA
    }

    ##7.cluster analysis bkbOnly using all cells
    message("\n\n\nCluster analysis with adjusted backbone markers for ALL cells...")
    if(cluster.analysis.bkb == TRUE){
    metadata <- cluster.analysis.bkbOnly(
    paths=paths,
    bkb.v=bkb.v,
    plots = plots.cluster.analysis.bkb)
    }
    message("\nCell group labels are saved in \"GP.denoised.bkb.allCells\" column...")
    
    message("\tCompleted!")
    return(list(logNormalised_bkb=norm.bkb, BKC_inf=bkc_inf, Imputation=imputation, Metadata=metadata))
    }

