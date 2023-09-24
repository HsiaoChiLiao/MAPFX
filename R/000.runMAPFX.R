#' MAssively Parallel Flow cytometry Xplorer (MAPFX) for analysing single-cell protein intensity from the Massively-Parallel Cytometry (MPC) Experiments - imputation using backbone markers
#'
#' @description This function is an end-to-end toolbox for analysing the data from the LEGENDScreen assay in FCS format, including data normalisation, imputation and cluster analysis.
#'
#' @param FCSpath Path to the input directory where input FCS files are stored (one file per well). Will look for FCS files recursively in that directory.
#' @param Outpath Path to the output directory where intermediate results and final results will be stored
#' @param file_meta If the file names of the FCS files are in the specified format, set file_meta="auto"; otherwise set file_meta="usr" and provide "filename_meta.csv" in FCSpath.
#' @param bkb.v A vector of the names of the backbone markers (MUST match to the names in the FCS file).
#' @param chans A vector of the names of the backbone markers (MUST match to the names in the FCS file).
#' @param yvar The name of the well-specific marker in the FCS files (e.g., "Legend").
#' @param control.wells The well label of the control wells, including the autofluorescence and the isotype controls (format: plate_well, e.g., P1_A01)
#' @param bkb.upper.quntile The cut-off (default = 0.9) for selecting cells used for estimating the parameter of signal.
#' @param bkb.lower.quntile The cut-off (default = 0.1) for selecting cells used for estimating the parameters of noise.
#' @param bkb.min.quntile The cut-off (default = 0.01) for omitting the cells with the smallest values to minimise the impact of outliers.
#' @param pe.mean.sd Selecting cells with the value larger than mean + ?sd (default = 3) to estimate the parameter of signal.
#' @param pe.lower.quntile The cut-off (default = 0.1) for selecting cells used for estimating the parameters of noise.
#' @param pe.min.quntile The cut-off (default = 0.01) for omitting the cells with the smallest values to minimise the impact of outliers.
#' @param trans.dat Method for transforming raw intensities (default: Logicle).
#' @param models.use A vector of the names of the models used for imputation (an example: c("LM", "LASSO3", "SVM", "XGBoost")). The length of the vector is arbitrary.
#' @param extra_args_regression_params A list of the lists of the parameters for running regression models. (default = list(list(degree = 1), list(nfolds = 10, degree = 3),
#' list(type = "nu-regression", cost = 8, nu = 0.5, kernel = "radial"),
#' list(nrounds = 1500, eta = 0.03)))
#' @param cores The number of cores used to perform parallel computation (default = 4L).
#'
#' @author Hsiao-Chi Liao, Agus Salim, and InfinityFlow (Becht et. al, 2021)
#'
#' @export
#' @return Normalised backbone measurements and imputed well-specific markers. Cluster analysis for both normalised backbones and the completed dataset. Graphs will be provided.
#'
#' @usage
#' runMAPFX(FCSpath="/PathToFCSfiles/",
#' Outpath="/PathToOutputFolder/Output/",
#' file_meta="auto",
#' bkb.v = c("FSC-H", "FSC-W", "SSC-H", "SSC-W", "CD69-CD301b", "MHCII", "CD4", "CD44",
#'           "CD8", "CD11c", "CD11b", "F480", "Ly6C", "Lineage", "CD45a488",
#'           "CD24", "CD103"),
#' chans = c("FSC-H", "FSC-W", "SSC-H", "SSC-W", "CD69-CD301b", "MHCII", "CD4", "CD44",
#'           "CD8", "CD11c", "CD11b", "F480", "Ly6C", "Lineage", "CD45a488",
#'           "CD24", "CD103"),
#' yvar="Legend", control.wells = c("P1_A01", "P2_A02", "P3_G02"),
#' bkb.upper.quntile=0.9, bkb.lower.quntile=0.1, bkb.min.quntile=0.01,
#' pe.mean.sd=3, pe.lower.quntile=0.1, pe.min.quntile=0.01,
#' trans.dat="cent.lgc",
#' models.use = c("XGBoost"),
#' extra_args_regression_params = list(list(nrounds = 1500, eta = 0.03)),
#' cores=2L)
#'
#'
runMAPFX <-
function(FCSpath,
                         Outpath,
                         file_meta="auto",
                         bkb.v,
                         chans,
                         yvar="Legend",
                         control.wells,
                         bkb.upper.quntile=0.9, bkb.lower.quntile=0.1, bkb.min.quntile=0.01,
                         pe.mean.sd=3, pe.lower.quntile=0.1, pe.min.quntile=0.01,
                         trans.dat="cent.lgc",
                         models.use = c("LM", "LASSO3", "SVM", "XGBoost"),
                         extra_args_regression_params = list(list(degree = 1), list(nfolds = 10, degree = 3),
                                                             list(type = "nu-regression", cost = 8, nu = 0.5, kernel = "radial"),
                                                             list(nrounds = 1500, eta = 0.03)),
                         cores=4L){
      
      ## Making sure that optional dependencies are installed if used. (not sure if needed yet)
      # lapply(
      #   regression_functions,
      #   function(fun){
      #     fun(x = NULL, params = NULL)
      #   }
      # )
      
      # if(length(extra_args_regression_params) != length(regression_functions)){
      #   stop("extra_args_regression_params and regression_functions should be lists of the same lengths")
      # }
      
      if(length(extra_args_regression_params) != length(models.use)){
        stop("'extra_args_regression_params' and 'models.use' should be lists of the same lengths")
      }
      
      ##/!\ Potentially add a check here to make sure parameters are consistent with FCS files
      
      settings <- initialize(
        path_to_fcs=FCSpath,
        path_to_output=Outpath,
        verbose=TRUE
      )
      paths <- settings$paths
      
      ###
      
      ##1.fcs to rds
      fcs_to_rds(paths, file_meta=file_meta, yvar=yvar)
      
      ##2-1.bkc bkb
      bkc_bkb(paths, bkb.v,
              bkb.upper.quntile=bkb.upper.quntile, #cells used for estimating parameter of signal
              bkb.lower.quntile=bkb.lower.quntile, #cells used for estimating parameters of noise
              bkb.min.quntile=bkb.min.quntile) #the lowest 1% of values will not be used to minimise the impact of outliers on sig
      
      ##2-2.bkc pe
      bkc_pe(paths,
             pe.mean.sd=pe.mean.sd, #mean+-3sd for extracting extremely large values
             # upper.quntile=0.9, #cells used for estimating parameter of signal
             pe.lower.quntile=pe.lower.quntile, #cells used for estimating parameters of noise
             pe.min.quntile=pe.min.quntile) #the lowest 1% of values will not be used to minimize the impact of outliers on sig
      
      ##3-0.initM
      initM(paths, bkb.v,
            trans.dat=trans.dat) #may give "cent.log.bkc" in the future
      
      ##3-1.rmWellEffect
      rmWellEffect(paths,
                   visualisation = TRUE)
      
      ##4.mkImputMT
      mkImputeMT(paths)
      
      ##5.imputation_bkb.pred
      imputation_bkb.predictors(paths,
                                chans,
                                yvar=yvar,
                                cores=cores,
                                models.use = models.use,
                                extra_args_regression_params = extra_args_regression_params,
                                prediction_events_downsampling=NULL,
                                impu.training = F)
      
      ##6.cluster analysis
      cluster.analysis(paths,
                       bkb.v,
                       yvar=yvar,
                       control.wells)
      
      
      # timings <- cbind(fit=timings_fit$timings,pred=timings_pred$timings)
      # rownames(timings) <- names(regression_functions)
      # saveRDS(timings,file.path(paths["rds"],"timings.Rds"))
      
      message("\tCompleted!")
    }
