#' Imputing the unmeasured well-specific markers with regression models
#'
#' This function has been designed to impute/predict the unmeasured well-specific markers with regression models.
#'
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param chans A vector of the names of the backbone markers (MUST match to the names in the FCS file).
#' @param yvar The name of the well-specific marker in the FCS files (has been changed to "Legend" in the first function).
#' @param cores The number of cores used to perform parallel computation (default = 8L).
#' @param models.use A vector of the names of the models used for imputation (an example: c("LM", "LASSO3", "SVM", "XGBoost")). The length of the vector is arbitrary.
#' @param extra_args_regression_params A list of the lists of the parameters for running regression models.
#' @param prediction_events_downsampling Default = NULL (not doing subsampling). How many cells per well you want to have the imputation? (must be less than or equal to a half as we won't get the prediction for cells in the training set).
#'
#' @author Hsiao-Chi Liao and InfinityFlow (Becht et. al, 2021)
#'
#' @import glmnetUtils
#' @import e1071
#' @import xgboost
#' @import foreach
#' @import doParallel
#' @importFrom stats lm predict
#' @importFrom utils getS3method
#' @importFrom parallel makeCluster clusterExport clusterEvalQ stopCluster
#' @importFrom pbapply pblapply
#' @import ggplot2
#' @importFrom reshape melt
#' @importFrom gtools combinations
#'
#' @return Imputing the unmeasured well-specific markers and save to predictions.Rds file. Visualising the result with the boxplots (r-sq).
#'
imputation_bkb.predictors <-
function(
    paths,
    chans,
    yvar="Legend",
    cores=cores,
    models.use = models.use,
    extra_args_regression_params = extra_args_regression_params,
    prediction_events_downsampling=NULL,
    impu.training = F
    ){
      # if(verbose){
      #   message("\tFitting regression models")
      # }
      
      # chans <- c("FSC-H", "FSC-W", "SSC-H", "SSC-W",
      #           "CD69-CD301b", "MHCII", "CD4", "CD44", "CD8",
      #           "CD11c", "CD11b", "F480", "Ly6C", "Lineage", "CD45a488",
      #           "CD24", "CD103")
      # yvar <- "Legend"
      # models.use <- "XGBoost" #c("LM", "LASSO3", "SVM", "XGBoost")
      # extra_args_regression_params <- list(list(nrounds = 1500, eta = 0.03))
      # models.use <- "LM" #c("LM", "LASSO3", "SVM", "XGBoost")
      # extra_args_regression_params <- list(list(degree = 1))
      # models.use <- c("LM", "LASSO3", "SVM", "XGBoost")
      # extra_args_regression_params <- list(list(degree = 1), list(nfolds = 10, degree = 3),
      #                                      list(type = "nu-regression", cost = 8, nu = 0.5, kernel = "radial"),
      #                                      list(nrounds = 1500, eta = 0.03))
      
      # cores=8L
      # prediction_events_downsampling=NULL
      
      # RDSpath="/Users/chelsea/Desktop/PhD/ProteinExpression/InFlow/Programming/infinity/R_CODE/72.mk_pkg/data/rds/"
      # Graphpath="/Users/chelsea/Desktop/PhD/ProteinExpression/InFlow/Programming/infinity/R_CODE/72.mk_pkg/data/graph/Imputation/"
      xp <- readRDS(file = file.path(paths["intermediary"], "impu.input_log.mt.rds"))
      metadata.cell <- readRDS(file = file.path(paths["intermediary"], "fcs_metadata_df.rds"))
      params <- extra_args_regression_params
      
      ### borrowing functions from inflow pkg ###
      ## fitters (may add arbitrary parameters in the future)
      {
        polynomial_formula <- function(variables,degree){
          n <- length(variables)
          polys <- lapply(
            seq_len(degree),
            function(deg){
              res <- apply(gtools::combinations(n,deg,variables,repeats.allowed=TRUE),1,paste,collapse="*")
              res <- paste0("I(",res,")")
            }
          )
          paste(do.call(c,lapply(polys,paste,sep="+")),collapse="+")
        }
        
        ## linear fitter
        fitter_linear <- function(x = NULL, params = NULL){
          if(!is.null(x) & !is.null(params)){
            w <- x[,"train_set"]==1
            fmla <- paste0(make.names(yvar),"~",polynomial_formula(variables=chans,degree=params$degree))
            model <- lm(formula=fmla,data=as.data.frame(x[w,c(chans,yvar)]))
            pred <- predict(model,as.data.frame(x[,chans]))
            rm(list=setdiff(ls(),c("pred","model")))
            model$model=NULL ## Trim down for slimmer objects
            model$qr$qr=NULL ## Trim down for slimmer objects
            return(list(pred=pred,model=model))
          }
        }
        ## lasso3 fitter
        fitter_glmnet <- function(x = NULL, params = NULL){
          if(!requireNamespace("glmnetUtils", quietly = TRUE)){
            stop("Please run install.packages(\"glmnetUtils\")")
          }
          if(!is.null(x) & !is.null(params)){
            w <- x[,"train_set"] == 1
            fmla <- paste0(make.names(yvar), "~", polynomial_formula(variables = chans, degree = params$degree))
            flma <- as.formula(fmla)
            params <- params[setdiff(names(params), "degree")]
            params <- c(
              params,
              list(
                formula = fmla,
                data = as.data.frame(x[w, c(chans, yvar)]),
                use.model.frame = TRUE
              )
            )
            
            fun <- getS3method("cv.glmnet", "formula", envir = asNamespace("glmnetUtils"))
            model <- do.call(fun, params)
            model$call <- NULL ## Slimming down object
            model$glmnet.fit$call <- NULL ## Slimming down object
            attributes(model$terms)[[".Environment"]] <- NULL ## Slimming down object
            pred <- predict(model, as.data.frame(x[, chans]), s = model$lambda.min)
            
            rm(list = setdiff(ls(), c("pred", "model")))
            return(list(pred = pred, model = model))
          }
        }
        ## svm fitter
        fitter_svm <- function(x = NULL, params = NULL){
          if(!requireNamespace("e1071", quietly = TRUE)){
            stop("Please run install.packages('e1071')")
          }
          if(!is.null(x) & !is.null(params)){
            w <- x[,"train_set"]==1
            model <- do.call(function(...){e1071::svm(...,x=x[w,chans],y=x[w,yvar])},params)
            pred <- predict(model,x[,chans])
            rm(list=setdiff(ls(),c("pred","model")))
            return(list(pred=pred,model=model))
          }
        }
        ## xgboost fitter
        fitter_xgboost <- function(x = NULL, params = NULL){
          if(!requireNamespace("xgboost", quietly = TRUE)){
            stop("Please run install.packages(\"xgboost\")")
          }
          
          if(!is.null(x) & !is.null(params)){
            w <- x[,"train_set"]==1
            args <- c(list(data = x[w, chans], label = x[w, yvar], nthread = 1L, verbose = 0), params)
            model <- do.call(xgboost::xgboost, args)
            pred <- predict(model, x[, chans])
            rm(list=setdiff(ls(),c("pred","model")))
            return(list(pred=pred,model=model))
          }
        }
        
        # regression_functions <- list(LM = fitter_linear, LASSO3 = fitter_glmnet, SVM = fitter_svm, XGBoost = fitter_xgboost)
        # extra_args_regression_params = list(list(degree = 1), list(nfolds = 10, degree = 3),
        #                                     list(type = "nu-regression", cost = 8, nu = 0.5, kernel = "radial"), list(nrounds = 1500, eta = 0.03))
      }
      ## predict_wrapper
      {
        ## xp is the matrix of predictors, x is the built model
        predict_wrapper <- function(x){
          if(is(x, "lm")){
            xp <- as.data.frame(xp)
          }
          if(is(x, "cv.glmnet")){
            requireNamespace("glmnetUtils")
            xp <- as.data.frame(xp)
            return(predict(x,xp,s=x$lambda.min)[,1])
          }
          # if(is(x, "raw")){ #for NN
          #   requireNamespace("keras")
          #   x = keras::unserialize_model(x)
          # }
          if(is(x, "xgb.Booster")){
            requireNamespace("xgboost")
            x <- xgboost::xgb.Booster.complete(x)
            xgboost::xgb.parameters(x) <- list(nthread = 1)
          }
          if(is(x, "svm")){
            requireNamespace("e1071")
          }
          res <- predict(x, xp)
        }
      }
      ## split function for selecting training data
      {
        split_matrix <- function (mat, vector, byrow = TRUE)
        {
          if (byrow & nrow(mat) != length(vector)) {
            stop("if byrow=TRUE, vector's length should have length nrow(mat)")
          }
          else if (!byrow & ncol(mat) != length(vector)) {
            !byrow & ncol(mat) != length(vector)
            stop("if byrow=FALSE, vector's length should have length ncol(mat)")
          }
          if (byrow) {
            levels  <-  split(seq_len(nrow(mat)), vector) #rows in mat should match order in vector
            res  <-  lapply(levels, function(x) {
              mat[x, , drop = FALSE]
            })
          }
          else {
            levels  <-  split(seq_len(ncol(mat)), vector)
            res  <-  lapply(levels, function(x) {
              mat[, x, drop = FALSE]
            })
          }
          res
        }
      }
      
      ## constructing the list of regression_functions
      regression_functions <- list()
      for(m in 1:length(models.use)){
        if(models.use[m] == "LM"){
          regression_functions[[m]] <- fitter_linear
        }else if(models.use[m] == "LASSO3"){
          regression_functions[[m]] <- fitter_glmnet
        }else if(models.use[m] == "SVM"){
          regression_functions[[m]] <- fitter_svm
        }else{
          regression_functions[[m]] <- fitter_xgboost
        }
      }
      names(regression_functions) <- models.use
      
      
      message("\tFitting regression models...")
      chans <- make.names(chans)
      yvar <- make.names(yvar)
      colnames(xp) <- make.names(colnames(xp))
      events.code <- metadata.cell$Well.lab
      
      #?
      requireNamespace("parallel")
      
      cl <- parallel::makeCluster(min(cores,length(unique(events.code))))
      
      ## random number generation
      RNGkind("L'Ecuyer-CMRG")
      
      ## OS related
      {
        if(.Platform$OS.type == "windows") {
          mc.reset.stream <- function() return(invisible(NULL))
        } else {
          mc.reset.stream <- parallel::mc.reset.stream
        }
        
        mc.reset.stream()
        
        env <- environment()
      }
      
      ## random selecting a half of the input (stratified by wells) for training
      message("\tRandomly selecting 50% of the cells in each well for model training...")
      {
        d.e <- split_matrix(xp, events.code) #list(matrix for cells from well 1, ... well 2, ...)
        d.e <- lapply( # add another column: train_set~ 0: testing, 1: training
          d.e,
          function(x){
            set.seed(123)
            w <- sample(rep(c(TRUE,FALSE),times=c(floor(nrow(x)/2),nrow(x)-floor(nrow(x)/2))))
            x <- cbind(x,train_set=ifelse(w,1,0))
            x
          }
        )
        
        train_set <- matrix(
          ncol=1,
          dimnames=list(NULL,"train_set"),
          do.call(
            c,
            lapply(
              d.e,
              function(x)
              {
                x[,"train_set"]
              }
            )
          )
        )
        # train_set
        # [1,]         0
        # [2,]         1
        # [3,]         1
      }
      
      ## not sure what these used for... about parallel computing??
      {
        clusterExport(
          cl,
          c("yvar","chans", "regression_functions", "polynomial_formula"),
          # c("yvar","chans","neural_networks_seed", "regression_functions", "fitter_nn"),
          envir=env
        )
        
        clusterEvalQ(
          cl,
          {
            chans <- make.names(chans)
            yvar <- make.names(yvar)
            ## if(any(sapply(regression_functions, function(x){identical(x, fitter_nn)}))){
            ##     if(requireNamespace("keras", quietly = TRUE) & requireNamespace("tensorflow", quietly = TRUE)){
            ##         if(!is.null(neural_networks_seed)){
            ##             tensorflow::use_session_with_seed(neural_networks_seed) ## This will make results reproducible, disable GPU and CPU parallelism (which is good actually). Source: https://keras.rstudio.com/articles/faq.html#how-can-i-obtain-reproducible-results-using-keras-during-development
            ##         } else {
            ##             tensorflow::tf$reset_default_graph()
            ##             config <- list()
            ##             config$intra_op_parallelism_threads <- 1L
            ##             config$inter_op_parallelism_threads <- 1L
            ##             session_conf <- do.call(tensorflow::tf$ConfigProto, config)
            ##             sess <- tensorflow::tf$Session(graph = tensorflow::tf$get_default_graph(), config = session_conf)
            ##             tensorflow:::call_hook("tensorflow.on_use_session", sess, TRUE)
            ##         }
            ##     }
            ## }
          }
        )
      }
      
      ## fitting models
      message("\t\tFitting...")
      {
        models <- list()
        timings <- numeric()
        
        for(i in seq_along(regression_functions)){
          cat("\t\t",names(regression_functions)[i],"\n\n",sep="")
          t0 <- Sys.time()
          fun <- regression_functions[[i]]
          environment(fun) <- environment() ## Fixing issue with scoping when cores = 1L
          models[[i]] <- pblapply(
            X=d.e,
            FUN=fun,
            params=params[[i]],
            cl=cl
          )
          t1 <- Sys.time()
          dt <- difftime(t1,t0,units="secs")
          cat("\t",dt," seconds","\n",sep="")
          timings <- c(timings,dt)
        }
        
        names(models) <- names(regression_functions)
        stopCluster(cl)
        saveRDS(models, file=file.path(paths["intermediary"], "regression_models.Rds"))
        saveRDS(train_set, file=file.path(paths["intermediary"], "train_set.Rds"))
        saveRDS(list(timings=timings),file.path(paths["intermediary"],"timings_train.Rds"))
      }
      metadata.cell$train_set <- train_set
      saveRDS(metadata.cell, file = file.path(paths["intermediary"], "fcs_metadata_df.rds"))
      
      ## imputation
      message("\tImputing unmeasured well-specific markers...")
      {
        xp_orig <- xp
        xp <- xp[, chans]
        
        requireNamespace("parallel")
        cl <- parallel::makeCluster(cores)
        
        ## OS related
        {
          if(.Platform$OS.type == "windows") {
            mc.reset.stream <- function() return(invisible(NULL))
          } else {
            mc.reset.stream <- parallel::mc.reset.stream
          }
          mc.reset.stream()
          
          env <- environment()
        }
        
        #impu.training = F: no prediction for training set
        if(impu.training == F){
          ## sub-sampling for prediction
          message("\t\tRandomly drawing events to predict from the test set (if it's been asked)")
          {
            spl <- split(train_set[,1], events.code) #split data by wells, train_set[,1]= a 0, 1, ... vector
            spl <- lapply(
              spl,
              function(x){ #randomly select test samples
                res <- rep(FALSE,length(x))
                w <- x==0 # if the cell is in test pool, then this returns 'true' to w
                ## sample within test
                res[w][sample(seq_len(sum(w)),min(prediction_events_downsampling,sum(w)))] <- TRUE
                res
              }
            )
            
            pred_set <- which(do.call(c,spl)) #combine testing cells from different wells
          }
          
          
          xp <- xp[pred_set,]
          xp_orig <- xp_orig[pred_set, ]
          
          ## not sure what these used for... about parallel computing??
          {
            clusterExport(
              cl,
              c("xp","chans","regression_functions","polynomial_formula"),
              # c("xp","chans","neural_networks_seed", "regression_functions", "fitter_nn"),
              envir=env
            )
            invisible(clusterEvalQ(
              cl,
              {
                colnames(xp) <- make.names(colnames(xp))
                xp <- xp[,make.names(chans)]
                ## if(any(sapply(regression_functions, function(x){identical(x, fitter_nn)}))){
                ##     if(requireNamespace("keras", quietly = TRUE) & requireNamespace("tensorflow", quietly = TRUE)){
                ##         if(!is.null(neural_networks_seed)){
                ##             tensorflow::use_session_with_seed(neural_networks_seed) ## This will make results reproducible, disable GPU and CPU parallelism (which is good actually). Source: https://keras.rstudio.com/articles/faq.html#how-can-i-obtain-reproducible-results-using-keras-during-development
                ##         }  else {
                ##             tensorflow::tf$reset_default_graph()
                ##             config <- list()
                ##             config$intra_op_parallelism_threads <- 1L
                ##             config$inter_op_parallelism_threads <- 1L
                ##             session_conf <- do.call(tensorflow::tf$ConfigProto, config)
                ##             sess <- tensorflow::tf$Session(graph = tensorflow::tf$get_default_graph(), config = session_conf)
                ##             tensorflow:::call_hook("tensorflow.on_use_session", sess, TRUE)
                ##         }
                ##     }
                ## }
              }
            ))
          }
          
          for(i in seq_along(models)){
            models[[i]] <- lapply(models[[i]],"[[",2)
          }
          
          ## imputing the testing set
          message("\t\tImputing...")
          {
            preds <- list()
            timings <- numeric()
            fun <- predict_wrapper
            environment(fun) <- environment() ## Fixing issue with scoping when cores = 1L
            for(i in seq_along(models)){
              cat("\t\t",names(models)[i],"\n\n",sep="")
              t0 <- Sys.time()
              preds[[i]] <- do.call( #length of list depends on types of models assigned by user
                cbind,
                ##parLapplyLB(
                pblapply(
                  cl=cl,# no. core
                  X=models[[i]], # used model
                  FUN=fun # function do prediction
                )
              )
              t1 <- Sys.time()
              dt <- difftime(t1,t0,units="secs")
              cat("\t",dt," seconds","\n",sep="")
              timings <- c(timings,dt)
              colnames(preds[[i]]) <- names(models[[i]])
            }
            stopCluster(cl)
          }
          
          
          message("\t\tConcatenating predictions...")
          {
            preds <- lapply(preds,function(x){cbind(xp_orig,x)}) #cbind(xp_orig: original in data; x: imputed variables)
            preds <- lapply(preds, as.matrix)
            names(preds) <- names(models)
          }
          
          message("\t\tWriting to disk...")
          {
            saveRDS(preds, file=file.path(paths["downstream"],"predictions.Rds"))
            saveRDS(pred_set, file=file.path(paths["intermediary"],"sampling_preds.Rds")) #should be the same as from inflow.logi.pipe
            saveRDS(list(timings=timings),file.path(paths["intermediary"],"timings_pred.Rds"))
          }
        }else{
          ##impu training
          {
            ## sub-sampling for prediction
            message("\t\tRandomly drawing events to predict from the test set (if it's been asked)")
            {
              spl <- split(train_set[,1], events.code) #split data by wells, train_set[,1]= a 0, 1, ... vector
              spl <- lapply(
                spl,
                function(x){ #randomly select test samples
                  res <- rep(FALSE,length(x))
                  w <- x==1 # if the cell is in test pool, then this returns 'true' to w
                  ## sample within test
                  res[w][sample(seq_len(sum(w)),min(prediction_events_downsampling,sum(w)))] <- TRUE
                  res
                }
              )
              
              pred_set <- which(do.call(c,spl)) #combine testing cells from different wells
            }
            
            
            xp <- xp[pred_set,]
            xp_orig <- xp_orig[pred_set, ]
            
            ## not sure what these used for... about parallel computing??
            {
              clusterExport(
                cl,
                c("xp","chans","regression_functions","polynomial_formula"),
                # c("xp","chans","neural_networks_seed", "regression_functions", "fitter_nn"),
                envir=env
              )
              invisible(clusterEvalQ(
                cl,
                {
                  colnames(xp) <- make.names(colnames(xp))
                  xp <- xp[,make.names(chans)]
                  ## if(any(sapply(regression_functions, function(x){identical(x, fitter_nn)}))){
                  ##     if(requireNamespace("keras", quietly = TRUE) & requireNamespace("tensorflow", quietly = TRUE)){
                  ##         if(!is.null(neural_networks_seed)){
                  ##             tensorflow::use_session_with_seed(neural_networks_seed) ## This will make results reproducible, disable GPU and CPU parallelism (which is good actually). Source: https://keras.rstudio.com/articles/faq.html#how-can-i-obtain-reproducible-results-using-keras-during-development
                  ##         }  else {
                  ##             tensorflow::tf$reset_default_graph()
                  ##             config <- list()
                  ##             config$intra_op_parallelism_threads <- 1L
                  ##             config$inter_op_parallelism_threads <- 1L
                  ##             session_conf <- do.call(tensorflow::tf$ConfigProto, config)
                  ##             sess <- tensorflow::tf$Session(graph = tensorflow::tf$get_default_graph(), config = session_conf)
                  ##             tensorflow:::call_hook("tensorflow.on_use_session", sess, TRUE)
                  ##         }
                  ##     }
                  ## }
                }
              ))
            }
            
            for(i in seq_along(models)){
              models[[i]] <- lapply(models[[i]],"[[",2)
            }
            
            ## imputing the testing set
            message("\t\tImputing...")
            {
              preds <- list()
              timings <- numeric()
              fun <- predict_wrapper
              environment(fun) <- environment() ## Fixing issue with scoping when cores = 1L
              for(i in seq_along(models)){
                cat("\t\t",names(models)[i],"\n\n",sep="")
                t0 <- Sys.time()
                preds[[i]] <- do.call( #length of list depends on types of models assigned by user
                  cbind,
                  ##parLapplyLB(
                  pblapply(
                    cl=cl,# no. core
                    X=models[[i]], # used model
                    FUN=fun # function do prediction
                  )
                )
                t1 <- Sys.time()
                dt <- difftime(t1,t0,units="secs")
                cat("\t",dt," seconds","\n",sep="")
                timings <- c(timings,dt)
                colnames(preds[[i]]) <- names(models[[i]])
              }
              stopCluster(cl)
            }
            
            
            message("\t\tConcatenating predictions...")
            {
              preds <- lapply(preds,function(x){cbind(xp_orig,x)}) #cbind(xp_orig: original in data; x: imputed variables)
              preds <- lapply(preds, as.matrix)
              names(preds) <- names(models)
            }
            
            message("\t\tWriting to disk...")
            {
              saveRDS(preds, file=file.path(paths["downstream"],"predictions_training.Rds"))
              saveRDS(pred_set, file=file.path(paths["intermediary"],"sampling_preds_training.Rds")) #should be the same as from inflow.logi.pipe
              saveRDS(list(timings=timings),file.path(paths["intermediary"],"timings_pred_training.Rds"))
            }
          }
          ##impu testing
          {
            ## sub-sampling for prediction
            message("\t\tRandomly drawing events to predict from the test set (if it's been asked)")
            {
              spl <- split(train_set[,1], events.code) #split data by wells, train_set[,1]= a 0, 1, ... vector
              spl <- lapply(
                spl,
                function(x){ #randomly select test samples
                  res <- rep(FALSE,length(x))
                  w <- x==0 # if the cell is in test pool, then this returns 'true' to w
                  ## sample within test
                  res[w][sample(seq_len(sum(w)),min(prediction_events_downsampling,sum(w)))] <- TRUE
                  res
                }
              )
              
              pred_set <- which(do.call(c,spl)) #combine testing cells from different wells
            }
            
            
            xp <- xp[pred_set,]
            xp_orig <- xp_orig[pred_set, ]
            
            ## not sure what these used for... about parallel computing??
            {
              clusterExport(
                cl,
                c("xp","chans","regression_functions","polynomial_formula"),
                # c("xp","chans","neural_networks_seed", "regression_functions", "fitter_nn"),
                envir=env
              )
              invisible(clusterEvalQ(
                cl,
                {
                  colnames(xp) <- make.names(colnames(xp))
                  xp <- xp[,make.names(chans)]
                  ## if(any(sapply(regression_functions, function(x){identical(x, fitter_nn)}))){
                  ##     if(requireNamespace("keras", quietly = TRUE) & requireNamespace("tensorflow", quietly = TRUE)){
                  ##         if(!is.null(neural_networks_seed)){
                  ##             tensorflow::use_session_with_seed(neural_networks_seed) ## This will make results reproducible, disable GPU and CPU parallelism (which is good actually). Source: https://keras.rstudio.com/articles/faq.html#how-can-i-obtain-reproducible-results-using-keras-during-development
                  ##         }  else {
                  ##             tensorflow::tf$reset_default_graph()
                  ##             config <- list()
                  ##             config$intra_op_parallelism_threads <- 1L
                  ##             config$inter_op_parallelism_threads <- 1L
                  ##             session_conf <- do.call(tensorflow::tf$ConfigProto, config)
                  ##             sess <- tensorflow::tf$Session(graph = tensorflow::tf$get_default_graph(), config = session_conf)
                  ##             tensorflow:::call_hook("tensorflow.on_use_session", sess, TRUE)
                  ##         }
                  ##     }
                  ## }
                }
              ))
            }
            
            for(i in seq_along(models)){
              models[[i]] <- lapply(models[[i]],"[[",2)
            }
            
            ## imputing the testing set
            message("\t\tImputing...")
            {
              preds <- list()
              timings <- numeric()
              fun <- predict_wrapper
              environment(fun) <- environment() ## Fixing issue with scoping when cores = 1L
              for(i in seq_along(models)){
                cat("\t\t",names(models)[i],"\n\n",sep="")
                t0 <- Sys.time()
                preds[[i]] <- do.call( #length of list depends on types of models assigned by user
                  cbind,
                  ##parLapplyLB(
                  pblapply(
                    cl=cl,# no. core
                    X=models[[i]], # used model
                    FUN=fun # function do prediction
                  )
                )
                t1 <- Sys.time()
                dt <- difftime(t1,t0,units="secs")
                cat("\t",dt," seconds","\n",sep="")
                timings <- c(timings,dt)
                colnames(preds[[i]]) <- names(models[[i]])
              }
              stopCluster(cl)
            }
            
            
            message("\t\tConcatenating predictions...")
            {
              preds <- lapply(preds,function(x){cbind(xp_orig,x)}) #cbind(xp_orig: original in data; x: imputed variables)
              preds <- lapply(preds, as.matrix)
              names(preds) <- names(models)
            }
            
            message("\t\tWriting to disk...")
            {
              saveRDS(preds, file=file.path(paths["downstream"],"predictions.Rds"))
              saveRDS(pred_set, file=file.path(paths["intermediary"],"sampling_preds.Rds")) #should be the same as from inflow.logi.pipe
              saveRDS(list(timings=timings),file.path(paths["intermediary"],"timings_pred.Rds"))
            }
          }
        }
        
        message("\t\tVisualising the accuracy of the predictions... (testing set)")
        {
          rsq.ls <- list()
          for(i in 1:length(preds)){
            preds.here <- preds[[i]]
            dat.for.rsq <- data.frame(Well.lab=metadata.cell$Well.lab[which(metadata.cell$train_set == 0)],
                                      train_set=metadata.cell$train_set[which(metadata.cell$train_set == 0)],
                                      preds.here[,-match(chans, colnames(preds.here))])
            uniq.well <- unique(dat.for.rsq$Well.lab)
            rsq.v <- c()
            for(j in 1:length(uniq.well)){
              dat.here <- dat.for.rsq[which(dat.for.rsq$Well.lab == uniq.well[j]),]
              obs <- dat.here$Legend
              pred <- dat.here[,(3+j)]
              
              residual <- (obs - pred)
              
              SSE <- sum(residual^2)
              SSTO <- sum((obs - mean(obs))^2)
              
              rsq.v[j] <- (1-(SSE/SSTO))
            }
            names(rsq.v) <- uniq.well
            rsq.ls[[i]] <- rsq.v
          }
          names(rsq.ls) <- names(preds)
          saveRDS(rsq.ls, file = file.path(paths["intermediary"], "Rsq.Rds"))
          
          ## plotting r-sq
          rsq.df <- do.call(cbind, rsq.ls)
          {
            if(i > 1){
              graph.dat <- melt(data.frame(rsq.df))
              colnames(graph.dat) <- c("Model", "Rsq")
            }else{
              graph.dat <- data.frame(rsq.df, Model=rep(colnames(rsq.df), nrow(rsq.df)))
              colnames(graph.dat)[1] <- "Rsq"
            }
            jpeg(filename = file.path(paths["graph"], paste0("Rsq_boxplot_", length(rsq.ls),"model_500x800.jpeg")), height = 500, width = 800, res = 80)
            p <- ggplot(graph.dat, aes(x=Model, y=Rsq)) +
              # ggtitle("Becht et al.'s Data") +
              geom_boxplot() +
              # scale_fill_manual(values=col.op) +
              coord_flip() +
              theme(legend.position="top",
                    axis.text=element_text(size=28),
                    axis.title=element_text(size=32,face="bold"),
                    legend.key.size = unit(2, 'cm'), #change legend key size
                    legend.key.height = unit(2, 'cm'), #change legend key height
                    legend.key.width = unit(2, 'cm'), #change legend key width
                    legend.title = element_text(size=15, face="bold"), #change legend title font size
                    legend.text = element_text(size=14.5, face="bold"))
            print(p)
            dev.off()
          }
        }
      }
    }
