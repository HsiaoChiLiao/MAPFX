#' Imputing the unmeasured well-specific markers with MICE framework
#'
#' This function has been designed to impute/predict the unmeasured well-specific markers with MICE framework.
#'
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param chans A vector of the names of the backbone markers (MUST match to the names in the FCS file).
#' @param yvar The name of the well-specific marker in the FCS files (has been changed to "Legend" in the first function).
#' @param cores The number of cores used to perform parallel computation (default = 8L).
#' @param iter The number of iterations for MICE to update the imputations
#' @param models.use A vector of the names of the models used for imputation (an example: c("LM", "LASSO3", "SVM", "XGBoost")). The length of the vector is arbitrary.
#' @param extra_args_regression_params A list of the lists of the parameters for running regression models.
#' @param control.wells The well label of the control wells, including the autofluorescence and the isotype controls (format: plate_well, e.g., P1_A01)
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
imputation_mice <-
function(
    paths,
    chans,
    yvar="Legend",
    cores=cores,
    iter=5,
    models.use = models.use,
    extra_args_regression_params = extra_args_regression_params,
    control.wells=control.wells,
    prediction_events_downsampling=NULL
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
      
      # cores=2L
      # control.wells = c("P1_A01", "P2_A01", "P3_A01",
      #                   "P1_A02", "P1_E06", "P1_F01", "P1_G02", "P1_G06", "P1_G08",
      #                   "P2_A12", "P3_B09", "P3_F05", "P3_F08", "P3_F11")
      #
      # paths <- c(intermediary="/Users/chelsea/Desktop/PhD/ProteinExpression/InFlow/Programming/infinity/R_CODE/72.mk_pkg/test_pkg/Output/intermediary/")
      # prediction_events_downsampling=NULL
      
      # RDSpath="/Users/chelsea/Desktop/PhD/ProteinExpression/InFlow/Programming/infinity/R_CODE/72.mk_pkg/data/rds/"
      # Graphpath="/Users/chelsea/Desktop/PhD/ProteinExpression/InFlow/Programming/infinity/R_CODE/72.mk_pkg/data/graph/Imputation/"
      
      chans <- make.names(chans)
      yvar <- make.names(yvar)
      ## impu. posterior dist
      {
        ##model
        bkb.mdl <- readRDS(file = file.path(paths["intermediary"], "regression_models.rds"))
        # bkb.mdl <- readRDS(file = file.path(paths["intermediary"], "logiclez_regression_models.Rds"))
        ##imputation
        bkb.impu <- readRDS(file = file.path(paths["intermediary"], "predictions_training.Rds"))
        # bkb.impu <- readRDS(file = file.path(paths["intermediary"], "logiclez_predictions_train.Rds"))
        bkb.impu <- bkb.impu[[1]] #LM results
        
        bkb.meas <- bkb.impu[,colnames(bkb.impu) %in% chans]; colnames(bkb.meas) <- paste0("I(",colnames(bkb.meas),")")
        
        ##obtain imputation from the posterior dist.
        bkb.impu.post <- bkb.impu
        mean.v <- c()
        sd.v <- c()
        for(i in 1:(ncol(bkb.impu)-length(chans)-1)){
          coef.intcpt <- bkb.mdl[[1]][[i]]$model$coefficients[1]
          coef.v.nointcpt <- bkb.mdl[[1]][[i]]$model$coefficients[-1]
          
          ord.bkb.meas <- bkb.meas[,match(names(coef.v.nointcpt), colnames(bkb.meas))]
          
          print(identical(colnames(ord.bkb.meas), names(coef.v.nointcpt)))
          norm.mean <- ((ord.bkb.meas %*% coef.v.nointcpt) + coef.intcpt)
          norm.sd <- sd(bkb.mdl[[1]][[i]]$model$residuals)
          set.seed(123)
          bkb.impu.post[,((length(chans)+1)+i)] <- rnorm(length(norm.mean), mean = norm.mean, sd = norm.sd)
          # > rnorm(3, mean=c(5,100,1000), sd=1)
          # [1]   4.273491 100.353463 997.852461
          
          mean.v[i] <- norm.mean
          sd.v[i] <- norm.sd
        }
        names(mean.v) = names(sd.v) <- colnames(bkb.impu)[-(1:(length(chans)+1))]
        
        saveRDS(bkb.impu.post, file=file.path(paths["downstream"],"predictions_LMposterior.Rds"))
        save(mean.v, sd.v, file=file.path(paths["downstream"],"FINDpredictions_LMposterior_mean.sd.normal.RData"))
      }
      
      ## doing MICE
      #whole file see: ~/Desktop/PhD/ProteinExpression/InFlow/Programming/infinity/R_CODE/56.MICE_4dat/no_controls_sep21/rev_coherent_train/rev_coherent_train_impu_posterior/pred_obs_testingGoodInitial/final_forPKG/56-4.MICE_rev2.convl.bio.R
      {
        ## impu. posterior dist for training set
        {
          # bkb.impu.post <- readRDS(file=file.path(paths["downstream"],"predictions_LMposterior.Rds"))
        }
        
        # dat <- "rev2.convl.bio"
        # MODEL <- "xgboost"
        {
          
          ## impu from bkb for testing set
          {
            bkb.impu.testing <- readRDS(file=file.path(paths["downstream"],"predictions.Rds"))
            bkb.impu.testing <- bkb.impu.testing[[2]] #xgb results
          }
          
          # inpath="/data/gpfs/projects/punim1447/InFlow/realdat2.66m/MICE_noCntrl/input/"
          # load(file = paste0(inpath, "rev2.convl.bkc.norm.adj_rev.backbone_convl.bkclegend.RData"))
          # [1] "NO.in.all"   "NO.in.well"  "Plate"       "Column"      "Row"         "Well"        "Well.lab"    "train_set"
          # [9] "FSC.H"       "FSC.W"       "SSC.H"       "SSC.W"       "CD69.CD301b" "MHCII"       "CD4"         "CD44"
          # [17] "CD8"         "CD11c"       "CD11b"       "F480"        "Ly6C"        "Lineage"     "CD45a488"    "Legend"
          # [25] "CD24"        "CD103"       "Time"
          # load(file = paste0("/data/gpfs/projects/punim1447/InFlow/realdat2.66m/MICE_noCntrl/input/bkb_impu/meta.pred.train_Set0.RData"))
          # [1] "NO.in.all"  "NO.in.well" "Plate"      "Column"     "Row"        "Well"       "Well.lab"   "train_set"
          ###
          
          # bkb.impu <- readRDS(file = paste0(inpath, "rev2.convl.bkc.norm.adj_rev.backbone_convl.bkclegend_predictions_train.Rds"))
          bkb.impu <- bkb.impu.post
          
          ## re-ordering the data first (this is not a problem for my pkg - have already fix this)
          # dat.here <- rev2.convl.bkc.norm.adj_rev.backbone_convl.bkclegend[with(rev2.convl.bkc.norm.adj_rev.backbone_convl.bkclegend, order(Well.lab, NO.in.all)),]
          # head(unique(dat.here$Well.lab))
          metadata.cell <- readRDS(file = file.path(paths["intermediary"], "fcs_metadata_df.rds"))
          impu.input <- readRDS(file=file.path(paths["downstream"],"impu.input_log.mt.rds"))
          
          dat.here <- data.frame(metadata.cell, impu.input)
          ## check if measurement and pred.bkb are matched
          # identical(as.numeric(bkb.impu[,1]), dat.here[which(dat.here$train_set == 0), 9]) #FSC.H
          # identical(as.numeric(bkb.impu[,5]), dat.here[which(dat.here$train_set == 0), 13]) #cd69.cd301b
          
          dat.here.no.ctrl <- dat.here[-which(dat.here$Well.lab %in% control.wells),]
          bkb.mt <- dat.here.no.ctrl[,which(colnames(dat.here.no.ctrl) %in% chans)]
          well.ord <- unique(dat.here.no.ctrl$Well.lab)
          legend.ls <- split(dat.here.no.ctrl$Legend, dat.here.no.ctrl$Well.lab) #I've renamed it for every input dataset
          head(names(legend.ls))
          
          ## no row and column controls
          ## training set
          bkb.impu.here.no.ctrl <- bkb.impu[-which(dat.here[which(dat.here$train_set == 1),]$Well.lab %in% control.wells), #training set only (rows)
                                            -which(colnames(bkb.impu) %in% control.wells)] #removing controls (cols)
          ## testing set
          bkb.impu.testing.here.no.ctrl <- bkb.impu.testing[-which(dat.here[which(dat.here$train_set == 0),]$Well.lab %in% control.wells), #testing set only (rows)
                                                            -which(colnames(bkb.impu.testing) %in% control.wells)] #removing controls (cols)
          
          training.for.init.aug <- bkb.impu.here.no.ctrl[,-which(colnames(bkb.impu.here.no.ctrl) %in% c(chans, yvar))]
          testing.for.init.aug <- bkb.impu.testing.here.no.ctrl[,-which(colnames(bkb.impu.testing.here.no.ctrl) %in% c(chans, yvar))]
          
          inf.mt <- matrix(NA,
                           nrow = (nrow(training.for.init.aug)+nrow(testing.for.init.aug)),
                           ncol = ncol(training.for.init.aug))
          colnames(inf.mt) <- names(legend.ls)
          
          inf.mt[which(dat.here.no.ctrl$train_set == 1),] <- training.for.init.aug
          inf.mt[which(dat.here.no.ctrl$train_set == 0),] <- testing.for.init.aug
          
          # prev not efficient one
          {
            # # initial: measurements, impu from bkb for training set, mean for testing set
            # inf.mt <- matrix(rep(do.call(c, lapply(legend.ls, mean)), each = 2520000), nrow = 2520000, ncol = 252)
            # colnames(inf.mt) <- names(legend.ls)
            # start.row <- seq(from = 1, to = 2510001, by = 10000)
            # end.row <- seq(from = 10000, to = 2520000, by = 10000)
            # for(i in 1:252){
            #   inf.mt[start.row[i]:end.row[i], i] <- legend.ls[[i]]
            # }
            # identical(colnames(inf.mt), names(legend.ls))
            # ## check data order & column order
            # identical(dat.here.no.ctrl$FSC.H[which(dat.here.no.ctrl$train_set == 1)], as.numeric(bkb.impu.here.no.ctrl[,1]))
            # identical(colnames(bkb.impu.here.no.ctrl)[-(1:18)], colnames(inf.mt))
            # meta.252markers <- dat.here[-which(dat.here$Well.lab %in% c("P1_A01", "P2_A01", "P3_A01",
            #                                                             "P3_F04", "P3_F05", "P3_F06", "P3_F07", "P3_F08",
            #                                                             "P3_F09", "P3_F10", "P3_F11", "P3_F12", "P3_G01", "P3_G02")),
            #                             1:8]
            # #impu from bkb for training set
            # start.row.impu.bkb <- seq(from = 1, to = 1255001, by = 5000)
            # end.row.impu.bkb <- seq(from = 5000, to = 1260000, by = 5000)
            # for(i in 1:252){
            #   meas.rows <- -(start.row[i]:end.row[i]) ##keep measurements
            #   bkb.impu.testing.rows <- -which(meta.252markers$train_set == 0) #skip testing cells
            #   drop.rows <- sort(unique(c(meas.rows, bkb.impu.testing.rows)), decreasing = T)
            #   ##insert into bkb impu rows (not for prev.testing & measurements)
            #   inf.mt[drop.rows, i] <- bkb.impu.here.no.ctrl[-(start.row.impu.bkb[i]:end.row.impu.bkb[i]), (18+i)]
            #   #not into meas well and testing cells <- not using bkb.impu as initial values as there're measurements
            # }
            # #for testing set
            # for(i in 1:252){
            #   meas.rows <- -(start.row[i]:end.row[i]) ##keep measurements
            #   bkb.impu.training.rows <- -which(meta.252markers$train_set == 1) #skip training cells
            #   drop.rows <- sort(unique(c(meas.rows, bkb.impu.training.rows)), decreasing = T)
            #   ##insert into bkb impu rows (not for prev.testing & measurements)
            #   inf.mt[drop.rows, i] <- bkb.impu.testing.here.no.ctrl[-(start.row.impu.bkb[i]:end.row.impu.bkb[i]), (18+i)]
            #   #not into meas well and testing cells <- not using bkb.impu as initial values as there're measurements
            # }
          }
          
          bkb.mt <- as.matrix(bkb.mt)
          save(bkb.mt, inf.mt, file = file.path(paths["intermediary"], "MICE_bkb_init.aug.RData"))
          
          inf.nam <- colnames(inf.mt)
          
          # for(m in 1:seq_along(regression_functions)){ # of models
          
          ## borrowing functions from inflow pkg ##
          ## fitters (only use xgb)
          {
            ## xgboost fitter
            fitter_xgboost <- function(x = NULL, params = NULL, chans = NULL, yvar = NULL){
              #training in one well
              #must extract training cells inside this fucntion and not taken legend for training
              if(!requireNamespace("xgboost", quietly = TRUE)){
                stop("Please run install.packages(\"xgboost\")")
              }
              
              if(!is.null(x) & !is.null(params)){
                args <- c(list(data = x[, chans], label = x[, yvar], nthread = 1L, verbose = 0), params)
                model <- do.call(xgboost::xgboost, args)
                pred <- predict(model, x[, chans])
                rm(list=setdiff(ls(),c("pred","model")))
                return(list(pred=pred,model=model))
                # return(list(model=model))
              }
            }
            
            #below would re-write the input
            regression_functions <- list(XGBoost = fitter_xgboost)
            extra_args_regression_params = list(list(nrounds = 1500, eta = 0.03))
          }
          ## predict_wrapper
          {
            ## xp is the matrix of predictors, x is the built model
            predict_wrapper <- function(x, xp = NULL){
              if(is(x, "lm")){
                xp <- as.data.frame(xp)
              }
              if(is(x, "cv.glmnet")){
                requireNamespace("glmnetUtils")
                xp <- as.data.frame(xp)
                return(predict(x,xp,s=x$lambda.min)[,1])
              }
              if(is(x, "raw")){
                requireNamespace("keras")
                x = keras::unserialize_model(x)
              }
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
          
          ## selecting cells for training
          extr.samp.all <- which(dat.here.no.ctrl$train_set == 1)
          well.train <- dat.here.no.ctrl[,which(colnames(dat.here.no.ctrl) %in% c("Well.lab","train_set"))]
          
          train.rsq.cy = test.rsq.cy <- list()
          train.pred.ls.cy = train.obs.ls.cy = test.pred.ls.cy = test.obs.ls.cy <- list()
          aug.mt <- cbind(inf.mt, bkb.mt)
          start.row <- as.integer(which(!duplicated(dat.here.no.ctrl$Well.lab))) #as.integer(seq(from = 1, to = 2510001, by = 10000))
          end.row <- as.integer(c(start.row-1, length(dat.here.no.ctrl$Well.lab))) #as.integer(seq(from = 10000, to = 2520000, by = 10000))
          end.row <- end.row[-1]
          ##
          
          well <- unique(well.train$Well.lab) # start.row.impu.bkb <- seq(from = 1, to = 1255001, by = 5000); end.row.impu.bkb <- seq(from = 5000, to = 1260000, by = 5000)
          
          # coef.ls.ls <- list()
          for(c in 1:iter){ # of cycle
            train.rsq.v = test.rsq.v <- c()
            train.pred.ls = train.obs.ls = test.pred.ls = test.obs.ls <- list()
            A <- Sys.time()
            # coef.ls <- list()
            for(i in 1:ncol(inf.mt)){ # of inf.
              # spec.inf <- aug.mt[start.row[i]:end.row[i],] #all data
              
              extr.samp <- which((well.train$Well.lab == well[i]) & (well.train$train_set == 1)) #extr.samp.all[start.row.impu.bkb[i]:end.row.impu.bkb[i]]
              training.dat <- aug.mt[extr.samp,]
              this.win <- (start.row[i]:end.row[i])
              # testing.dat.w.obs <- aug.mt[this.win[-which(this.win %in% extr.samp)],] #'pos' for calculating r-sq (don't use dat from aug.mt as that wasn't true values)
              
              {
                ## set up yvar and chans~preds for building model and making prediction
                ## predicting process
                {
                  # markers as predictors
                  chans <- colnames(aug.mt)[-i]
                  names(chans) <- chans
                  
                  ## parameters
                  {
                    # paths="~/trained_model/"
                    # cores=8L
                    chans=make.names(chans)
                    yvar=make.names(colnames(training.dat)[i])
                    xp=aug.mt[,-which(colnames(aug.mt) == yvar)]
                    
                    verbose=TRUE
                  }
                  
                  if(verbose){
                    message("Fitting regression models")
                  }
                  
                  # timings <- numeric()
                  cat("\t\t", names(regression_functions)[1], "\n\n", sep="")
                  t0 <- Sys.time()
                  fun <- regression_functions[[1]] #xgb
                  built.model <- fun(x = training.dat, params = extra_args_regression_params[[1]],
                                     chans = chans, yvar = yvar) #***
                  # if(m == 1){ #for LM
                  #   coef.ls[[i]] <- built.model$model$coefficients
                  # }
                  t1 <- Sys.time()
                  dt <- difftime(t1,t0,units="secs")
                  cat("\t",dt," seconds","\n",sep="")
                  # timings <- c(timings,dt)
                  
                  if(verbose){
                    message("\tImputing...")
                  }
                  
                  print(head(xp))
                  # timings <- numeric()
                  t0 <- Sys.time()
                  Ypred <- predict_wrapper(x = built.model[[2]], xp = xp) #input aug matrix, output should be the target column (length = 2.66m)
                  t1 <- Sys.time()
                  dt <- difftime(t1,t0,units="secs")
                  cat("\t",dt," seconds","\n",sep="")
                  # timings <- c(timings,dt)
                }
              }
              
              ## updating the aug.mt
              aug.mt[,i] <- Ypred
              aug.mt[start.row[i]:end.row[i], i] <- legend.ls[[i]] #put measurement back!
              
              # calculating r-sq with 5,000 cells which have observation in testing dataset
              {
                ## extracting prediction for testing set for calculating r-sq (oj)
                ##using testing set only
                this.win <- (start.row[i]:end.row[i]) #size=10,000 /// to 2520000
                ## test-rsq ##
                test.pred.Y <- Ypred[this.win[-which(this.win %in% extr.samp)]] #not using cells for training (ps: train_set=1 is the training here)
                test.true.Y <- dat.here.no.ctrl$Legend[this.win[-which(this.win %in% extr.samp)]] #not using cells for training
                
                test.pred.ls[[i]] <- test.pred.Y
                test.obs.ls[[i]] <- test.true.Y
                
                test.residual <- (test.true.Y - test.pred.Y)
                
                test.SSE <- sum(test.residual^2)
                test.SSTO <- sum((test.true.Y - mean(test.true.Y))^2)
                
                #rsq
                test.rsq.v[i] <- (1-(test.SSE/test.SSTO))
                ###
                ## train-rsq ##
                train.pred.Y <- Ypred[this.win[which(this.win %in% extr.samp)]] #using cells for training
                train.true.Y <- dat.here.no.ctrl$Legend[this.win[which(this.win %in% extr.samp)]] #using cells for training
                
                train.pred.ls[[i]] <- train.pred.Y
                train.obs.ls[[i]] <- train.true.Y
                
                train.residual <- (train.true.Y - train.pred.Y)
                
                train.SSE <- sum(train.residual^2)
                train.SSTO <- sum((train.true.Y - mean(train.true.Y))^2)
                
                #rsq
                train.rsq.v[i] <- (1-(train.SSE/train.SSTO))
              }
            }
            # if(m == 1){
            #   coef.ls.ls[[c]] <- coef.ls
            # }
            B <- Sys.time()
            print(paste0("model: ", names(regression_functions)[1]))
            print(paste0("cycle: ", c))
            B-A
            names(train.rsq.v) = names(test.rsq.v) <- inf.nam
            test.rsq.cy[[c]] <- test.rsq.v
            train.rsq.cy[[c]] <- train.rsq.v
            
            train.pred.ls.cy[[c]] <- train.pred.ls
            train.obs.ls.cy[[c]] <- train.obs.ls
            test.pred.ls.cy[[c]] <- test.pred.ls
            test.obs.ls.cy[[c]] <- test.obs.ls
            
            ## output - each cycle
            if(c == iter){
              saveRDS(aug.mt, file=file.path(paths["downstream"], "MICE_XGBimpu_aug_FINALcyc.", c, ".RData"))
            }else{
              saveRDS(aug.mt, file=file.path(paths["intermediary"], "MICE_XGBimpu_aug_cyc.", c, ".RData"))
            }
          }
          names(train.rsq.cy) = names(test.rsq.cy) <- paste0("cycle",1:iter)
          # rsq.md[[m]] <- rsq.cy
          save(train.rsq.cy, test.rsq.cy, file=file.path(paths["intermediary"], "MICE_XGBimpu_rsq.RData"))
          save(train.pred.ls.cy, train.obs.ls.cy, test.pred.ls.cy, test.obs.ls.cy, file=file.path(paths["intermediary"], "MICE_XGBimpu_obs.pred.RData"))
          # if(m == 1){
          #   save(coef.ls.ls, file = paste0(outpath, "/coefLM/posterior/", dat, "_mdl.", names(regression_functions)[m], "_impu.MICE.RData"))
          # }
        }
      }
    }
