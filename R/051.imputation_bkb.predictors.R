#' Imputing the unmeasured well-specific markers with regression models
#'
#' This function has been designed to impute/predict the unmeasured well-specific markers with regression models.
#'
#' @param paths a vector of characters of paths to store intput, intermediary results, outputs...
#' @param chans a vector of the names of the backbone markers (MUST match to the names in the FCS file).
#' @param yvar the name of the well-specific marker in the FCS files (has been changed to "Legend" in the first function).
#' @param cores the number of cores used to perform parallel computation (default = 8L).
#' @param models.use a vector of the names of the models used for imputation (an example: c("LM", "LASSO3", "SVM", "XGBoost")). The length of the vector is arbitrary.
#' @param extra_args_regression_params a list of the lists of the parameters for running regression models.
#' @param prediction_events_downsampling default = NULL (not doing subsampling). How many cells per well you want to have the imputation? (must be less than or equal to a half as we won't get the prediction for cells in the training set).
#' @param impu.training logical; if FALSE (default), not impute the training set (the dataset used to train the imputation models).
#' @param plots logical; if TRUE (default), visualise the distribution of R-sq from each infinity marker.
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
#' @import cowplot
#' @importFrom reshape2 melt
#' @importFrom gtools combinations
#' 
#' @importFrom grDevices dev.off jpeg
#' @importFrom graphics abline
#' @importFrom methods is
#' @importFrom stats as.formula contr.sum contrasts<- dexp dnorm median model.matrix optim pexp pnorm quantile sd setNames
#' @importFrom utils head read.csv
#'
#' @return Imputing the unmeasured well-specific markers and save to predictions.Rds file. Visualising the result with the boxplots (r-sq).
#'
imputation_bkb.predictors <-
function(
    paths,
    chans,
    yvar="Legend",
    cores=4L,
    models.use,
    extra_args_regression_params,
    prediction_events_downsampling=NULL,
    impu.training,
    plots=TRUE
    ){
    
    #binding the variable locally to the function
    Model <- Rsq <- NULL
    
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
    model$model <- NULL ## Trim down for slimmer objects
    model$qr$qr <- NULL ## Trim down for slimmer objects
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
    else{
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
    for(m in seq_len(length(models.use))){
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
    # set.seed(123)
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
    envir=env
    )
    
    clusterEvalQ(
    cl,
    {
    chans <- make.names(chans)
    yvar <- make.names(yvar)
    }
    )
    }
    
    ## fitting models
    message("\t\tFitting...")
    {
    models <- list()
    timings <- numeric()
    
    for(i in seq_along(regression_functions)){
    message("\t\t",names(regression_functions)[i],"\n\n")
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
    message("\t",dt," seconds","\n")
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
    if(impu.training == FALSE){
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
    
    {
    clusterExport(
    cl,
    c("xp","chans","regression_functions","polynomial_formula"),
    envir=env
    )
    invisible(clusterEvalQ(
    cl,
    {
    colnames(xp) <- make.names(colnames(xp))
    xp <- xp[,make.names(chans)]
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
    message("\t\t",names(models)[i],"\n\n")
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
    message("\t",dt," seconds","\n")
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
    message("\t\tRandomly drawing events (if specified) from the training set for prediction.")
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
    envir=env
    )
    invisible(clusterEvalQ(
    cl,
    {
    colnames(xp) <- make.names(colnames(xp))
    xp <- xp[,make.names(chans)]
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
    message("\t\t",names(models)[i],"\n\n")
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
    message("\t",dt," seconds","\n")
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
    message("\t\tRandomly drawing events (if specified) from the testing set for prediction.")
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
    
    {
    clusterExport(
    cl,
    c("xp","chans","regression_functions","polynomial_formula"),
    envir=env
    )
    invisible(clusterEvalQ(
    cl,
    {
    colnames(xp) <- make.names(colnames(xp))
    xp <- xp[,make.names(chans)]
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
    message("\t\t",names(models)[i],"\n\n")
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
    message("\t",dt," seconds","\n")
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
    
    
    {
    rsq.ls <- list()
    for(i in seq_along(preds)){
    preds.here <- preds[[i]]
    dat.for.rsq <- data.frame(
    Well.lab=metadata.cell$Well.lab[which(metadata.cell$train_set == 0)],
    train_set=metadata.cell$train_set[which(metadata.cell$train_set == 0)],
    preds.here[,-match(chans, colnames(preds.here))])
    uniq.well <- unique(dat.for.rsq$Well.lab)
    rsq.v <- c()
    for(j in seq_along(uniq.well)){
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
    
    if(plots == TRUE){
    message("\t\tVisualising the accuracy of the predictions... (using testing set)")
    ## plotting r-sq
    rsq.df <- do.call(cbind, rsq.ls)
    {
    if(length(rsq.ls) > 1){
    graph.dat <- melt(data.frame(rsq.df))
    colnames(graph.dat) <- c("Model", "Rsq")
    }else{
    graph.dat <- data.frame(rsq.df, Model=rep(colnames(rsq.df), nrow(rsq.df)))
    colnames(graph.dat)[1] <- "Rsq"
    }
    
    if(length(rsq.ls) == 1){
    plt1 <- ggplot(data=graph.dat, mapping = aes(x=Model, y=Rsq)) +
    geom_boxplot(fill = "magenta3", color = "black") + 
    coord_flip() +
    xlab("") +
    theme(axis.text.y=element_blank(),
    axis.ticks.y=element_blank()) +
    theme_classic()
    
    plt2 <- ggplot(data=graph.dat, mapping = aes(x=Rsq)) +
    labs(title = "The Performance of Imputation",
    subtitle = paste0("Imputation Models: ", paste(names(rsq.ls), collapse = ","))) +
    geom_histogram(fill = "magenta", color = "black") +
    ylab("Counts") +
    theme_classic()
    
    jpeg(filename = file.path(paths["graph"], paste0("Rsq_boxplot_", length(rsq.ls),"model_500x800.jpeg")), height = 500, width = 800, res = 80)
    p <- plot_grid(plt2, plt1, 
    ncol = 1, rel_heights = c(2, 1),
    align = 'v', axis = 'lr')
    plot(p)
    dev.off()
    }else{ #more than one imputation method (boxplot only)
    jpeg(filename = file.path(paths["graph"], paste0("Rsq_boxplot_", length(rsq.ls),"model_500x800.jpeg")), height = 500, width = 800, res = 80)
    p <- ggplot(data=graph.dat, mapping = aes(x=Model, y=Rsq, fill=Model)) +
    geom_boxplot() +
    coord_flip() +
    labs(title = "The Performance of Imputation",
    subtitle = paste0("Imputation Models: ", paste(names(rsq.ls), collapse = ", "))) +
    theme(legend.position="none",
    plot.title = element_text(size=20),
    plot.subtitle = element_text(size=16),
    axis.text=element_text(size=14, face = "bold"),
    axis.title=element_text(size=18,face="bold"),
    # legend.key.size = unit(1, 'cm'), #change legend key size
    # legend.key.height = unit(1, 'cm'), #change legend key height
    # legend.key.width = unit(1, 'cm'), #change legend key width
    # legend.title = element_text(size=10, face="bold"), #change legend title font size
    # legend.text = element_text(size=8, face="bold")
    )
    plot(p)
    dev.off()
    }
    }
    }
    }
    }
    }
    