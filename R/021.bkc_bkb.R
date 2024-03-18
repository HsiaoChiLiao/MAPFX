#' Background correction for the backbone markers
#' 
#' This function has been designed to do background correction for the backbone markers by using normal-exponential convolution model.
#' 
#' @param paths a vector of characters of paths to store intput, intermediary results, outputs...
#' @param bkb.v a vector of the names of the backbone markers (MUST match to the names in the FCS file).
#' @param bkb.upper.quantile the cut-off (default = 0.9) for selecting cells used for estimating the parameter of signal.
#' @param bkb.lower.quantile the cut-off (default = 0.1) for selecting cells used for estimating the parameters of noise.
#' @param bkb.min.quantile the cut-off (default = 0.01) for omitting the cells with the smallest values to minimise the impact of outliers.
#' @param plots logical; if TRUE (default), produce scatter plots for pre- and post- background adjusted backbone markers (calibrated values on y-axis and raw values on x-axis).
#' 
#' @author Hsiao-Chi Liao and Agus Salim
#' 
#' @importFrom grDevices dev.off jpeg
#' @importFrom graphics abline
#' @importFrom methods is
#' @importFrom stats as.formula contr.sum contrasts<- dexp dnorm median model.matrix optim pexp pnorm quantile sd setNames
#' @importFrom utils head read.csv
#' 
#' @return Generating the calibrated measurements and save to medpara_bkc.bkb_no.bkcPhy_mt.rds file. Visualising the result with the scatter plots.
#' 
bkc_bkb <-
function(
    paths,
    bkb.v,
    bkb.upper.quantile=0.9, #cells used for estimating parameter of signal
    bkb.lower.quantile=0.1, #cells used for estimating parameters of noise
    bkb.min.quantile=0.01, #the lowest 1% of values will not be used to minimise the impact of outliers on sig
    plots
    ){
    ## for later use!
    upper.quantile <- bkb.upper.quantile
    lower.quantile <- bkb.lower.quantile
    min.quantile <- bkb.min.quantile
    
    rawInten <- readRDS(file = file.path(paths["intermediary"], "fcs_rawInten_mt.rds"))
    metadata.cell <- readRDS(file = file.path(paths["intermediary"], "fcs_metadata_df.rds"))
    sel.bkb.raw <- rawInten[,match(bkb.v, colnames(rawInten))]
    
    ## functions for estimating parameters
    {
    # nlogl for mu and sig parameters of noise distribution
    # xx = the observed intensity used for estimating mu and sig
    # the lowest 1% of values will not be used to minimize the impact of outliers on sig
    nlogl.norm.v2 <- function(p, xx, n, rm.min.quantile = min.quantile){
    xx <- sort(xx)[-seq_len(round(n*rm.min.quantile))]
    loglik <- sum(dnorm(xx,mean=p[1],sd=exp(p[2]),log=TRUE)) + 
    (n-length(xx))*pnorm(max(xx),mean=p[1],sd=exp(p[2]),log.p=TRUE,lower.tail=FALSE) +
    round(n*rm.min.quantile)*pnorm(min(xx),mean=p[1],sd=exp(p[2]),log.p=TRUE)
    -loglik
    }
    
    #nlogl for alpha parameter of the true signal
    # xx = the observed intensity used for estimating alpha
    nlogl.exp <- function(p, xx, n){
    # exp(p[1]) = 1/exponential mean = 1/alpha
    loglik <- sum(dexp(xx,rate=exp(p[1]),log=TRUE)) + 
    (n-length(xx))*pexp(min(xx),rate=exp(p[1]),log.p=TRUE)
    -loglik
    }
    
    # input the obs intensity and return the calibrated intensity
    calib <- function(xx, alpha, mu, sig.2){
    mu.sx <- xx-mu-sig.2/alpha
    log.2ndterm <- dnorm(0,mean=mu.sx,sd=sqrt(sig.2),log=TRUE)-pnorm(0,mean=mu.sx,sd=sqrt(sig.2),lower.tail=FALSE,log.p=TRUE)
    calib.xx <- mu.sx + sig.2*exp(log.2ndterm)
    calib.xx
    }
    }
    
    ## estimating parameters (estimate phy's para but don't use them)
    {
    message("\tEstimating parameters for calibration...")
    para.ls <- list()
    for(m in seq_len(length(bkb.v))){
    para.df <- data.frame(Well.lab = unique(metadata.cell$Well.lab), mu = NA, sig = NA, alpha = NA)
    for(w in seq_len(length(unique(metadata.cell$Well.lab)))){
    dat.here <- sel.bkb.raw[which(metadata.cell$Well.lab == para.df$Well.lab[w]), m]
    ordering.bkb <- dat.here[order(dat.here)] #for finding quantiles
    
    ## estimate mean and sd of noise (10% smallest values)
    {
    first.quantile <- ordering.bkb[seq_len(length(ordering.bkb)*lower.quantile)]
    # initial values
    ini.mu <- mean(first.quantile)
    ini.log.sd <- log(sd(first.quantile))
    p.init <- c(ini.mu, ini.log.sd)
    
    suppressWarnings(
    out.1 <- tryCatch(optim(
    par=p.init,
    fn=nlogl.norm.v2,
    xx=first.quantile,
    n=length(ordering.bkb)),
    error = function(e){NA}))
    mle.mean <- out.1$par[1]
    mle.sd <- exp(out.1$par[2])
    
    para.df$mu[w] <- mle.mean
    para.df$sig[w] <- mle.sd
    }
    ## estimate log alpha using tryCatch, so it will not stop if some markers return error
    {
    last.quantile <- ordering.bkb[which(ordering.bkb > quantile(ordering.bkb, probs = upper.quantile))] #or +2*mle.sd
    p.init <- ifelse(mean(last.quantile) < 0, log(0.000001), log(mean(last.quantile)))
    suppressWarnings(
    out.2 <- tryCatch(optim(
    par=p.init,
    fn=nlogl.exp,
    xx=last.quantile,
    n=length(ordering.bkb)),
    error = function(e){NA}))

    # exp(p[1]) = 1/exponential mean = 1/alpha
    alpha <- (1/exp(out.2$p))
    para.df$alpha[w] <- alpha
    }
    # print(paste0('well=',w,',alpha=',round(alpha[length(alpha)],1)))
    }
    para.ls[[m]] <- para.df
    message('bkb=',m)
    }
    names(para.ls) <- colnames(sel.bkb.raw)
    saveRDS(para.ls, file = file.path(paths["intermediary"], paste0("para.ls_",ncol(sel.bkb.raw),"bkb_",length(unique(metadata.cell$Well.lab)),"well.rds")))
    message("\tEstimation of parameters... Completed!")
    }
    
    ## calibrating values
    {
    message("\tCalibrating backbone markers (except for physical measurements)...")
    calib.dat <- sel.bkb.raw
    {
    skip.phy.meas <- which(!(colnames(sel.bkb.raw) %in% c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W")))
    for(m in skip.phy.meas){
    ## prepare 'median' parameters
    med.para <- apply(para.ls[[m]][,-1], 2, median)
    med.mle.mean <- med.para[1]
    med.mle.sd <- med.para[2]
    med.alpha <- med.para[3]
    
    calib.temp <- list()
    for(w in seq_len(length(unique(metadata.cell$Well.lab)))){
    dat.here <- sel.bkb.raw[which(metadata.cell$Well.lab == para.df$Well.lab[w]), m]
    
    ## calibrated values
    calib.temp[[w]] <- calib(dat.here, med.alpha, med.mle.mean, sig.2 = med.mle.sd^2)
    }
    #concatenating calib values for each well
    calib.dat[,m] <- do.call(c, calib.temp)
    }
    saveRDS(calib.dat, file = file.path(paths["intermediary"], "medpara_bkc.bkb_no.bkcPhy_mt.rds"))
    
    ## scater plot: raw vs calib
    if(plots == TRUE){
    marker.nam <- colnames(calib.dat)
    for(i in seq_along(marker.nam)){
    jpeg(filename = file.path(paths["graph"], paste0("medpara_sub0.01.", marker.nam[i], "_para.", upper.quantile,"up.", lower.quantile,"lo_500x500.jpeg")), width = 500, height = 500, res = 80)
    #subsample for graphing
    # set.seed(123)
    extr.cell <- sample(seq_len(nrow(calib.dat)), nrow(calib.dat)*0.01)
    #can consider str. samp in the future
    if(i %in% skip.phy.meas){
    plot(
    x = sel.bkb.raw[extr.cell,i], y = calib.dat[extr.cell,i],
    xlab = "Raw", ylab = "Calib.", main = marker.nam[i], col = "blue",
    sub = paste0("Parameters Estimated from ", upper.quantile," (hi); ", lower.quantile," (lo)"))
    }else{
    plot(
    x = sel.bkb.raw[extr.cell,i], y = calib.dat[extr.cell,i],
    xlab = "Raw", ylab = "Raw.", main = marker.nam[i], col = "blue",
    sub = paste0("No need to calibrate the physical measurements."))
    }
    abline(a=0, b=1)
    dev.off()
    }
    }
    }
    message("\tCalibration of backbone markers... Completed!")
    }
    }
