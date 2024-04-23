#' Background correction for the well-specific markers (PE)
#' 
#' This function has been designed to do background correction for the well-specific markers (PE) by using normal-exponential convolution model.
#' 
#' @param paths a vector of characters of paths to store intput, intermediary results, outputs...
#' @param pe.lower.quantile the cut-off (default = 0.1) for selecting cells used for estimating the parameters of noise for infinity markers.
#' @param pe.min.quantile the cut-off (default = 0.01) for omitting the cells with the smallest values to minimise the impact of outliers for infinity markers.
#' @param plots logical; if TRUE (default), produce scatter plots for pre- and post- background adjusted infinity markers (calibrated values on y-axis and raw values on x-axis).
#' 
#' @author Hsiao-Chi Liao and Agus Salim
#' 
#' @importFrom grDevices dev.off jpeg
#' @importFrom graphics abline
#' @importFrom methods is
#' @importFrom stats as.formula contr.sum contrasts<- dexp dnorm median model.matrix optim pexp pnorm quantile sd setNames
#' @importFrom utils head read.csv write.csv
#' 
#' @return Background noise corrected infinity markers and graphs if specified
#' 
#' @details
#' Generating the calibrated measurements and save to bkc.pe_mt.rds file, and visualising the result with the scatter plots in the output directory.
#' 
bkc_pe <-
function(
    paths,
    pe.lower.quantile=0.1, #cells used for estimating parameters of noise
    pe.min.quantile=0.01, #the lowest 1% of values will not be used to minimize the impact of outliers on sig
    plots=TRUE
    ){
    ## for later use!
    min.quantile <- pe.min.quantile
    lower.quantile <- pe.lower.quantile
    
    rawInten <- readRDS(file = file.path(paths["intermediary"], "fcs_rawInten_mt.rds"))
    metadata.cell <- readRDS(file = file.path(paths["intermediary"], "fcs_metadata_df.rds"))
    pe.raw <- as.matrix(rawInten[,"Legend"])
    colnames(pe.raw) <- "Legend"
    
    ## functions for estimating parameters
    {
    # nlogl for mu and sig parameters of noise distribution
    # xx is the observed intensity used for estimating mu and sig
    # the lowest 1% of values will not be used to minimize the impact of outliers on sig
    nlogl.norm.v2 <- function(p, xx, n, rm.min.quantile = min.quantile){
    xx <- sort(xx)[-seq_len(round(n*rm.min.quantile))]
    loglik <- sum(dnorm(xx,mean=p[1],sd=exp(p[2]),log=TRUE)) + 
    (n-length(xx))*pnorm(max(xx),mean=p[1],sd=exp(p[2]),log.p=TRUE,lower.tail=FALSE) +
    round(n*rm.min.quantile)*pnorm(min(xx),mean=p[1],sd=exp(p[2]),log.p=TRUE)
    -loglik
    }
    
    #nlogl for alpha parameter of the true signal
    # xx is the observed intensity used for estimating alpha
    nlogl.exp <- function(p, xx, n) {
    # exp(p[1]) = 1/exponential mean = 1/alpha
    loglik <- sum(dexp(xx,rate=exp(p[1]),log=TRUE)) + (n-length(xx))*pexp(min(xx),rate=exp(p[1]),log.p=TRUE)
    -loglik
    }
    
    # input the obs intensity and return the calibrated intensity
    calib <- function(xx, alpha, mu, sig.2) {
    mu.sx <- xx-mu-sig.2/alpha
    log.2ndterm <- dnorm(0,mean=mu.sx,sd=sqrt(sig.2),log=TRUE)-pnorm(0,mean=mu.sx,sd=sqrt(sig.2),lower.tail=FALSE,log.p=TRUE)
    calib.xx <- mu.sx + sig.2*exp(log.2ndterm)
    calib.xx
    }
    }
    
    ## estimating parameters and calibrating values
    message("\tEstimating parameters for calibration AND calibrating infinity markers...")
    para.df <- data.frame(Well.lab = unique(metadata.cell$Well.lab), mu = NA, sig = NA, alpha = NA)
    
    calib.dat.ls <- list()
    samp.alpha <- c()
    criteria.alpha <- c()
    few.cell.well.ls <- list()
    for(i in seq_along(para.df[,1])){
    well.here <- para.df$Well.lab[i]
    dat.here <- as.matrix(pe.raw[which(metadata.cell$Well.lab == well.here),1])
    
    #ordering according to size
    ordering.legend <- dat.here[,1][order(dat.here[,1])]
    
    ## estimate mean and var of noise (10% smallest values)
    first.quantile <- ordering.legend[seq_len(length(ordering.legend)*lower.quantile)]
    # initial values
    ini.mu <- mean(first.quantile)
    ini.log.sd <- log(sd(first.quantile))
    p.init <- c(ini.mu, ini.log.sd)
    
    suppressWarnings(
    out.1 <- tryCatch(optim(
    par=p.init,
    fn=nlogl.norm.v2,
    xx=first.quantile,
    n=length(ordering.legend)),
    error = function(e){NA})
    )
    
    mle.mean <- out.1$par[1]
    mle.sd <- exp(out.1$par[2])
    
    para.df$mu[i] <- mle.mean
    para.df$sig[i] <- mle.sd
    
    ## estimate log alpha using tryCatch, so it will not stop if some markers return error
    ## orig with large no cell: mle.mean+3*mle.sd -> mle.mean+1*mle.sd
    last.mu.3sd <- ordering.legend[which(ordering.legend > (mle.mean+3*mle.sd))]
    if(length(last.mu.3sd) < 10){
    last.mu.3sd <- ordering.legend[which(ordering.legend >= head(ordering.legend,10))]
    criteria.alpha[i] <- "largest10"
    few.cell.well.ls[[i]] <- well.here
    }else{
    criteria.alpha[i] <- "mean+3sd"
    }
    samp.alpha[i] <- length(last.mu.3sd)
    p.init <- ifelse(mean(last.mu.3sd) < 0, log(0.000001), log(mean(last.mu.3sd)))
    suppressWarnings(
    out.2 <- tryCatch(optim(
    par=p.init,
    fn=nlogl.exp,
    xx=last.mu.3sd,
    n=length(ordering.legend),
    method = "Brent", lower = -20, upper = 0),
    error = function(e){NA})
    )
    
    alpha <- (1/exp(out.2$p))
    para.df$alpha[i] <- alpha
    
    ## calibrating values
    calib.dat <- dat.here
    calib.dat[,1] <- calib(dat.here, alpha, mle.mean, sig.2 = mle.sd^2)
    calib.dat.ls[[i]] <- calib.dat
    
    }
    message("Could not find enough cells (>=10) when used \"mle.mean+3*mle.sd\", so estimated alpha with the top 10 cells with \"the largest values\":\n", length(do.call(c, few.cell.well.ls)), " wells applied this strategy\nSee Wellname_largest10.csv in the intermediary directory for details.")
    out_largest10 <- data.frame(do.call(c, few.cell.well.ls))
    colnames(out_largest10) <- "Well"
    write.csv(out_largest10, file = file.path(paths["intermediary"], "Wellname_largest10.csv"))
    
    para.df <- cbind(para.df, samp.alpha, criteria.alpha)
    saveRDS(para.df, file = file.path(paths["intermediary"], "para.266pe.rds"))
    
    calib.dat <- do.call(rbind, calib.dat.ls); colnames(calib.dat) <- "Legend"
    
    saveRDS(calib.dat, file = file.path(paths["intermediary"], "bkc.pe_mt.rds"))
    
    ## scater plot: raw vs calib
    if(plots==TRUE){
    for(i in seq_along(para.df$Well.lab)){
    jpeg(filename = file.path(paths["graph"], paste0("pe_", para.df$Well.lab[i], "_para.", para.df$criteria.alpha[i], "up.", lower.quantile,"lo_500x500.jpeg")), width = 500, height = 500, res = 80)
    #can consider str. samp in the future
    extr.cell <- which(metadata.cell$Well.lab == para.df$Well.lab[i])
    plot(
    x = pe.raw[extr.cell,1], y = calib.dat[extr.cell,1],
    xlab = "Raw", ylab = "Calib.", main = para.df$Well.lab[i], col = "darkgoldenrod4", 
    sub = paste0("Parameters Estimated from ", criteria.alpha[i], " (", samp.alpha[i]," cells (hi)); ", lower.quantile," (lo)"))
    abline(a=0, b=1)
    dev.off()
    }
    }
    message("\tCalibration of infinity markers... Completed!")
    
    return(calib.dat)
    }
