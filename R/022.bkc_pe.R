#' Background correction for the well-specific markers (PE)
#' 
#' This function has been designed to do background correction for the well-specific markers (PE) by using normal-exponential convolution model.
#' 
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param pe.mean.sd Selecting cells with the value larger than mean + ?sd (default = 3) to estimate the parameter of signal.
#' @param pe.lower.quntile The cut-off (default = 0.1) for selecting cells used for estimating the parameters of noise.
#' @param pe.min.quntile The cut-off (default = 0.01) for omitting the cells with the smallest values to minimise the impact of outliers.
#' 
#' @author Hsiao-Chi Liao and Agus Salim
#' 
#' @return Generating the calibrated measurements and save to bkc.pe_mt.rds file. Visualising the result with the scatter plots.
#' 
bkc_pe <-
function(
    paths,
    pe.mean.sd=3, #mean+-3sd for extracting extremely large values
    # upper.quntile=0.9, #cells used for estimating parameter of signal
    pe.lower.quntile=0.1, #cells used for estimating parameters of noise
    pe.min.quntile=0.01 #the lowest 1% of values will not be used to minimize the impact of outliers on sig
    ){
      ## for later use!
      min.quntile <- pe.min.quntile
      lower.quntile <- pe.lower.quntile
      
      # RDSpath="/Users/chelsea/Desktop/PhD/ProteinExpression/InFlow/Programming/infinity/R_CODE/72.mk_pkg/data/rds/"
      # # upper.quntile=0.9
      # mean.sd=3
      # lower.quntile=0.1
      # min.quntile=0.01 #the lowest 1% of values will not be used to minimize the impact of outliers on sig
      # Graphpath="/Users/chelsea/Desktop/PhD/ProteinExpression/InFlow/Programming/infinity/R_CODE/72.mk_pkg/data/graph/bkc/"
      
      rawInten <- readRDS(file = file.path(paths["intermediary"], "fcs_rawInten_mt.rds"))
      metadata.cell <- readRDS(file = file.path(paths["intermediary"], "fcs_metadata_df.rds"))
      pe.raw <- as.matrix(rawInten[,"Legend"]); colnames(pe.raw) <- "Legend"
      
      ## functions for estimating parameters
      {
        # nlogl for mu and sig parameters of noise distribution
        # xx = the observed intensity used for estimating mu and sig
        # the lowest 1% of values will not be used to minimize the impact of outliers on sig
        nlogl.norm.v2 <- function(p, xx, n, rm.min.quntile = min.quntile) {
          xx <- sort(xx)[-1:-round(n*rm.min.quntile)]
          loglik <- sum(dnorm(xx,mean=p[1],sd=exp(p[2]),log=TRUE)) + (n-length(xx))*pnorm(max(xx),mean=p[1],sd=exp(p[2]),log.p=TRUE,lower.tail=FALSE) +
            round(n*rm.min.quntile)*pnorm(min(xx),mean=p[1],sd=exp(p[2]),log.p=TRUE)
          -loglik
        }
        
        #nlogl for alpha parameter of the true signal
        # xx = the observed intensity used for estimating alpha
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
      message("\tEstimating parameters for calibration AND calibrating PE markers...")
      para.df <- data.frame(Well.lab = unique(metadata.cell$Well.lab), mu = NA, sig = NA, alpha = NA)
      
      calib.dat.ls <- list()
      samp.alpha <- c()
      criteria.alpha <- c()
      for(i in 1:nrow(para.df)){
        well.here <- para.df$Well.lab[i]
        dat.here <- as.matrix(pe.raw[which(metadata.cell$Well.lab == well.here),1])
        
        #ordering according to size
        ordering.legend <- dat.here[,1][order(dat.here[,1])]
        
        ## estimate mean and var of noise (10% smallest values)
        first.quntile <- ordering.legend[1:(length(ordering.legend)*lower.quntile)]
        # initial values
        ini.mu <- mean(first.quntile)
        ini.log.sd <- log(sd(first.quntile))
        p.init <- c(ini.mu, ini.log.sd)
        
        out.1 <- tryCatch(optim(p=p.init,
                                fn=nlogl.norm.v2,
                                xx=first.quntile,
                                n=length(ordering.legend)),
                          error = function(e){NA})
        mle.mean <- out.1$par[1]
        mle.sd <- exp(out.1$par[2])
        
        para.df$mu[i] <- mle.mean
        para.df$sig[i] <- mle.sd
        
        ## estimate log alpha using tryCatch, so it will not stop if some markers return error
        last.mu.3sd <- ordering.legend[which(ordering.legend > (mle.mean+3*mle.sd))]
        if(length(last.mu.3sd) == 0){
          last.mu.3sd <- ordering.legend[which(ordering.legend > (mle.mean+1*mle.sd))]
          message(paste0("problem marker: ", well.here, "; ",
                         "no. cell for alpha (mean+1sd): ", length(last.mu.3sd)))
          criteria.alpha[i] <- "mean+1sd"
        }else{
          criteria.alpha[i] <- "mean+3sd"
        }
        samp.alpha[i] <- length(last.mu.3sd)
        p.init <- ifelse(mean(last.mu.3sd) < 0, log(0.000001), log(mean(last.mu.3sd)))
        out.2 <- tryCatch(optim(p=p.init,
                                fn=nlogl.exp,
                                xx=last.mu.3sd,
                                n=length(ordering.legend)),
                          error = function(e){NA})
        # 50: In optim(p = p.init, fn = nlogl.exp, xx = last.mu.3sd,  ... :
        #                one-dimensional optimization by Nelder-Mead is unreliable:
        #                use "Brent" or optimize() directly
        
        # exp(p[1]) = 1/exponential mean = 1/alpha
        alpha <- (1/exp(out.2$p))
        para.df$alpha[i] <- alpha
        
        ## calibrating values
        calib.dat <- dat.here
        calib.dat[,1] <- calib(dat.here, alpha, mle.mean, sig.2 = mle.sd^2)
        calib.dat.ls[[i]] <- calib.dat
        
        # print(paste0('Iteration=',i,', alpha=',round(alpha[length(alpha)],1)))
      }
      para.df <- cbind(para.df, samp.alpha, criteria.alpha)
      saveRDS(para.df, file = file.path(paths["intermediary"], "para.266pe.rds"))
      
      calib.dat <- do.call(rbind, calib.dat.ls); colnames(calib.dat) <- "Legend"
      
      saveRDS(calib.dat, file = file.path(paths["intermediary"], "bkc.pe_mt.rds"))
      
      ## scater plot: raw vs calib
      for(i in 1:length(para.df$Well.lab)){
        jpeg(filename = file.path(paths["graph"], paste0("pe_", para.df$Well.lab[i], "_para.", para.df$criteria.alpha[i], "up.", lower.quntile,"lo_500x500.jpeg")), width = 500, height = 500, res = 80)
        #can consider str. samp in the future
        extr.cell <- which(metadata.cell$Well.lab == para.df$Well.lab[i])
        plot(x = pe.raw[extr.cell,1], y = calib.dat[extr.cell,1],
             xlab = "Raw", ylab = "Calib.", main = para.df$Well.lab[i], col = "darkgoldenrod4", 
             sub = paste0("Parameters Estimated from ", criteria.alpha[i], " (", samp.alpha[i]," cells (hi)); ", lower.quntile," (lo)"))
        abline(a=0, b=1)
        dev.off()
      }
      message("\tCalibration of Legend markers... Completed!")
    }
