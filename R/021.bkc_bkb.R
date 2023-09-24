#' Background correction for the backbone markers
#' 
#' This function has been designed to do background correction for the backbone markers by using normal-exponential convolution model.
#' 
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param bkb.v A vector of the names of the backbone markers (MUST match to the names in the FCS file).
#' @param bkb.upper.quntile The cut-off (default = 0.9) for selecting cells used for estimating the parameter of signal.
#' @param bkb.lower.quntile The cut-off (default = 0.1) for selecting cells used for estimating the parameters of noise.
#' @param bkb.min.quntile The cut-off (default = 0.01) for omitting the cells with the smallest values to minimise the impact of outliers.
#' 
#' @author Hsiao-Chi Liao and Agus Salim
#' 
#' @return Generating the calibrated measurements and save to medpara_bkc.bkb_no.bkcPhy_mt.rds file. Visualising the result with the scatter plots.
#' 
bkc_bkb <-
function(
    paths,
    bkb.v,
    bkb.upper.quntile=0.9, #cells used for estimating parameter of signal
    bkb.lower.quntile=0.1, #cells used for estimating parameters of noise
    bkb.min.quntile=0.01 #the lowest 1% of values will not be used to minimise the impact of outliers on sig
    ){
      ## for later use!
      upper.quntile <- bkb.upper.quntile
      lower.quntile <- bkb.lower.quntile
      min.quntile <- bkb.min.quntile
      
      # RDSpath="/Users/chelsea/Desktop/PhD/ProteinExpression/InFlow/Programming/infinity/R_CODE/72.mk_pkg/data/rds/"
      # bkb.v <- c("FSC-H", "FSC-W", "SSC-H", "SSC-W",
      #            "CD69-CD301b", "MHCII", "CD4", "CD44", "CD8",
      #            "CD11c", "CD11b", "F480", "Ly6C", "Lineage", "CD45a488",
      #            "CD24", "CD103")
      # upper.quntile=0.9
      # lower.quntile=0.1
      # min.quntile=0.01 #the lowest 1% of values will not be used to minimize the impact of outliers on sig
      # Graphpath="/Users/chelsea/Desktop/PhD/ProteinExpression/InFlow/Programming/infinity/R_CODE/72.mk_pkg/data/graph/bkc/"
      
      rawInten <- readRDS(file = file.path(paths["intermediary"], "fcs_rawInten_mt.rds"))
      metadata.cell <- readRDS(file = file.path(paths["intermediary"], "fcs_metadata_df.rds"))
      sel.bkb.raw <- rawInten[,match(bkb.v, colnames(rawInten))]
      
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
      
      ## estimating parameters (estimate phy's para but don't use them)
      {
        message("\tEstimating parameters for calibration...")
        para.ls <- list()
        for(m in 1:length(bkb.v)){
          para.df <- data.frame(Well.lab = unique(metadata.cell$Well.lab), mu = NA, sig = NA, alpha = NA)
          for(w in 1:length(unique(metadata.cell$Well.lab))){
            dat.here <- sel.bkb.raw[which(metadata.cell$Well.lab == para.df$Well.lab[w]), m]
            ordering.bkb <- dat.here[order(dat.here)] #for finding quantiles
            
            ## estimate mean and sd of noise (10% smallest values)
            {
              first.quntile <- ordering.bkb[1:(length(ordering.bkb)*lower.quntile)]
              # initial values
              ini.mu <- mean(first.quntile)
              ini.log.sd <- log(sd(first.quntile))
              p.init <- c(ini.mu, ini.log.sd)
              
              out.1 <- tryCatch(optim(p=p.init,
                                      fn=nlogl.norm.v2,
                                      xx=first.quntile,
                                      n=length(ordering.bkb)),
                                error = function(e){NA})
              mle.mean <- out.1$par[1]
              mle.sd <- exp(out.1$par[2])
              
              para.df$mu[w] <- mle.mean
              para.df$sig[w] <- mle.sd
            }
            ## estimate log alpha using tryCatch, so it will not stop if some markers return error
            {
              last.quntile <- ordering.bkb[which(ordering.bkb > quantile(ordering.bkb, probs = upper.quntile))] #or +2*mle.sd
              p.init <- ifelse(mean(last.quntile) < 0, log(0.000001), log(mean(last.quntile)))
              out.2 <- tryCatch(optim(p=p.init,
                                      fn=nlogl.exp,
                                      xx=last.quntile,
                                      n=length(ordering.bkb)),
                                error = function(e){NA})
              # 50: In optim(p = p.init, fn = nlogl.exp, xx = last.0.1, n = length(ordering.bkb)) :
              #   one-dimensional optimization by Nelder-Mead is unreliable:
              #   use "Brent" or optimize() directly
              
              # exp(p[1]) = 1/exponential mean = 1/alpha
              alpha <- (1/exp(out.2$p))
              para.df$alpha[w] <- alpha
            }
            # print(paste0('well=',w,',alpha=',round(alpha[length(alpha)],1)))
          }
          para.ls[[m]] <- para.df
          message(paste0('bkb=',m))
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
          skip.phy.meas <- which(!(colnames(sel.bkb.raw) %in% c("FSC-A", "FSC-H", "FSC-W", 
                                                                "SSC-A", "SSC-H", "SSC-W")))
          for(m in skip.phy.meas){
            ## prepare 'median' parameters
            med.para <- apply(para.ls[[m]][,-1], 2, median)
            med.mle.mean <- med.para[1]
            med.mle.sd <- med.para[2]
            med.alpha <- med.para[3]
            
            calib.temp <- list()
            for(w in 1:length(unique(metadata.cell$Well.lab))){
              dat.here <- sel.bkb.raw[which(metadata.cell$Well.lab == para.df$Well.lab[w]), m]
              
              ## calibrated values
              calib.temp[[w]] <- calib(dat.here, med.alpha, med.mle.mean, sig.2 = med.mle.sd^2)
            }
            #concatenating calib values for each well
            calib.dat[,m] <- do.call(c, calib.temp)
          }
          saveRDS(calib.dat, file = file.path(paths["intermediary"], "medpara_bkc.bkb_no.bkcPhy_mt.rds"))
          
          ## scater plot: raw vs calib
          marker.nam <- colnames(calib.dat)
          for(i in 1:length(marker.nam)){
            jpeg(filename = file.path(paths["graph"], paste0("medpara_sub0.01.", marker.nam[i], "_para.", upper.quntile,"up.", lower.quntile,"lo_500x500.jpeg")), width = 500, height = 500, res = 80)
            #subsample for graphing
            set.seed(123)
            extr.cell <- sample(1:nrow(calib.dat), nrow(calib.dat)*0.01)
            #can consider str. samp in the future
            if(i %in% skip.phy.meas){
              plot(x = sel.bkb.raw[extr.cell,i], y = calib.dat[extr.cell,i],
                   xlab = "Raw", ylab = "Calib.", main = marker.nam[i], col = "blue", 
                   sub = paste0("Parameters Estimated from ", upper.quntile," (hi); ", lower.quntile," (lo)"))
            }else{
              plot(x = sel.bkb.raw[extr.cell,i], y = calib.dat[extr.cell,i],
                   xlab = "Raw", ylab = "Raw.", main = marker.nam[i], col = "blue", 
                   sub = paste0("No need to calibrate the physical measurements."))
            }
            abline(a=0, b=1)
            dev.off()
          }
        }
        message("\tCalibration of backbone markers... Completed!")
      }
    }
