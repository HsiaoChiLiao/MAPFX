library(flowCore)
## Using the Example Datasets in MAPFX Package for the Vignette
### MPC
# This dataset is a subset of the single-cell murine lung data at steady state downloaded from 
# [FlowRepository](https://flowrepository.org/id/FR-FCM-Z2LP) provided by Etienne Becht (Nov 2020). 
# The raw protein intensities and the corresponding metadata were saved in 
# `data/ord.fcs.raw.mt_mpc.rda` and `data/ord.fcs.raw.meta.df.out_mpc.rda` 
# which were generated from 266 .FCS files from 266 wells with 50 cells in each file.

## MPC - generating example data (random 50 per file)
{
  inpath="/path-to-MPCfcs/"
  subsampath="/path-to-subsamp/"
  
  ## Subsampling
  input_events_downsampling=50
  extra_args_read_FCS=NULL
  
  files <- list.files(inpath,full.names=TRUE,recursive=TRUE,pattern=".fcs")
  invisible(
    lapply(
      files,
      function(file){
        res <- do.call(read.FCS,c(list(filename=file),extra_args_read_FCS))
        w <- sort(sample(seq_len(nrow(res)),min(input_events_downsampling,nrow(res))))
        res <- res[w,]
        write.FCS(res,file.path(subsampath, basename(file)))
      }
    )
  )
}

### FFC
# This mice splenocytes dataset contains 50 cells (sorted CD4+ and CD8+ T cells) in each .FCS files which 
# was down-sampled from the data provided by Jalal Alshaweesh (Oct 2023) on [FlowRepository](http://flowrepository.org/id/FR-FCM-Z6UG). 
# The raw protein intensities and the corresponding metadata were saved in
# `data/ord.fcs.raw.mt_ffc.rda` and `data/ord.fcs.raw.meta.df.out_ffc.rda`.

## FFC - generating example data (random 50 per file)
{
  inpath="/path-to-FFCfcs/"
  subsampath="/path-to-subsamp/"
  
  ## Subsampling
  input_events_downsampling=50
  extra_args_read_FCS=NULL
  
  files <- list.files(inpath,full.names=TRUE,recursive=TRUE,pattern=".fcs")
  invisible(
    lapply(
      files,
      function(file){
        res <- do.call(read.FCS,c(list(filename=file),extra_args_read_FCS))
        w <- sort(sample(seq_len(nrow(res)),min(input_events_downsampling,nrow(res))))
        res <- res[w,]
        write.FCS(res,file.path(subsampath, basename(file)))
      }
    )
  )
}

