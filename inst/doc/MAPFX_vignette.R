## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- required packages, eval = FALSE-----------------------------------------
#  library(flowCore)
#  library(Biobase)
#  library(stringr)
#  library(uwot)
#  library(Rphenograph)
#  library(igraph)
#  library(ggplot2)
#  library(RColorBrewer)
#  library(Rfast)
#  library(ggrepel)
#  library(ComplexHeatmap)
#  library(circlize)
#  library(glmnetUtils)
#  library(e1071)
#  library(xgboost)
#  library(foreach)
#  library(doParallel)
#  library(parallel)
#  library(pbapply)
#  library(reshape)
#  library(gtools)

## ---- install the package, eval = FALSE---------------------------------------
#  if (!requireNamespace("devtools", quietly = TRUE)){
#      install.packages("devtools")
#  }
#  devtools::install_github("HsiaoChiLiao/MAPFX", build_vignettes=FALSE)

## ---- load_package, eval = FALSE----------------------------------------------
#  library(MAPFX)

## ---- running the pipeline with runMAPFX, eval = FALSE------------------------
#  # running with build-in FCS files:
#  runMAPFX(
#  FCSpath="/PathToFCSfiles/",
#  Outpath="/OutputPath/Output/",  #the folder for saving output should be named 'Output'
#  file_meta="auto",
#  batch_key="plate", ##
#  bkb.v = c("FSC-H", "FSC-W", "SSC-H", "SSC-W", "CD69-CD301b", "MHCII", "CD4", "CD44", "CD8", "CD11c", "CD11b",
#  "F480", "Ly6C", "Lineage", "CD45a488", "CD24", "CD103"),
#  chans = c("FSC-H", "FSC-W", "SSC-H", "SSC-W", "CD69-CD301b", "MHCII", "CD4", "CD44", "CD8", "CD11c", "CD11b",
#  "F480", "Ly6C", "Lineage", "CD45a488", "CD24", "CD103"),
#  yvar="Legend", control.wells = c("P1_A01", "P2_A02", "P3_G02"),
#  bkb.upper.quntile=0.9, bkb.lower.quntile=0.1, bkb.min.quntile=0.01,
#  pe.mean.sd=3, pe.lower.quntile=0.1, pe.min.quntile=0.01,
#  trans.dat="cent.lgc",
#  models.use = c("XGBoost"),
#  extra_args_regression_params = list(list(nrounds = 1500, eta = 0.03)),
#  cores=4L)
#  
#  # check the usage
#  help(runMAPFX, package = "MAPFX")

## ---- running the pipeline with runMAPFX_ImpuMICE, eval = FALSE---------------
#  # running with build-in FCS files:
#  runMAPFX_ImpuMICE(
#  FCSpath="/PathToFCSfiles/",
#  Outpath="/OutputPath/Output/",  #the folder for saving output should be named 'Output'
#  file_meta="auto",
#  batch_key="plate", ##
#  bkb.v = c("FSC-H", "FSC-W", "SSC-H", "SSC-W", "CD69-CD301b", "MHCII", "CD4", "CD44", "CD8", "CD11c", "CD11b",
#  "F480", "Ly6C", "Lineage", "CD45a488", "CD24", "CD103"),
#  chans = c("FSC-H", "FSC-W", "SSC-H", "SSC-W", "CD69-CD301b", "MHCII", "CD4", "CD44", "CD8", "CD11c", "CD11b",
#  "F480", "Ly6C", "Lineage", "CD45a488", "CD24", "CD103"),
#  yvar="Legend", control.wells = c("P1_A01", "P2_A02", "P3_G02"),
#  bkb.upper.quntile=0.9, bkb.lower.quntile=0.1, bkb.min.quntile=0.01,
#  pe.mean.sd=3, pe.lower.quntile=0.1, pe.min.quntile=0.01,
#  trans.dat="cent.lgc",
#  models.use = c("XGBoost"),
#  extra_args_regression_params = list(list(nrounds = 1500, eta = 0.03)),
#  cores=4L)
#  
#  # check the usage
#  help(runMAPFX, package = "MAPFX")

## ---- running the pipeline with runMAPFX_bkbOnly, eval = FALSE----------------
#  # running with build-in FCS files:
#  runMAPFX_bkbOnly(
#  FCSpath="/PathToFCSfiles/",
#  Outpath="/OutputPath/Output/",  #the folder for saving output should be named 'Output'
#  batch_key="plate", ##
#  bkb.v = c("FSC-H", "FSC-W", "SSC-H", "SSC-W", "CD69-CD301b", "MHCII", "CD4", "CD44", "CD8", "CD11c", "CD11b",
#  "F480", "Ly6C", "Lineage", "CD45a488", "CD24", "CD103"),
#  bkb.upper.quntile=0.9, bkb.lower.quntile=0.1, bkb.min.quntile=0.01,
#  trans.dat="cent.lgc")
#  
#  # check the usage
#  help(runMAPFX_bkbOnly, package = "MAPFX")

## -----------------------------------------------------------------------------
sessionInfo()

