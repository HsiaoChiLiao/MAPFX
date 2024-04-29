# MAPFX: MAssively Parallel Flow cytometry Xplorer
### Author: Hsiao-Chi Liao
### Date: 30 April 2024
<br>

# Introduction
## Motivation
Massively-Parallel Cytometry (MPC) experiments allow cost-effective quantification of more than 200 surface proteins at single-cell resolution. [The Inflow protocol](https://www.science.org/doi/10.1126/sciadv.abg0505) is the pioneer of the pipeline for analysing MPC data, and the Bioconductor's `infinityFlow` package was developed for comprehensive analyses. However, the methods for background correction and removal of unwanted variation implemented in the package can be improved. We develop the `MAPFX` package as an alternative that has a more thoughtful strategy to clean up the raw protein intensities. Unique features of our package compared to the `infinityFlow` pipeline include performing background correction prior to imputation and removing unwanted variation from the data at the cell-level, while explicitly accounting for the potential association between biology and unwanted factors. We benchmarked our pipeline against the `infinityFlow` pipeline and demonstrated that our approach is better at preserving biological signals, removing unwanted variation, and imputing unmeasured infinity markers [MAPFX pipeline](https://doi.org/10.1101/2024.02.28.582452). Two user friendly functions `MapfxMPC` and `MapfxFFC` are included in the `MAPFX` package that were designed for data from either MPC or FFC experiments (see below sections for details). 

## Analysing Data from MPC Experiments
This package implemented an end-to-end toolbox for analysing raw data from MPC experiments. More details on the methodology can be found in [the MAPFX paper](https://doi.org/10.1101/2024.02.28.582452). The `MapfxMPC` function is designed for running through the whole pipeline. The pipeline starts by performing background correction on raw intensities to remove the noise from electronic baseline restoration and fluorescence compensation by adapting a normal-exponential convolution model. Unwanted technical variation, from sources such as well effects, is then removed using a log-normal model with plate, column, and row factors, after which infinity markers are imputed using the informative backbone markers as predictors with machine learning models. Cluster analysis and visualisation with UMAP two-dimensional representations can then be carried out if desired. Users can set `MapfxMPC(..., impute=FALSE)` if the imputation is not needed. 

## Analysing Data from the Fluorescence Flow Cytometry (FFC) Experiments
For the protein intensities from FFC experiments, the function `MapfxFFC` is used to carry out normalisation steps which include background correction and removal of unwanted variation, and the function can further perform cluster analysis and visualisation with UMAP two-dimensional representations if specified.

## Preparing Data for the Analysis - the Folder Diagram
```
# FCSpath
└───FCSpath
│   └───fcs
│       │   Plate1_A01.fcs
│       │   Plate1_A02.fcs
│       │   ...
│   └───meta
│       │   filename_meta.csv

# Outpath
└───Outpath
│   └───intermediary
│   └───downstream
│   └───graph

## Note: the sub-folders `intermediary`, `downstream`, and `graph` will 
## be generated automatically by MAPFX.
```

## Notes on Metadata
#### For MPC (the plate-based) Experiments
When set `file_meta = "auto"` for `MapfxMPC`, the file identifier keyword (GUID) of the FCS files MUST contain the following information and in the specified format:<br>
Plate information: Plate1, Plate2, …, Plate9<br>
Well information: A1, A2, …, A12, B1, …, H1, …, H12<br>

When set `file_meta = "usr"`, prepare `filename_meta.csv` in the following format and save the CSV file under `FCSpath/meta/`.<br>
An example:<br>
|  Filenam  |  Plate  |  Well  |  Column  |  Row  |  Well.lab  |
|:---------:|:-------:|:------:|:--------:|:-----:|:----------:|
|p1_a12.fcs |  Plate1 |   A12  |  Col.12  | Row.01|   P1_A12   |
|p2_d08.fcs |  Plate2 |   D08  |  Col.08  | Row.04|   P2_D08   |
|p3_g1.fcs  |  Plate3 |   G01  |  Col.01  | Row.07|   P3_G01   |
<br>
Note that the "Filenam" column refers to the GUID (file name) of each FCS file in the `FCSpath/fcs/`.
<br>

#### For FFC Experiments from Different Batches
Prepare `filename_meta.csv` in the following format and save the CSV file in `FCSpath/meta/`.<br>
An example:
|  Filenam  |  Batch  |
|:---------:|:-------:|
|090122.fcs |  Batch1 |
|070122.fcs |  Batch2 |
|010122.fcs |  Batch3 |
<br>


# Analysing Data with the MAPFX Package
## Installation
The MAPFX package can be installed using the code below.
``` r
if (!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools")
}
devtools::install_github("HsiaoChiLiao/MAPFX", build_vignettes=FALSE)

library(MAPFX)
```

Along with the MAPFX package, we also load the following packages required for running functions in MAPFX.
``` r
library(flowCore)
library(Biobase)
library(stringr)
library(uwot)
library(iCellR)
library(igraph)
library(ggplot2)
library(RColorBrewer)
library(Rfast)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(glmnetUtils)
library(e1071)
library(xgboost)
library(foreach)
library(doParallel)
library(parallel)
library(pbapply)
library(reshape2)
library(gtools)
library(utils)
library(stats)
library(cowplot)
```

## The Example Datasets in the MAPFX Package
### MPC
This dataset is a subset of the single-cell murine lung data at steady state downloaded from [FlowRepository](https://flowrepository.org/id/FR-FCM-Z2LP) provided by Etienne Becht (Nov 2020). The raw protein intensities and the corresponding metadata were saved in the objects `ord.fcs.raw.mt_mpc` and `ord.fcs.raw.meta.df.out_mpc` which were generated from 266 .FCS files from 266 wells with 50 cells in each file.

### FFC
This mice splenocytes dataset contains 50 cells (sorted CD4+ and CD8+ T cells) in each .FCS files which was down-sampled from the data provided by Jalal Alshaweesh (Oct 2023) on [FlowRepository](http://flowrepository.org/id/FR-FCM-Z6UG). The raw protein intensities and the corresponding metadata were saved in the objects `ord.fcs.raw.mt_ffc` and `ord.fcs.raw.meta.df.out_ffc`.

## `MapfxMPC(..., impute=TRUE)` - analysing data from MPC experiments
For users who would like to perform all of the following steps: background correction, removal of unwanted variation (well effects), **imputation**, and cluster analysis.

``` r
# import built-in data
data(ord.fcs.raw.meta.df.out_mpc)
data(ord.fcs.raw.mt_mpc)

# create an Output directory in the current working directory for the argument 'Outpath' of the MapfxMPC function
dir.create(file.path(tempdir(), "MPC_impu_Output"))

# usage
# when impute = TRUE, randomly selecting 50% of the cells in each well for model training
set.seed(123) 
MapfxMPC_impu_obj <- MapfxMPC(
    runVignette = TRUE, #set FALSE if not running this Vignette
    runVignette_meta = ord.fcs.raw.meta.df.out_mpc, #set NULL if not running this Vignette
    runVignette_rawInten = ord.fcs.raw.mt_mpc, #set NULL if not running this Vignette
    FCSpath = NULL, # users specify their own input path
    Outpath = file.path(tempdir(), "MPC_impu_Output"), # or users specify their own output path
    file_meta = "auto",
    bkb.v = c(
    "FSC-H", "FSC-W", "SSC-H", "SSC-W", "CD69-CD301b", "MHCII", 
    "CD4", "CD44", "CD8", "CD11c", "CD11b", "F480", 
    "Ly6C", "Lineage", "CD45a488", "CD24", "CD103"),
    yvar = "Legend", 
    control.wells = c(
    "P1_A01", "P2_A01", "P3_A01",
    "P3_F04", "P3_F05", "P3_F06", "P3_F07", "P3_F08", 
    "P3_F09", "P3_F10", "P3_F11", "P3_F12",
    "P3_G01", "P3_G02"),
    bkb.upper.quantile = 0.9, 
    bkb.lower.quantile = 0.1, 
    bkb.min.quantile = 0.01,
    inf.lower.quantile = 0.1, 
    inf.min.quantile = 0.01, 
    plots.bkc.bkb = TRUE, plots.bkc.inf = TRUE, 
    plots.initM = TRUE,
    plots.rmWellEffect = TRUE,
    impute = TRUE,
    models.use = c("XGBoost"),
    extra_args_regression_params = list(list(nrounds = 1500, eta = 0.03)),
    prediction_events_downsampling = NULL,
    impu.training = FALSE,
    plots.imputation = TRUE,
    cluster.analysis.bkb = TRUE, plots.cluster.analysis.bkb = TRUE,
    cluster.analysis.all = TRUE, plots.cluster.analysis.all = TRUE,
    cores = 4L)
    
# check the details
help(MapfxMPC, package = "MAPFX")
```

All the output will be stored in `file.path(tempdir(), "MPC_impu_Output")` (users can specify their own path: `/Outpath/`).

## `MapfxMPC(..., impute=FALSE)` - normalising data from MPC experiments
For users who would like to perform the following steps: background correction, removal of unwanted variation (well effects), and cluster analysis using backbones only.
``` r
# import built-in data
data(ord.fcs.raw.meta.df.out_mpc)
data(ord.fcs.raw.mt_mpc)

# create an Output directory in the current working directory for the argument 'Outpath' of the MapfxMPC function
dir.create(file.path(tempdir(), "MPC_NOimpu_Output"))

# usage
MapfxMPC_NOimpu_obj <- MapfxMPC(
    runVignette = TRUE, #set FALSE if not running this Vignette
    runVignette_meta = ord.fcs.raw.meta.df.out_mpc, #set NULL if not running this Vignette
    runVignette_rawInten = ord.fcs.raw.mt_mpc, #set NULL if not running this Vignette
    FCSpath = NULL, # users specify their own input path
    Outpath = file.path(tempdir(), "MPC_NOimpu_Output"), # or users specify their own output path
    file_meta="auto",
    bkb.v = c(
    "FSC-H", "FSC-W", "SSC-H", "SSC-W", "CD69-CD301b", "MHCII", 
    "CD4", "CD44", "CD8", "CD11c", "CD11b", "F480", 
    "Ly6C", "Lineage", "CD45a488", "CD24", "CD103"),
    yvar="Legend", 
    control.wells = c(
    "P1_A01", "P2_A01", "P3_A01",
    "P3_F04", "P3_F05", "P3_F06", "P3_F07", "P3_F08", 
    "P3_F09", "P3_F10", "P3_F11", "P3_F12",
    "P3_G01", "P3_G02"),
    bkb.upper.quantile = 0.9, 
    bkb.lower.quantile = 0.1, 
    bkb.min.quantile = 0.01,
    inf.lower.quantile = 0.1, 
    inf.min.quantile = 0.01, 
    plots.bkc.bkb = TRUE, plots.bkc.inf = TRUE, 
    plots.initM = TRUE,
    plots.rmWellEffect = TRUE,
    impute = FALSE,
    cluster.analysis.bkb = TRUE, plots.cluster.analysis.bkb = TRUE,
    cores = 4L)

# check the details
help(MapfxMPC, package = "MAPFX")
```

All the output will be stored in `file.path(tempdir(), "MPC_NOimpu_Output")` (users can specify their own path: `/Outpath/`).

## `MapfxFFC` - normalising data from FFC experiments
For users who would like to perform the following steps: background correction, removal of unwanted variation (batch effects), and cluster analysis.
``` r
# import built-in data
data(ord.fcs.raw.meta.df.out_ffc)
data(ord.fcs.raw.mt_ffc)

# create an Output directory in the current working directory for the argument 'Outpath' of the MapfxFFC function
dir.create(file.path(tempdir(), "FFCnorm_Output"))

MapfxFFC_obj <- MapfxFFC(
    runVignette = TRUE, #set FALSE if not running this Vignette
    runVignette_meta = ord.fcs.raw.meta.df.out_ffc, #set NULL if not running this Vignette
    runVignette_rawInten = ord.fcs.raw.mt_ffc, #set NULL if not running this Vignette
    FCSpath = NULL, # users specify their own input path
    Outpath = file.path(tempdir(), "FFCnorm_Output"), # or users specify their own output path
    protein.v = c("CD3","CD4","CD8","CD45"),
    protein.upper.quantile = 0.9, 
    protein.lower.quantile = 0.1, 
    protein.min.quantile = 0.01,
    plots.bkc.protein = TRUE,
    plots.initM = TRUE,
    plots.rmBatchEffect = TRUE,
    cluster.analysis.protein = TRUE, plots.cluster.analysis.protein = TRUE)
    
# check the details
help(MapfxFFC, package = "MAPFX")
```

All the output will be stored in `file.path(tempdir(), "FFCnorm_Output")` (users can specify their own path: `/Outpath/`).

## Description of the output
Three folders will be automatically generated in `/Outpath/`.<br>
1. `intermediary`:<br>
Intermediary results will be saved in the `.rds` and `.RData` formats and will be stored here.<br>
2. `downstream`:<br>
Final results will be saved in the `.rds` format and will be stored here. The results include normalised backbone measurements (on both linear and log scale: `bkc.adj.bkb_linearScale_mt.rds` and `bkc.adj.bkb_logScale_mt.rds`), the completed dataset with imputed infinity (exploratory) markers (`predictions.Rds`), UMAP coordinates derived from both normalised backbones (`ClusterAnalysis_umap_#bkb.rds`) and the completed dataset (`ClusterAnalysis_ImpuMtd_umap_#bkb.#impuPE.rds`), and metadata (`fcs_metadata_df.rds`) for cells including cluster labels derived from both normalised backbones and the completed data matrix.<br>
3. `graph`:<br>
Figures will be stored here, including **scatter plots** for comparing background corrected and raw protein intensities for each protein marker, **heatmaps** for presenting the biological and unwanted effects in the data before and after removal of unwanted variation with *mapfx.norm*, **boxplots** (for imputations from multiple models) and **a boxplot and a histogram** (for imputations from a single model) of R-sq values for visualising the accuracy of imputed infinity (exploratory) markers, and **UMAP plots** for showing the cluster structure.<br>

## The "Larger" Example Datasets in the MapfxData Package
**Soon will be available from the Bioconductor Experiment Data Package - MapfxData.**

# THANK YOU

Thanks for carrying out analyses with our `MAPFX` package. Please feel free to raise any issues and/or questions on our GitHub page and/or send them to Hsiao-Chi [hsiaochi.liao@student.unimelb.edu.au](mailto:hsiaochi.liao@student.unimelb.edu.au).

## Information about the R session when this vignette was built

``` r
sessionInfo()
## R Under development (unstable) (2024-03-18 r86148)
## Platform: aarch64-apple-darwin20
## Running under: macOS Sonoma 14.4
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
## 
## Random number generation:
##  RNG:     L'Ecuyer-CMRG 
##  Normal:  Inversion 
##  Sample:  Rejection 
##  
## locale:
## [1] C/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8
## 
## time zone: Australia/Melbourne
## tzcode source: internal
## 
## attached base packages:
## [1] parallel  grid      stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] cowplot_1.1.3         gtools_3.9.5          reshape2_1.4.4       
##  [4] pbapply_1.7-2         doParallel_1.0.17     iterators_1.0.14     
##  [7] foreach_1.5.2         xgboost_1.7.7.1       e1071_1.7-14         
## [10] glmnetUtils_1.1.9     circlize_0.4.16       ComplexHeatmap_2.19.0
## [13] ggrepel_0.9.5         Rfast_2.1.0           RcppParallel_5.1.7   
## [16] RcppZiggurat_0.1.6    Rcpp_1.0.12           RColorBrewer_1.1-3   
## [19] igraph_2.0.3          iCellR_1.6.7          plotly_4.10.4        
## [22] ggplot2_3.5.0         uwot_0.1.16           Matrix_1.6-5         
## [25] stringr_1.5.1         Biobase_2.63.0        BiocGenerics_0.49.1  
## [28] flowCore_2.15.2       MAPFX_0.99.0          knitr_1.45           
## [31] BiocStyle_2.31.0     
## 
## loaded via a namespace (and not attached):
##   [1] ggdendro_0.2.0       rstudioapi_0.16.0    jsonlite_1.8.8      
##   [4] shape_1.4.6.1        magrittr_2.0.3       farver_2.1.1        
##   [7] rmarkdown_2.26       GlobalOptions_0.1.2  vctrs_0.6.5         
##  [10] base64enc_0.1-3      rstatix_0.7.2        htmltools_0.5.8     
##  [13] progress_1.2.3       broom_1.0.5          Formula_1.2-5       
##  [16] sass_0.4.9           bslib_0.7.0          htmlwidgets_1.6.4   
##  [19] plyr_1.8.9           cachem_1.0.8         mime_0.12           
##  [22] lifecycle_1.0.4      pkgconfig_2.0.3      R6_2.5.1            
##  [25] fastmap_1.1.1        shiny_1.8.1          clue_0.3-65         
##  [28] digest_0.6.35        reshape_0.8.9        colorspace_2.1-0    
##  [31] S4Vectors_0.41.3     irlba_2.3.5.1        Hmisc_5.1-2         
##  [34] ggpubr_0.6.0         labeling_0.4.3       cytolib_2.15.2      
##  [37] fansi_1.0.6          httr_1.4.7           abind_1.4-5         
##  [40] compiler_4.4.0       proxy_0.4-27         withr_3.0.0         
##  [43] bit64_4.0.5          htmlTable_2.4.2      backports_1.4.1     
##  [46] carData_3.0-5        ggsignif_0.6.4       MASS_7.3-60.2       
##  [49] rjson_0.2.21         scatterplot3d_0.3-44 tools_4.4.0         
##  [52] foreign_0.8-86       ape_5.7-1            httpuv_1.6.15       
##  [55] nnet_7.3-19          glue_1.7.0           nlme_3.1-164        
##  [58] promises_1.2.1       checkmate_2.3.1      Rtsne_0.17          
##  [61] cluster_2.1.6        generics_0.1.3       hdf5r_1.3.10        
##  [64] gtable_0.3.4         class_7.3-22         tidyr_1.3.1         
##  [67] data.table_1.15.2    hms_1.1.3            car_3.1-2           
##  [70] utf8_1.2.4           RcppAnnoy_0.0.22     RANN_2.6.1          
##  [73] pillar_1.9.0         later_1.3.2          splines_4.4.0       
##  [76] dplyr_1.1.4          lattice_0.22-6       FNN_1.1.4           
##  [79] survival_3.5-8       bit_4.0.5            RProtoBufLib_2.15.0 
##  [82] tidyselect_1.2.1     gridExtra_2.3        bookdown_0.38       
##  [85] IRanges_2.37.1       stats4_4.4.0         xfun_0.43           
##  [88] matrixStats_1.2.0    pheatmap_1.0.12      stringi_1.8.3       
##  [91] lazyeval_0.2.2       yaml_2.3.8           evaluate_0.23       
##  [94] codetools_0.2-19     NbClust_3.0.1        tibble_3.2.1        
##  [97] BiocManager_1.30.22  cli_3.6.2            rpart_4.1.23        
## [100] xtable_1.8-4         munsell_0.5.0        jquerylib_0.1.4     
## [103] png_0.1-8            prettyunits_1.2.0    glmnet_4.1-8        
## [106] viridisLite_0.4.2    scales_1.3.0         purrr_1.0.2         
## [109] crayon_1.5.2         GetoptLong_1.0.5     rlang_1.1.3
```
