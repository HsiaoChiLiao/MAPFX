# Checks MapfxMPC
data(ord.fcs.raw.meta.df.out_mpc)
data(ord.fcs.raw.mt_mpc)
# create an Output directory in the current working directory for the argument 'Outpath' of the MapfxMPC function
dir.create(file.path(tempdir(), "MPC_impu_Output"))

# When `impute = TRUE`, randomly selecting 50% of the cells in each well for model training
set.seed(123) 
MapfxMPC_impu_obj <- MapfxMPC(
   runVignette = TRUE, #set FALSE if not running this example
   runVignette_meta = ord.fcs.raw.meta.df.out_mpc, #set NULL if not running this example
   runVignette_rawInten = ord.fcs.raw.mt_mpc, #set NULL if not running this example
   FCSpath = NULL, # users specify their own input path if not running this example
   Outpath = file.path(tempdir(), "MPC_impu_Output"),
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
   cores = 2L)

test_that("MapfxMPC returns an object", {
  obj <- MapfxMPC_impu_obj
  expect_true(validObject(obj))
})

test_that("MapfxMPC returns a list", {
  obj <- MapfxMPC_impu_obj
  expect_true(inherits(obj, "list"))
})

test_that("MapfxMPC returns a list with size 4", {
  obj <- MapfxMPC_impu_obj
  expect_length(obj, 4)
})


# When `impute = FALSE`, randomly selecting 50% of the cells in each well for model training
set.seed(123)
MapfxMPC_NOimpu_obj <- MapfxMPC(
  runVignette = TRUE, #set FALSE if not running this example
  runVignette_meta = ord.fcs.raw.meta.df.out_mpc, #set NULL if not running this example
  runVignette_rawInten = ord.fcs.raw.mt_mpc, #set NULL if not running this example
  FCSpath = NULL, # users specify their own input path if not running this example
  Outpath = file.path(tempdir(), "MPC_impu_Output"),
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
  plots.bkc.bkb = FALSE, plots.bkc.inf = FALSE, 
  plots.initM = FALSE,
  plots.rmWellEffect = FALSE,
  impute = FALSE,
  models.use = c("XGBoost"),
  extra_args_regression_params = list(list(nrounds = 1500, eta = 0.03)),
  prediction_events_downsampling = NULL,
  impu.training = FALSE,
  plots.imputation = FALSE,
  cluster.analysis.bkb = FALSE, plots.cluster.analysis.bkb = FALSE,
  cluster.analysis.all = FALSE, plots.cluster.analysis.all = FALSE,
  cores = 2L)

test_that("MapfxMPC returns an object", {
  obj <- MapfxMPC_NOimpu_obj
  expect_true(validObject(obj))
})

test_that("MapfxMPC returns a list", {
  obj <- MapfxMPC_NOimpu_obj
  expect_true(inherits(obj, "list"))
})

test_that("MapfxMPC returns a list with size 4", {
  obj <- MapfxMPC_NOimpu_obj
  expect_length(obj, 4)
})

test_that("MapfxMPC(...,impute = FALSE) does not return imputations", {
  obj <- MapfxMPC_NOimpu_obj
  expect_equal(obj$Imputation, NA)
})
