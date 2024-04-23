# Checks MapfxFFC
data(ord.fcs.raw.meta.df.out_ffc)
data(ord.fcs.raw.mt_ffc)

# create an Output directory in the current working directory for the argument 'Outpath' of the MapfxFFC function
dir.create(file.path(tempdir(), "FFCnorm_Output"))

MapfxFFC_obj <- MapfxFFC(
  runVignette = TRUE, #set FALSE if not running this example
  runVignette_meta = ord.fcs.raw.meta.df.out_ffc, #set NULL if not running this example
  runVignette_rawInten = ord.fcs.raw.mt_ffc, #set NULL if not running this example
  FCSpath = NULL, # users specify their own input path if not running this example
  Outpath = file.path(tempdir(), "FFCnorm_Output"),
  protein.v = c("CD3","CD4","CD8","CD45"),
  protein.upper.quantile = 0.9,
  protein.lower.quantile = 0.1,
  protein.min.quantile = 0.01,
  plots.bkc.protein = TRUE,
  plots.initM = TRUE,
  plots.rmBatchEffect = TRUE,
  cluster.analysis.protein = TRUE, plots.cluster.analysis.protein = TRUE)

test_that("MapfxFFC returns an object", {
  obj <- MapfxFFC_obj
  expect_true(validObject(obj))
})

test_that("MapfxFFC returns a list", {
  obj <- MapfxFFC_obj
  expect_true(inherits(obj, "list"))
})

test_that("MapfxFFC returns a list with size 2", {
  obj <- MapfxFFC_obj
  expect_length(obj, 2)
})
