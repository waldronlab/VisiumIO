sample_dir <- system.file(
    "extdata", package = "VisiumIO", mustWork = TRUE
)

tvh <- TENxVisiumHD(
    spacerangerOut = sample_dir, bin_size = "002", images = "lowres"
)

expect_true(
    is(tvh, "TENxVisiumHD")
)
expect_true(
    validObject(tvh)
)

spe <- import(tvh)

expect_true(
    is(spe, "SpatialExperiment")
)

expect_identical(
    colnames(spe), rownames(colData(spe))
)

expect_identical(
    colnames(spe), colData(spe)[["barcode"]]
)

expect_identical(
    rownames(colData(spe)), colData(spe)[["barcode"]]
)

expect_identical(
    rownames(rowData(spe)), rowData(spe)[["ID"]]
)

expect_identical(
    SpatialExperiment::spatialCoordsNames(spe),
    c("pxl_col_in_fullres", "pxl_row_in_fullres")
)

# loadImage tests ---------------------------------------------------------

tvh_no_load <- TENxVisiumHD(
    spacerangerOut = sample_dir, bin_size = "002", images = "lowres",
    loadImage = FALSE
)
expect_false(
    tvh_no_load@loadImage
)

speh_no_load <- import(tvh_no_load)
imgh_no_load <- SpatialExperiment::imgData(speh_no_load)$data[[1L]]
expect_true(
    is(imgh_no_load, "StoredSpatialImage")
)

tvh_load <- TENxVisiumHD(
    spacerangerOut = sample_dir, bin_size = "002", images = "lowres",
    loadImage = TRUE
)
expect_true(
    tvh_load@loadImage
)

speh_load <- import(tvh_load)
imgh_load <- SpatialExperiment::imgData(speh_load)$data[[1]]
expect_true(
    is(imgh_load, "LoadedSpatialImage")
)

# test segmented_outputs input --------------------------------------------

library(sf)
sample_dir <- system.file(
    "extdata", "segmented_outputs", package = "VisiumIO", mustWork = TRUE
)

tvh <- TENxVisiumHD(
    segmented_outputs = sample_dir,
    format = "h5",
    processing = "filtered",
    images = "lowres"
)

expect_true(
    is(tvh, "TENxVisiumHD")
)

expect_true(
    validObject(tvh)
)

spe <- import(tvh)

expect_true(
    is(spe, "SpatialExperiment")
)

expect_identical(
    colnames(spe), rownames(colData(spe))
)

expect_identical(
    colnames(spe), rownames(colData(spe))
)


expect_identical(
    rownames(rowData(spe)), rowData(spe)[["ID"]]
)

expect_true(
    "cellseg" %in% names(metadata(spe))
)

expect_true(
    is.data.frame(metadata(spe)$cellseg)
)
expect_true(
    is(metadata(spe)$cellseg, "sf")
)
expect_identical(
    nrow(metadata(spe)$cellseg), ncol(spe)
)

