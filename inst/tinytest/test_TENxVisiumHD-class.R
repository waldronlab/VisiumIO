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

