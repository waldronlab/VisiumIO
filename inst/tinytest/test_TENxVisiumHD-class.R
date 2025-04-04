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
