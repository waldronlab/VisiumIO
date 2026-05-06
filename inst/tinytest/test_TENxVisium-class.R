outs_dir <- system.file(
    file.path("extdata", "10xVisium", "section1", "outs"),
    package = "VisiumIO"
)

tv <- TENxVisium(
    spacerangerOut = outs_dir, processing = "raw", images = "lowres"
)

expect_true(
    is(tv, "TENxVisium")
)
expect_true(
    validObject(tv)
)

expect_silent(
    TENxVisium(
        spacerangerOut = dirname(outs_dir),
        processing = "raw",
        images = "hires"
    )
)

tv <- TENxVisium(
    resources = TENxFileList(
        file.path(outs_dir, "raw_feature_bc_matrix")
    ),
    spatialResource = TENxSpatialList(
        resources = file.path(outs_dir, "spatial"),
        sample_id = "sample01",
        images = "lowres"
    ),
    processing = "raw",
    images = "lowres"
)

expect_true(
    validObject(tv)
)
expect_true(
    is(import(tv), "SpatialExperiment")
)

tv <- TENxVisium(
    spacerangerOut = TENxFileList(outs_dir),
    processing = "raw",
    images = "lowres"
)
expect_true(
    validObject(tv)
)
expect_true(
    is(import(tv), "SpatialExperiment")
)

expect_error(
    TENxVisium(
        spacerangerOut = TENxFileList(file.path(outs_dir, "spatial")),
        processing = "raw",
        images = "lowres",
        sample_id = "sample01"
    )
)

tv <- TENxVisium(
    resources = file.path(outs_dir, "raw_feature_bc_matrix"),
    spatialResource = file.path(outs_dir, "spatial"),
    images = "lowres",
    sample_id = "sample01"
)

expect_true(
    validObject(tv)
)
expect_true(
    is(import(tv), "SpatialExperiment")
)
expect_identical(
    dim(import(tv)),
    c(50L, 50L)
)

expect_error(
    VisiumIO:::.TENxVisium(
        resources = file.path(outs_dir, "raw_feature_bc_matrix"),
        spatialList =  TENxSpatialList(file.path(outs_dir, "spatial")),
        coordNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
        sampleId = "sample01"
    )
)

expect_error(
    VisiumIO:::.TENxVisium(
        resources = TENxFileList(
            file.path(outs_dir, "raw_feature_bc_matrix")
        ),
        spatialList =  file.path(outs_dir, "spatial"),
        coordNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
        sampleId = "sample01"
    )
)

# loadImage tests ---------------------------------------------------------

tv_no_load <- TENxVisium(
    spacerangerOut = outs_dir, processing = "raw", images = "lowres",
    loadImage = FALSE
)
expect_false(tv_no_load@loadImage)

spe_no_load <- import(tv_no_load)
img_no_load <- SpatialExperiment::imgData(spe_no_load)$data[[1L]]
expect_true(
    is(img_no_load, "StoredSpatialImage")
)

tv_load <- TENxVisium(
    spacerangerOut = outs_dir, processing = "raw", images = "lowres",
    loadImage = TRUE
)
expect_true(tv_load@loadImage)

spe_load <- import(tv_load)
img_load <- SpatialExperiment::imgData(spe_load)$data[[1L]]
expect_true(
    is(img_load, "LoadedSpatialImage")
)
