sample_dir <- system.file(
    file.path("extdata", "10xVisium", "section1"),
    package = "SpatialExperiment"
)

tv <- TENxVisium(
    spacerangerOut = sample_dir, processing = "raw", images = "lowres"
)

expect_true(
    is(tv, "TENxVisium")
)
expect_true(
    validObject(tv)
)

expect_error(
    TENxVisium(
        spacerangerOut = dirname(sample_dir),
        processing = "raw",
        images = "hires"
    )
)

tv <- TENxVisium(
    resources = TENxFileList(
        file.path(sample_dir, "outs/raw_feature_bc_matrix")
    ),
    spatialResource = TENxSpatialList(
        resources = file.path(sample_dir, "outs/spatial"),
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
    spacerangerOut = TENxFileList(sample_dir),
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
        spacerangerOut = TENxFileList(file.path(sample_dir, "outs/spatial")),
        processing = "raw",
        images = "lowres",
        sample_id = "sample01"
    )
)

tv <- TENxVisium(
    resources = file.path(sample_dir, "outs/raw_feature_bc_matrix"),
    spatialResource = file.path(sample_dir, "outs/spatial"),
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
        resources = file.path(sample_dir, "outs/raw_feature_bc_matrix"),
        spatialList =  TENxSpatialList(file.path(sample_dir, "outs/spatial")),
        coordNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
        sampleId = "sample01"
    )
)

expect_error(
    VisiumIO:::.TENxVisium(
        resources = TENxFileList(
            file.path(sample_dir, "outs/raw_feature_bc_matrix")
        ),
        spatialList =  file.path(sample_dir, "outs/spatial"),
        coordNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
        sampleId = "sample01"
    )
)
