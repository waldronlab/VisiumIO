spatial_dir <- system.file(
    file.path("extdata", "10xVisium", "section1", "outs", "spatial"),
    package = "SpatialExperiment"
)

txsl <- TENxSpatialList(resources = spatial_dir)

expect_true(
    validObject(txsl)
)
expect_true(
    is(txsl, "TENxSpatialList")
)

# When none are specified all images are listed
expect_identical(
    txsl@images,
    c("lowres", "hires", "detected", "aligned")
)
expect_identical(
    txsl@scaleJSON,
    grep("scalefactors", dir(spatial_dir), value = TRUE)
)
expect_identical(
    txsl@tissuePos,
    grep("tissue_positions", dir(spatial_dir), value = TRUE)
)
expect_identical(
    txsl@sampleId,
    "sample01"
)

txsl <- TENxSpatialList(resources = spatial_dir, sample_id = "sample02")
expect_identical(
    txsl@sampleId,
    "sample02"
)

txsl <- TENxSpatialList(resources = spatial_dir, images = "lowres")
expect_identical(
    txsl@images,
    "lowres"
)

txsl <- TENxSpatialList(resources = spatial_dir, images = "aligned")
expect_error(
    import(txsl),
    "The 'aligned' image was not found.*"
)

expect_error(
    TENxSpatialList(
        resources = spatial_dir, jsonFile = "scalefactors.json"
    )
)

txsl <- TENxSpatialList(resources = spatial_dir, images = "lowres")
spatlist <- import(txsl)
expect_true(
    is.list(
        spatlist
    )
)
expect_identical(
    names(spatlist),
    c("imgData", "colData")
)
expect_true(
    is(spatlist$imgData, "DataFrame")
)
expect_true(
    is(spatlist$colData, "DataFrame")
)

expect_error(
    TENxSpatialList(
        resources = spatial_dir, images = "lowres", compressed = TRUE
    )
)
