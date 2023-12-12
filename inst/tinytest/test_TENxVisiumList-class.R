sample_dirs <- list.dirs(
    system.file(
        file.path("extdata", "10xVisium"),
        package = "SpatialExperiment"
    ),
    recursive = FALSE, full.names = TRUE
)

tvl <- TENxVisiumList(
    sampleFolders = sample_dirs,
    sample_ids = c("sample01", "sample02"),
    processing = "raw",
    images = "lowres"
)

expect_true(
    is(tvl, "TENxVisiumList")
)
expect_true(
    validObject(tvl)
)

expect_true(
    is(import(tvl), "SpatialExperiment")
)
expect_true(
    validObject(import(tvl))
)

tvl <- TENxVisiumList(
    sampleFolders = sample_dirs,
    processing = "raw",
    images = "lowres"
)

expect_identical(
    vapply(tvl@VisiumList, function(x) x@sampleId, character(1L)),
    basename(sample_dirs)
)
