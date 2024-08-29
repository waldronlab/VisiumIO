#' @include TENxVisiumList-class.R

#' @docType class
#'
#' @title A class to represent and import multiple Visium HD samples
#'
#' @description This class contains a `SimpleList` of [TENxVisiumHD] objects
#'   each corresponding to one sample.
#'
#' @details Typically, the user will provide a path to a directory containing
#'  the output of the `spaceranger count` command. The `spaceranger count`
#'  command outputs a folder containing the "raw" or "filtered"
#'  `()_feature_bc_matrix`.
#'
#' @inheritParams TENxVisiumList-class
#'
#' @return A [SpatialExperiment] object
#'
#' @exportClass TENxVisiumHD
.TENxVisiumHD <- setClass(Class = "TENxVisiumHD", contains = "TENxVisiumList")

.getSpatialPath <- function(path, bin_size) {
    outputs <- list.dirs(path, recursive = FALSE, full.names = TRUE)
    stopifnot(
        "The 'binned_outputs' directory was not found." =
            endsWith(outputs, "binned_outputs")
    )
    squaref <- paste0("square_", bin_size, "um")
    spatf <- file.path(outputs, squaref, "spatial")
    stopifnot(
        "The 'spatial' directory was not found." = all(dir.exists(spatf))
    )
    spatf
}

.find_convert_resources_hd <-
    function(path, processing, format, bin_size, ...)
{
    if (!is(path, "TENxFileList")) {
        squaref <- .getSpatialPath(path, bin_size) |> dirname()
        path <- vapply(
            squaref, .find_file_or_dir, character(1L), processing, format
        )
        fdirname <- paste0(processing, "_feature_bc_matrix")
        if (
            !all(dir.exists(file.path(squaref, fdirname))) ||
            !all(dir.exists(file.path(squaref, "spatial")))
        )
            stop(
                "The 'spatial' or '", fdirname, "' directory was not found.",
                "\n  Verify 'spacerangerOut' and 'processing' inputs.",
                call. = FALSE
            )
        path <- file.path(squaref, fdirname)
    } else {
        if (!all(.FEATURE_BC_MATRIX_FILES %in% names(path)))
            stop(
                "'TENxFileList' does not contain the expected files:\n  ",
                .FEATURE_BC_MATRIX_FILES_PRINT
            )
        path <- path[.FEATURE_BC_MATRIX_FILES]
    }
    FileFUN <- if (identical(unique(tools::file_ext(path)), "h5"))
            function(x) {TENxH5(x, ranges = NA_character_) }
        else if (identical(unique(tools::file_ext(path)), ""))
            TENxFileList

    mapply(
        FileFUN,
        path,
        MoreArgs = list(...),
        SIMPLIFY = FALSE
    )
}

.find_convert_spatial_hd <- function(path, bin_size, ...) {
    if (!is(path, "TENxFileList")) {
        path <- .getSpatialPath(path, bin_size)
    } else {
        path <- path[!names(path) %in% .FEATURE_BC_MATRIX_FILES]
    }
    mapply(
        TENxSpatialList,
        path,
        MoreArgs = list(...),
        SIMPLIFY = FALSE
    )
}

#' @rdname TENxVisiumHD-class
#'
#' @inheritParams TENxVisium
#' @inheritParams TENxVisiumList
#'
#' @param bin_size `character(1)` The bin size of the images to import. The
#'   default is `008`. It corresponds to the directory name `square_000um` where
#'   `000` is the bin value.
#'
#' @export
TENxVisiumHD <- function(
    resources,
    spatialResource,
    spacerangerOut,
    sample_id = "sample01",
    processing = c("filtered", "raw"),
    format = c("mtx", "h5"),
    images = c("lowres", "hires", "detected", "aligned_fiducials"),
    bin_size = c("008", "016", "002"),
    jsonFile = .SCALE_JSON_FILE,
    tissuePattern = "tissue_positions\\.parquet",
    spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    ...
) {
    images <- match.arg(images, several.ok = TRUE)
    processing <- match.arg(processing)
    bin_size <- match.arg(bin_size)

    if (!missing(spacerangerOut)) {
        if (isScalarCharacter(spacerangerOut))
            stopifnot(
                dir.exists(spacerangerOut)
            )
        resources <- .find_convert_resources_hd(
            spacerangerOut, processing, format, bin_size, ...
        )
        spatialResource <- .find_convert_spatial_hd(
            path = spacerangerOut, bin_size = bin_size, sample_id = sample_id,
            images = images, jsonFile = jsonFile, tissuePattern = tissuePattern
        )
    } else {
        stopifnot(
            (isScalarCharacter(resources) && file.exists(resources)) ||
                is(resources, "TENxFileList"),
            (isScalarCharacter(spatialResource) &&
                file.exists(spatialResource)) ||
                    is(spatialResource, "TENxSpatialList")
        )
        if (!is(resources, "TENxFileList"))
            resources <- TENxFileList(resources, ...)
        if (!is(spatialResource, "TENxSpatialList"))
            spatialResource <- TENxSpatialList(
                resources = spatialResource, sample_id = sample_id,
                images = images, jsonFile = jsonFile,
                tissuePattern = tissuePattern
            )
    }

    txvList <- mapply(
        .TENxVisium,
        resources = resources,
        spatialList = spatialResource,
        MoreArgs = list(
            coordNames = spatialCoordsNames,
            sampleId = sample_id
        )
    )

    .TENxVisiumHD(
        VisiumList = SimpleList(txvList)
    )
}

#' @describeIn TENxVisiumHD-class Import Visium HD data from multiple bin sizes
#'
#' @inheritParams TENxVisiumList
#'
#' @exportMethod import
setMethod("import", "TENxVisiumHD", function(con, format, text, ...) {
    SElist <- lapply(con@VisiumList, import)
    do.call(cbind, SElist)
})
