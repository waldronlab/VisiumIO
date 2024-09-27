#' @include TENxVisiumList-class.R

#' @docType class
#'
#' @title A class to represent and import multiple Visium HD samples
#'
#' @description This class contains a `SimpleList` of [TENxVisiumHD] objects
#'   each corresponding to one sample. The provided `spacerangerOut` folder
#'   should contain a `binned_outputs` folder where multiple `bin_size`
#'   subfolders are present, e.g., `square_002um`.
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
.TENxVisiumHD <- setClass(Class = "TENxVisiumHD", contains = "TENxVisium")

.getSpatialPath <- function(path, bin_size) {
    outputs <- file.path(path, "binned_outputs")
    stopifnot(
        "The 'binned_outputs' directory was not found." = dir.exists(outputs)
    )
    squaref <- paste0("square_", bin_size, "um")
    spatf <- file.path(outputs, squaref, "spatial")
    stopifnot(
        "The 'spatial' directory was not found." = all(dir.exists(spatf))
    )
    spatf
}

.filter_sort_mtx_files <- function(namesvec) {
    files <- .FEATURE_BC_MATRIX_FILES
    names(files) <- files
    res <- lapply(files, function(file) {
        namesvec[startsWith(namesvec, file)]
    })
    unlist(res)
}

.exclude_mtx_files <- function(filelist) {
    files <- .FEATURE_BC_MATRIX_FILES
    names(files) <- files
    res <- lapply(files, function(file) {
        startsWith(names(filelist), file)
    })
    filelist[!Reduce(`|`, res)]
}

.exclude_h5_files <- function(filelist) {
    filelist[tools::file_ext(names(filelist)) != "h5"]
}

.check_filter_mtx <- function(filelist) {
    afiles <- .filter_sort_mtx_files(names(filelist))
    if (!identical(names(afiles), .FEATURE_BC_MATRIX_FILES))
        stop(
            "'TENxFileList' does not contain the expected files:\n  ",
            .FEATURE_BC_MATRIX_FILES_PRINT
        )
    filelist[afiles]
}

.find_convert_resources_hd <-
    function(path, processing, format, bin_size, ...)
{
    if (!is(path, "TENxFileList")) {
        squaref <- .getSpatialPath(path, bin_size) |> dirname()
        path <-  .find_file_or_dir(squaref, processing, format)
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
    } else {
        path <- .check_filter_mtx(path)
    }
    path
}

.find_convert_spatial_hd <- function(path, bin_size, ...) {
    if (!is(path, "TENxFileList")) {
        path <- .getSpatialPath(path, bin_size)
    } else {
        path <- .exclude_mtx_files(path)
        path <- .exclude_h5_files(path)
    }
    TENxSpatialList(path, ...)
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
#' @examples
#'
#' vdir <- system.file(
#'     "extdata", package = "VisiumIO", mustWork = TRUE
#' )
#'
#' ## with spacerangerOut folder
#' TENxVisiumHD(spacerangerOut = vdir, bin_size = "002", images = "lowres")
#'
#' TENxVisiumHD(spacerangerOut = vdir, bin_size = "002", images = "lowres") |>
#'     import()
#'
#' ## indicate h5 format
#' TENxVisiumHD(
#'     spacerangerOut = vdir, bin_size = "002",
#'     images = "lowres", format = "h5"
#' )
#'
#' TENxVisiumHD(
#'     spacerangerOut = vdir, bin_size = "002",
#'     images = "lowres", format = "h5"
#' ) |>
#'     import()
#'
#' ## use resources and spatialResource arguments as file paths
#' TENxVisiumHD(
#'     resources = file.path(
#'         vdir, "binned_outputs", "square_002um",
#'         "filtered_feature_bc_matrix.h5"
#'     ),
#'     spatialResource = file.path(
#'         vdir, "binned_outputs", "square_002um",
#'         "spatial"
#'     ),
#'     bin_size = "002", processing = "filtered",
#'     images = "lowres", format = "h5"
#' ) |>
#'     import()
#'
#' ## provide the spatialResource argument as a TENxFileList
#' TENxVisiumHD(
#'     resources = file.path(
#'         vdir, "binned_outputs", "square_002um",
#'         "filtered_feature_bc_matrix.h5"
#'     ),
#'     spatialResource = TENxFileList(
#'         "~/gh/VisiumIO/inst/extdata/binned_outputs/square_002um/spatial/"
#'     ),
#'     bin_size = "002", images = "lowres", format = "h5"
#' ) |>
#'     import()
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
    format <- match.arg(format)

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
                is(resources, "TENxFileList_OR_TENxH5"),
            (isScalarCharacter(spatialResource) &&
                file.exists(spatialResource)) ||
                    is(spatialResource, "TENxFileList")
        )
        if (
            !is(resources, "TENxFileList_OR_TENxH5") &&
            identical(tools::file_ext(resources), "h5")
        )
            resources <- TENxH5(resources, ranges = NA_character_)
        else if (is.character(resources))
            resources <- TENxFileList(resources, ...)
        if (!is(spatialResource, "TENxSpatialList"))
            spatialResource <- TENxSpatialList(
                resources = spatialResource, sample_id = sample_id,
                images = images, jsonFile = jsonFile,
                tissuePattern = tissuePattern
            )
    }

    txv <- TENxVisium(
        resources = resources,
        spatialResource = spatialResource,
        sampleId = sample_id,
        processing = processing,
        format = format,
        images = images,
        jsonFile = jsonFile,
        tissuePattern = tissuePattern,
        spatialCoordsNames = spatialCoordsNames,
        ...
    )

    .TENxVisiumHD(txv)
}

# import TENxVisiumHD method ----------------------------------------------

#' @describeIn TENxVisiumHD-class Import Visium HD data from multiple bin sizes
#'
#' @inheritParams TENxVisiumList
#'
#' @exportMethod import
setMethod("import", "TENxVisiumHD", function(con, format, text, ...) {
    methods::callNextMethod()
})
