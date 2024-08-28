#' @docType class
#'
#' @title A class to represent and import a single Visium Sample
#'
#' @description This class is a composed class of [TENxFileList] which can
#'   contain a list of [TENxFile] objects and a [TENxSpatialList] object. It is
#'   meant to handle a single Visium sample from 10X Genomics.
#'
#' @details Typically, the user will not create an object of this class directly
#'   but rather use [TENxVisiumList] constructor function for multiple samples.
#'   Note that the `images`, `jsonFile`, `tissuePattern`, and
#'   `spatialCoordsNames` arguments are only considered when the
#'   `spacerangerOut` argument or both the `resources` and `spatialResource`
#'   arguments are paths to files.
#'
#' @slot resources A [TENxFileList] object containing the Visium data.
#'
#' @slot spatialList A [TENxSpatialList] object containing the spatial
#'
#' @slot coordNames `character()` A vector specifying the names
#'   of the columns in the spatial data containing the spatial coordinates.
#'
#' @slot sampleId `character(1)` A scalar specifying the sample identifier.
#'
#' @return A [SpatialExperiment] object
#'
#' @seealso <https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/overview>
#'
#' @exportClass TENxVisium
.TENxVisium <- setClass(
    Class = "TENxVisium",
    slots = c(
        resources = "TENxFileList",
        spatialList = "TENxSpatialList",
        coordNames = "character",
        sampleId = "character"
    )
)

.FEATURE_BC_MATRIX_FILES <- c("barcodes.tsv", "features.tsv", "matrix.mtx")

.FEATURE_BC_MATRIX_FILES_PRINT <- paste0(
    sQuote(.FEATURE_BC_MATRIX_FILES), collapse = ", "
)

#' @importFrom TENxIO TENxFileList
.find_convert_resources <- function(path, processing, ...) {
    if (!is(path, "TENxFileList")) {
        odir <- list.dirs(path, recursive = FALSE, full.names = TRUE)
        fdirname <- paste0(processing, "_feature_bc_matrix")
        featdir <- file.path(odir, fdirname)
        if (
            identical(length(odir), 1L) && endsWith(odir, "outs") &&
            identical(length(featdir), 1L) && dir.exists(featdir)
        )
            path <- featdir
        else
            stop(
                "The 'outs' or '", fdirname, "' directory was not found.",
                "\n  Verify 'spacerangerOut' and 'processing' inputs.",
                call. = FALSE
            )
    } else {
        if (!all(.FEATURE_BC_MATRIX_FILES %in% names(path)))
            stop(
                "'TENxFileList' does not contain the expected files:\n  ",
                .FEATURE_BC_MATRIX_FILES_PRINT
            )
        path <- path[.FEATURE_BC_MATRIX_FILES]
    }
    TENxFileList(path, ...)
}

.find_convert_spatial <- function(path, ...) {
    if (!is(path, "TENxFileList")) {
        odir <- list.dirs(path, recursive = FALSE, full.names = TRUE)
        if (endsWith(odir, "outs"))
            path <- file.path(odir, "spatial")
        else
            stop("The 'outs' directory was not found")
    } else {
        path <- path[!names(path) %in% .FEATURE_BC_MATRIX_FILES]
    }
    TENxSpatialList(path, ...)
}

#' @rdname TENxVisium-class
#'
#' @param resources A [TENxFileList] object or a file path to the tarball
#'   containing the matrix / assay data resources.
#'
#' @param spatialResource A [TENxSpatialList] object or a file path to the
#'   tarball containing the spatial data.
#'
#' @param spacerangerOut `character(1)` A single string specifying the path to
#'   the sample directory of `spaceranger count`. The directory must contain the
#'   `(processing)_feature_bc_matrix` and `spatial` sub directories in addition
#'   to the `outs` folder.
#'
#' @param sample_id `character(1)` A single string specifying the sample ID.
#'
#' @param processing `character(1)` A single string indicating the processing
#'   folder available e.g., "filtered_feature_barcode_matrix" in the
#'   `spacerangerOut` folder. It can be either "filtered" or "raw" (default
#'   "filtered"). Only used when `spacerangerOut` is specified.
#'
#' @param images `character()` A vector specifying the images to be imported;
#'   can be one or multiple of "lowres", "hires", "detected", "aligned".
#'
#' @param jsonFile `character(1)` A single string specifying the name of the
#'  JSON file containing the scale factors.
#'
#' @param tissuePattern `character(1)` A single string specifying the pattern
#'   to match the tissue positions file.
#'
#' @param spatialCoordsNames `character()` A vector of strings specifying the
#'  names of the columns in the spatial data containing the spatial coordinates.
#'
#' @param ... In the constructor, additional arguments passed to
#'   [TENxFileList()]; otherwise, not used.
#'
#' @importFrom BiocBaseUtils isScalarCharacter
#'
#' @examples
#'
#' sample_dir <- system.file(
#'     file.path("extdata", "10xVisium", "section1"),
#'     package = "SpatialExperiment"
#' )
#'
#' tv <- TENxVisium(
#'     spacerangerOut = sample_dir, processing = "raw", images = "lowres"
#' )
#'
#' import(tv)
#'
#' @export
TENxVisium <- function(
    resources,
    spatialResource,
    spacerangerOut,
    sample_id = "sample01",
    processing = c("filtered", "raw"),
    images = c("lowres", "hires", "detected", "aligned"),
    jsonFile = .SCALE_JSON_FILE,
    tissuePattern = "tissue_positions.*\\.csv",
    spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    ...
) {
    images <- match.arg(images, several.ok = TRUE)
    processing <- match.arg(processing)
    if (!missing(spacerangerOut)) {
        if (isScalarCharacter(spacerangerOut))
            stopifnot(
                dir.exists(spacerangerOut)
            )
        resources <- .find_convert_resources(spacerangerOut, processing, ...)
        spatialResource <- .find_convert_spatial(
            path = spacerangerOut, sample_id = sample_id, images = images,
            jsonFile = jsonFile, tissuePattern = tissuePattern
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

    .TENxVisium(
        resources = resources,
        spatialList = spatialResource,
        coordNames = spatialCoordsNames,
        sampleId = sample_id
    )
}

.validTENxVisium <- function(object) {
    isFL <- is(object@resources, "TENxFileList")
    isSL <- is(object@spatialList, "TENxSpatialList")
    if (isSL && isFL)
        TRUE
    else if (!isFL)
        "'TENxFileList' component is not of TENxFileList class"
    else
        "'TENxSpatialList' component is not of TENxSpatialList class"
}

S4Vectors::setValidity2("TENxVisium", .validTENxVisium)

#' @describeIn TENxVisium-class Import Visium data
#'
#' @inheritParams BiocIO::import
#'
#' @importFrom BiocIO import
#' @importFrom methods as
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom SummarizedExperiment assays rowData colData
#'
#' @exportMethod import
setMethod("import", "TENxVisium", function(con, format, text, ...) {
    sce <- import(con@resources)
    slist <- import(con@spatialList)
    img <- slist[["imgData"]]
    spd <- slist[["colData"]]
    is_tbl_df <- inherits(spd, "tbl_df")
    if (is_tbl_df)
        rownames <- spd[["barcode"]]
    else
        rownames <- rownames(spd)
    matches <- intersect(
        colnames(sce),
        rownames
    )

    if (!length(matches))
        stop("No matching barcodes found between spatial and expression data.")

    if (is_tbl_df)
        spd <- spd[spd[["barcode"]] %in% matches, ]
    else
        spd <- spd[matches, ]

    sce <- sce[, matches]

    SpatialExperiment::SpatialExperiment(
        assays = assays(sce),
        rowData = S4Vectors::DataFrame(Symbol = rowData(sce)[["Symbol"]]),
        sample_id = con@sampleId,
        colData = as(spd, "DataFrame"),
        spatialCoordsNames = con@coordNames,
        imgData = img
    )
})
