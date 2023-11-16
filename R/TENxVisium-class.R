#' @docType class
#'
#' @title A class to represent and import a Single Visium Sample
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

#' @importFrom TENxIO TENxFileList
.find_convert_resources <- function(path, processing, ...) {
    odir <- list.dirs(path, recursive = FALSE, full.names = TRUE)
    fdirname <- paste0(processing, "_feature_bc_matrix")
    featdir <- file.path(odir, fdirname)
    if (endsWith(odir, "outs") && dir.exists(featdir))
        path <- featdir
    else
        stop(
            "The 'outs' or '", fdirname, "' directory was not found.",
            "\n  Verify 'spacerangerOut' and 'processing' inputs.",
            call. = FALSE
        )

    if (!is(path, "TENxFileList"))
        resources <- TENxFileList(path, ...)
    else
        resources <- path

    resources
}

.find_convert_spatial <- function(path, ...) {
    odir <- list.dirs(path, recursive = FALSE, full.names = TRUE)
    if (endsWith(odir, "outs"))
        path <- file.path(odir, "spatial")
    else
        stop("The 'outs' directory was not found")

    if (!is(path, "TENxSpatialList"))
        spatialList <- TENxSpatialList(path, ...)
    else
        spatialList <- path

    spatialList
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
        stopifnot(
            isScalarCharacter(spacerangerOut), dir.exists(spacerangerOut)
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
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom SummarizedExperiment assays rowData colData
#'
#' @exportMethod import
setMethod("import", "TENxVisium", function(con, format, text, ...) {
    sce <- import(con@resources)
    slist <- import(con@spatialList)
    img <- slist[["imgData"]]
    spd <- slist[["colData"]]
    matches <- intersect(
        colnames(sce),
        rownames(spd)
    )
    spd <- spd[matches, ]
    sce <- sce[, matches]

    SpatialExperiment::SpatialExperiment(
        assays = assays(sce),
        rowData = S4Vectors::DataFrame(Symbol = rowData(sce)[["Symbol"]]),
        sample_id = con@sampleId,
        colData = spd,
        spatialCoordsNames = con@coordNames,
        imgData = img
    )
})

#' @docType class
#'
#' @title A class to represent and import multiple Visium samples
#'
#' @description This class contains a `SimpleList` of [TENxVisium] objects each
#'   corresponding to one sample.
#'
#' @details Typically, the user will provide a path to a directory containing
#'   the output of the `spaceranger count` command. The `spaceranger count`
#'   command outputs a folder containing the "raw" or "filtered"
#'   `()_feature_bc_matrix`.
#'
#' @inheritParams TENxVisium
#'
#' @return A [SpatialExperiment] object
#'
#' @seealso <https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/overview>
#'
#' @exportClass TENxVisium
.TENxVisiumList <- setClass(
    Class = "TENxVisiumList",
    slots = c(
        VisiumList = "SimpleList"
    )
)

#' @rdname TENxVisiumList-class
#'
#' @param sampleFolders `character()` A vector of strings specifying the
#'   directories containing the output of the `spaceranger count` command.
#'
#' @param sample_ids `character()` A vector of strings specifying the sample
#'   IDs. If not provided, the sample IDs will be the names of the
#'   `sampleFolders`. Therefore, the `sample_ids` must be the same length
#'   as `sampleFolders`.
#'
#' @importFrom S4Vectors SimpleList
#'
#' @examples
#'
#' sample_dirs <- list.dirs(
#'     system.file(
#'         file.path("extdata", "10xVisium"),
#'         package = "SpatialExperiment"
#'     ),
#'     recursive = FALSE, full.names = TRUE
#' )
#'
#' tvl <- TENxVisiumList(
#'     sampleFolders = sample_dirs,
#'     sample_ids = c("sample01", "sample02"),
#'     processing = "raw",
#'     images = "lowres"
#' )
#'
#' import(tvl)
#'
#' @export
TENxVisiumList <- function(
    sampleFolders,
    sample_ids,
    processing = c("filtered", "raw"),
    images = c("lowres", "hires", "detected", "aligned"),
    jsonFile = .SCALE_JSON_FILE,
    tissuePattern = "tissue_positions.*\\.csv",
    spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    ...
) {
    resources <- lapply(
        sampleFolders,
        .find_convert_resources,
        processing = processing, ...
    )

    if (missing(sample_ids))
        sample_ids <- basename(sampleFolders)

    spatialResources <- Map(
        f = .find_convert_spatial,
        path = sampleFolders,
        sample_id = sample_ids,
        MoreArgs = list(
            images = images,
            jsonFile = jsonFile,
            tissuePattern = tissuePattern
        )
    )
    txvList <- Map(
        f = .TENxVisium,
        resources = resources,
        spatialList = spatialResources,
        sampleId = sample_ids,
        MoreArgs = list(
            coordNames = spatialCoordsNames
        )
    )

    .TENxVisiumList(
        VisiumList = SimpleList(txvList)
    )
}

#' @describeIn TENxVisiumList-class Import multiple Visium samples
#'
#' @inheritParams BiocIO::import
#' @importFrom BiocGenerics cbind
#'
#' @exportMethod import
setMethod("import", "TENxVisiumList", function(con, format, text, ...) {
    SElist <- lapply(con@VisiumList, import)
    do.call(cbind, SElist)
})
