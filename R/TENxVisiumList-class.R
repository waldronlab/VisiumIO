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
#' @exportClass TENxVisiumList
.TENxVisiumList <- setClass(
    Class = "TENxVisiumList",
    slots = c(
        VisiumList = "SimpleList"
    )
)

#' @rdname TENxVisiumList-class
#'
#' @inheritParams TENxVisium
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
#'         package = "VisiumIO"
#'     ),
#'     recursive = FALSE, full.names = TRUE
#' )
#'
#' tvl <- TENxVisiumList(
#'     sampleFolders = sample_dirs,
#'     sample_ids = c("sample01", "sample02"),
#'     processing = "raw",
#'     images = "lowres",
#'     format = "mtx"
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
    format = c("mtx", "h5"),
    jsonFile = .SCALE_JSON_FILE,
    tissuePattern = "tissue_positions.*\\.csv",
    spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    ...
) {
    images <- match.arg(images, several.ok = TRUE)
    processing <- match.arg(processing)

    resources <- lapply(
        sampleFolders,
        .find_convert_resources,
        processing = processing,
        format = format,
        ...
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
