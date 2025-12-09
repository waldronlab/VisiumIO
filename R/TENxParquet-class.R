.MAPPING_POS_COLS <- c(
    "square_002um", "square_008um", "square_016um",
    "cell_id", "in_nucleus", "in_cell"
)

#' @docType class
#'
#' @title Represent and import Parquet data from 10X Genomics
#'
#' @description `TENxParquet` is a virtual class to represent and import Parquet
#'   files from 10X Genomics. It is a composed class of [TENxIO::TENxFile].
#'
#' @details Typically, the user will not create an object of this class directly
#'   but rather use the `TENxParquet` constructor function to create an object
#'   of either [TENxSpatialParquet] or [TENxMappingParquet] class in the
#'   background.
#'
#' @importClassesFrom TENxIO TENxFile
#' @importFrom methods new is
#'
#' @exportClass TENxParquet
.TENxParquet <- setClass(
    Class = "TENxParquet",
    contains = c("TENxFile", "VIRTUAL"),
    slots = c(
        colnames = "character"
    )
)

#' @rdname TENxParquet-class
#'
#' @inheritParams TENxIO::TENxFile
#'
#' @param type `character(1)` A scalar specifying the type of Parquet file to
#'   import, either "spatial" or "mapping". If not provided, the type is
#'   determined by whether or not there is a "mapping" keyword in the file name.
#'
#' @return `TENxParquet()`: An object of class [TENxSpatialParquet] or
#'   [TENxMappingParquet]
#'
#' @examples
#' sample_dir <- system.file(
#'     file.path("extdata", "binned_outputs", "square_002um", "spatial"),
#'     package = "VisiumIO"
#' )
#' spatial_dir <- Filter(
#'     function(x) endsWith(x, "spatial"), list.dirs(sample_dir)
#' )
#' parquetres <- file.path(spatial_dir, "tissue_positions.parquet")
#' TENxParquet(parquetres)
#' @export
TENxParquet <- function(resource, type = c("spatial", "mapping")) {
    if (missing(type)) {
        fname <- tolower(basename(resource))
        hasMapping <- grepl("mapping", fname, TRUE)
        type <- if (hasMapping) "mapping" else "spatial"
    } else {
        stopifnot(
            isScalarCharacter(type), type %in% c("spatial", "mapping")
        )
    }
    if (!is(resource, "TENxFile"))
        resource <- TENxFile(resource)

    CFUN <- switch(
        type,
        spatial = TENxSpatialParquet,
        mapping = TENxMappingParquet
    )
    CFUN(resource)
}

#' @rdname TENxParquet-class
#'
#' @description `TENxSpatialParquet` is a class to represent and import spatial
#'   Parquet files from 10X Genomics. It is a composed class of
#'   [TENxIO::TENxFile].
#'
#' @exportClass TENxSpatialParquet
.TENxMappingParquet <- setClass(
    Class = "TENxMappingParquet",
    contains = "TENxParquet"
)

#' @rdname TENxParquet-class
#'
#' @inheritParams TENxIO::TENxFile
#' @inheritParams TENxSpatialParquet
#'
#' @return `TENxMappingParquet()`: An object of class [TENxMappingParquet]
#'
#' @examples
#' parq_file <- system.file(
#'     "extdata", "barcode_mappings.parquet",
#'     package = "VisiumIO", mustWork = TRUE
#' )
#' TENxMappingParquet(parq_file)
#' @export
TENxMappingParquet <- function(resource, colnames = .MAPPING_POS_COLS) {
    if (!is(resource, "TENxFile"))
        resource <- TENxFile(resource)
    checkInstalled("arrow")
    cnames <- names(arrow::open_dataset(path(resource)))
    if (!all(colnames %in% cnames)) {
        warning(
            "The provided column names do not match all in the Parquet file.",
            call. = FALSE
        )
        colnames <- intersect(colnames, cnames)
    }
    .TENxMappingParquet(
        resource, colnames = colnames
    )
}

#' @rdname TENxParquet-class
#'
#' @exportMethod import
setMethod("import", "TENxMappingParquet", function(con, format, text, ...) {
    checkInstalled("arrow")
    res <- arrow::read_parquet(
        file = path(con),
        col_select = con@colnames,
        ...
    )
    attr(res, "metadata") <- metadata(con)
    res
})
