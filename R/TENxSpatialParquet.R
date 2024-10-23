#' @title Represent and import spatial Parquet data from 10X Genomics
#'
#' @description `TENxSpatialParquet` is a class to represent and import spatial
#'   Parquet files with specific column names. It is a composed class of
#'   [TENxIO::TENxFile] and contains additional slots for the column names and
#'   whether the Parquet is a list-type of file.
#'
#' @details Typically, the user will not create an object of this class directly
#'   but rather use the [TENxVisium()] constructor function to create an object
#'   of this class in the background. The column names are set to the default
#'   values of `c("barcode", "in_tissue", "array_row", "array_col",
#'   "pxl_row_in_fullres", "pxl_col_in_fullres")`. The column names can be
#'   changed by specifying the `colnames` argument in the constructor function.
#'
#' @slot colnames `character()` A vector specifying the column names of the
#'   Parquet
#'
#' @importClassesFrom TENxIO TENxFile
#' @importFrom methods new is
#'
#' @exportClass TENxSpatialParquet
.TENxSpatialParquet <- setClass(
    Class = "TENxSpatialParquet",
    contains = "TENxFile",
    slots = c(
        colnames = "character"
    )
)

#' @rdname TENxSpatialParquet-class
#'
#' @inheritParams TENxIO::TENxFile
#'
#' @param colnames `character()` A vector specifying the column names of the
#'   Parquet, defaults to `c("barcode", "in_tissue", "array_row", "array_col",
#'   "pxl_row_in_fullres", "pxl_col_in_fullres")`.
#'
#' @return `TENxSpatialParquet()`: An object of class [TENxSpatialParquet]
#'
#' @examples
#' sample_dir <- system.file(
#'     file.path("extdata", "binned_outputs", "square_002um", "spatial"),
#'     package = "VisiumIO"
#' )
#' spatial_dir <- Filter(
#'   function(x) endsWith(x, "spatial"), list.dirs(sample_dir)
#' )
#' parquetres <- file.path(spatial_dir, "tissue_positions.parquet")
#' TENxSpatialParquet(parquetres)
#' import(TENxSpatialParquet(parquetres))
#'
#' ## metadata in attributes
#' import(TENxSpatialParquet(parquetres)) |>
#'     attr("metadata") |>
#'     lapply(names)
#' @export
TENxSpatialParquet <- function(resource, colnames = .TISSUE_POS_COLS) {
    if (!is(resource, "TENxFile"))
        resource <- TENxFile(resource)
    .TENxSpatialParquet(
        resource, colnames = colnames
    )
}

#' @rdname TENxSpatialParquet-class
#'
#' @inheritParams BiocIO::import
#'
#' @importFrom BiocBaseUtils checkInstalled
#'
#' @return import-method: A `DataFrame` object containing the data from the
#'   Parquet file
#'
#' @exportMethod import
setMethod("import", "TENxSpatialParquet", function(con, format, text, ...) {
    checkInstalled("arrow")
    res <- arrow::read_parquet(
        path(con)
    )
    attr(res, "metadata") <- metadata(con)
    res
})
