#' @docType class
#'
#' @title Represent and import spatial CSV data from 10X Genomics
#'
#' @description `TENxSpatialCSV` is a class to represent and import spatial CSV
#'   files with specific column names. It is a composed class of
#'   [TENxIO::TENxFile] and contains additional slots for the column names and
#'   whether the CSV is a list-type of file.
#'
#' @details Typically, the user will not create an object of this class directly
#'   but rather use the [TENxVisium()] constructor function to create an object
#'   of this class in the background. The column names are set to the default
#'   values of `c("barcode", "in_tissue", "array_row", "array_col",
#'   "pxl_row_in_fullres", "pxl_col_in_fullres")`. The column names can be
#'   changed by specifying the `colnames` argument in the constructor function.
#'
#' @slot isList `logical(1)` A scalar specifying whether the CSV is a list-type
#'   of file
#'
#' @slot colnames `character()` A vector specifying the column names of the CSV
#'
#' @importClassesFrom TENxIO TENxFile
#' @importFrom methods new is
#'
#' @exportClass TENxSpatialCSV
.TENxSpatialCSV <- setClass(
    Class = "TENxSpatialCSV",
    contains = "TENxFile",
    slots = c(
        isList = "logical",
        colnames = "character"
    )
)

.TISSUE_POS_COLS <- c(
    "barcode", "in_tissue", "array_row", "array_col",
    "pxl_row_in_fullres", "pxl_col_in_fullres"
)

#' @rdname TENxSpatialCSV-class
#'
#' @inheritParams TENxIO::TENxFile
#'
#' @param colnames `character()` A vector specifying the column names of the
#'   CSV, defaults to `c("barcode", "in_tissue", "array_row", "array_col",
#'   "pxl_row_in_fullres", "pxl_col_in_fullres")`.
#'
#' @importFrom TENxIO TENxFile
#' @importFrom BiocGenerics path
#'
#' @return TENxSpatialCSV: An object of class [TENxSpatialCSV]
#'
#' @examples
#' sample_dir <- system.file(
#'     file.path("extdata", "10xVisium", "section1"),
#'     package = "SpatialExperiment"
#' )
#' spatial_dir <- Filter(
#'   function(x) endsWith(x, "spatial"), list.dirs(sample_dir)
#' )
#' csvresource <- file.path(spatial_dir, "tissue_positions_list.csv")
#' TENxSpatialCSV(csvresource)
#' head(import(TENxSpatialCSV(csvresource)), 4)
#'
#' @export
TENxSpatialCSV <- function(resource, colnames = .TISSUE_POS_COLS) {
    if (!is(resource, "TENxFile"))
        resource <- TENxFile(resource)
    isList <- grepl("_list", path(resource), fixed = TRUE)
    .TENxSpatialCSV(
        resource, isList = isList, colnames = colnames
    )
}

#' @rdname TENxSpatialCSV-class
#'
#' @inheritParams BiocIO::import
#'
#' @return import-method: A `DataFrame` object containing the data from the CSV
#'   file
#'
#' @importFrom S4Vectors DataFrame
#' @exportMethod import
setMethod("import", "TENxSpatialCSV", function(con, format, text, ...) {
    dat <- utils::read.csv(
        path(con),
        header = !con@isList,
        row.names = 1L,
        col.names = con@colnames
    )
    DataFrame(dat)
})
