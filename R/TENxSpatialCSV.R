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
#' @slot variant `character(1)` A scalar specifying the variant of the CSV file
#'   "positions", "cell_boundaries", or "other". The variant is determined by
#'   the name of the CSV file within the constructor function. Values include
#'   "positions", "cell_boundaries", and "other".
#'
#' @slot compressed `logical(1)` A scalar specifying whether the CSV is
#'   compressed (mainly with a `.gz` file extension).
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
        colnames = "character",
        variant = "character",
        compressed = "logical"
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
#'   "pxl_row_in_fullres", "pxl_col_in_fullres")`. Mainly used for the
#'   "positions" CSV type of file which does not include column names in the
#'   file.
#'
#' @details Set the option "VisiumIO.csvreader" to either "data.table" or
#'   "readr" to use the `data.table::fread` or `readr::read_csv` functions,
#'   respectively. These options are useful when the CSV file is relatively
#'   large and the user wants to use faster read-in options. Note that the
#'   outputs will still be converted to `DataFrame` when incorporated to the
#'   `SpatialExperiment` or `SingleCellExperiment` object.
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
    variant <- "other"
    filename <- basename(path(resource))
    if (grepl("positions", fixed = TRUE, filename))
        variant <- "positions"
    else if (grepl("cell_boundaries", fixed = TRUE, filename))
        variant <- "cell_boundaries"

    isCompressed <- endsWith(path(resource), "gz")
    isList <- grepl("_list", path(resource), fixed = TRUE)
    .TENxSpatialCSV(
        resource, isList = isList, variant = variant,
        compressed = isCompressed, colnames = colnames
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
    args <- list(file = path(con), header = !con@isList, row.names = 1L)
    if (identical(con@variant, "positions"))
        args <- c(args, list(col.names = con@colnames))
    else if (identical(con@variant, "cell_boundaries"))
        args <- args[names(args) != "row.names"]
    if (
        identical(getOption("VisiumIO.csvreader"), "data.table") &&
        checkInstalled("data.table")
    )
        dat <- do.call(data.table::fread, args[names(args) != "row.names"])
    else if (
        identical(getOption("VisiumIO.csvreader"), "readr") &&
        checkInstalled("readr")
    )
        dat <- do.call(
            readr::read_csv,
            list(
                col_names = args[["col.names"]],
                file = args[["file"]],
                show_col_types = FALSE
            )
        )
    else
        dat <- do.call(utils::read.csv, args)
    res <- DataFrame(dat)
    if ("barcode" %in% colnames(res)) {
        rownames(res) <- res[["barcode"]]
        res <- res[, names(res) != "barcode", drop = FALSE]
    }
    res
})
