#' @title Import 10X Genomics GeoJSON files
#'
#' @description `TENxGeoJSON` is a class to represent and import GeoJSON files
#'   from 10X Genomics. It is a composed class of [TENxIO::TENxFile].
#'
#' @details Typically, the user will not create an object of this class directly
#'   but rather use the [TENxVisium()] constructor function to create an object
#'   of this class in the background.
#'
#' @importClassesFrom TENxIO TENxFile
#' @importFrom methods new is
#'
#' @exportClass TENxGeoJSON
.TENxGeoJSON <- setClass(
    Class = "TENxGeoJSON",
    contains = "TENxFile"
)

#' @rdname TENxGeoJSON-class
#'
#' @inheritParams TENxIO::TENxFile
#'
#' @return `TENxGeoJSON()`: An object of class [TENxGeoJSON]
#'
#' @seealso <https://www.10xgenomics.com/support/software/xenium-ranger/3.0/analysis/segmentation-inputs>
#'
#' @examples
#' segout_folder <- system.file(
#'     file.path("extdata", "segmented_outputs"),
#'     package = "VisiumIO"
#' )
#' geojsonres <- file.path(segout_folder, "cell_segmentations.geojson")
#'
#' TENxGeoJSON(geojsonres)
#'
#' TENxGeoJSON(geojsonres) |>
#'   import()
#' @export
TENxGeoJSON <- function(resource) {
    if (!is(resource, "TENxFile"))
        resource <- TENxFile(resource)
    .TENxGeoJSON(
        resource
    )
}

#' @rdname TENxGeoJSON-class
#'
#' @inheritParams BiocIO::import
#'
#' @importFrom BiocBaseUtils checkInstalled
#'
#' @return import-method: An `sf` and `data.frame` with the GeoJSON data
#'
#' @exportMethod import
setMethod("import", "TENxGeoJSON", function(con, format, text, ...) {
    checkInstalled("sf")
    geo_data <- sf::st_read(
        dsn = path(con),
        quiet = TRUE,
        stringsAsFactors = FALSE,
        ...
    )
    sf::st_crs(geo_data) <- NA
    rownames(geo_data) <- geo_data[["cell_id"]]
    geo_data
})

