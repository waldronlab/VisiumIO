#' @exportClass TENxGeoJSON
.TENxGeoJSON <- setClass(
    Class = "TENxGeoJSON",
    contains = "TENxFile"
)

#' @export
TENxGeoJSON <- function(resource) {
    if (!is(resource, "TENxFile"))
        resource <- TENxFile(resource)
    .TENxGeoJSON(
        resource
    )
}

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

