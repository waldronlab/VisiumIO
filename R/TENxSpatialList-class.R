#' @docType class
#'
#' @title A class to represent and import spatial Visium data
#'
#' @description This class is a composed class of [TENxFileList], which can
#'   contain a list of [TENxFile] objects, and a [TENxSpatialList] object. It is
#'   meant to handle spatial Visium data from 10X Genomics.
#'
#' @details Typically, the user will not create an object of this class directly
#'   but rather use the [TENxVisium()] constructor function to create an object
#'   of this class.
#'
#' @inheritParams TENxVisium
#'
#' @exportClass TENxSpatialList
.TENxSpatialList <- setClass(
    "TENxSpatialList",
    contains = "TENxFileList",
    slots = c(
        images = "character",
        scaleJSON = "character",
        tissuePos = "character",
        sampleId = "character"
    )
)

.check_file_pattern <- function(obj, pattern) {
    fname <- switch(
        pattern,
        "positions.*\\.csv$" = "tissue positions",
        "scalefactors.*\\.json$" = "scalefactor JSON"
    )
    if (any(grepl(pattern, names(obj))))
        TRUE
    else
        paste0("The '", fname, "' file was not found")
}

.check_file <- function(obj, filename) {
    if (filename %in% names(obj))
        TRUE
    else
        paste0("The '", filename, "' file was not found")
}

.validTENxSpatialList <- function(object) {
    .check_file_pattern(object, "positions.*\\.csv$")
    .check_file_pattern(object, "scalefactors.*\\.json$")
    .check_file(object, object@scaleJSON)
}

S4Vectors::setValidity2("TENxSpatialList", .validTENxSpatialList)

.SCALE_JSON_FILE <- "scalefactors_json.json"

#' @rdname TENxSpatialList-class
#'
#' @inheritParams TENxVisium
#'
#' @importFrom BiocIO decompress
#'
#' @returns A `SpatialExperiment` object
#'
#' @examples
#'
#' spatial_dir <- system.file(
#'     file.path("extdata", "10xVisium", "section1", "outs", "spatial"),
#'     package = "SpatialExperiment"
#' )
#'
#' TENxSpatialList(resources = spatial_dir)
#'
#' @export
TENxSpatialList <- function(
    resources,
    sample_id = "sample01",
    images = c("lowres", "hires", "detected", "aligned"),
    jsonFile = .SCALE_JSON_FILE,
    tissuePattern = "tissue_positions.*\\.csv",
    ...
) {
    images <- match.arg(images, several.ok = TRUE)
    spatf <- TENxFileList(resources, ...)
    if (spatf@compressed)
        spatf <- decompress(con = spatf)
    tissuePos <- grep(tissuePattern, names(spatf), value = TRUE)
    if (!length(tissuePos))
        stop("No tissue positions file found with pattern: ", tissuePattern)

    .TENxSpatialList(
        spatf, images = images, scaleJSON = jsonFile,
        tissuePos = tissuePos, sampleId = sample_id
    )
}

#' @describeIn TENxSpatialList Import a `TENxSpatialList` object
#'
#' @param ... Parameters to pass to the format-specific method.
#'
#' @inheritParams BiocIO::import
#'
#' @exportMethod import
setMethod("import", "TENxSpatialList", function(con, format, text, ...) {
    jsonFile <- con@scaleJSON
    sampid <- con@sampleId
    sfs <- jsonlite::fromJSON(txt = path(con)[jsonFile])

    DFs <- lapply(con@images, function(image) {
        .getImgRow(con = con, sampleId = sampid, image = image, scaleFx = sfs)
    })
    list(
        imgData = DataFrame(
            do.call(rbind, DFs)
        ),
        colData = import(
            TENxSpatialCSV(
                path(con)[con@tissuePos]
            )
        )
    )
})

.getImgRow <- function(con, sampleId, image, scaleFx) {
    scfactor <- NA_integer_
    fileNames <- names(con)
    filePaths <- path(con)
    imgFile <- grep(image, fileNames, value = TRUE)
    imgPath <- filePaths[endsWith(filePaths, imgFile)]
    if (!length(imgPath))
        stop(
            "The '", image, "' image was not found in the list of file names.",
            call. = FALSE
        )
    spi <- SpatialExperiment::SpatialImage(imgPath)
    scaleName <- grep(image, names(scaleFx), value = TRUE)
    if (length(scaleName))
        scfactor <- unlist(scaleFx[scaleName])

    S4Vectors::DataFrame(
        sample_id = sampleId,
        image_id = image,
        data = I(list(spi)),
        scaleFactor = scfactor
    )
}
