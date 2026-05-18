#' @docType class
#'
#' @title A class to represent and import spatial Visium data
#'
#' @description This class is a composed class of
#'   [TENxFileList][TENxIO::TENxFileList-class], which can contain a list of
#'   [TENxFile][TENxIO::TENxFile-class] objects, and a [TENxSpatialList] object.
#'   It is meant to handle spatial Visium data from 10X Genomics.
#'
#' @details Typically, the user will not create an object of this class directly
#'   but rather use the [TENxVisium()] constructor function to create an object
#'   of this class.
#'
#' @inheritParams TENxVisium
#'
#' @slot images `character()` The image name(s) to use with `grep` and include
#'   in the list of files. Can be one of "lowres", "hires",  "lowres", "hires",
#'   "detected", "aligned", "aligned_fiducials", or "cytassist".
#'
#' @slot scaleJSON `character(1)` The file name of the scale factors JSON file,
#'   defaults to 'scalefactors_json.json'.
#'
#' @slot tissuePos `character(1)` An optional slot indicating the file name of
#'   the tissue positions file; typically a `.parquet` or `.csv` file. To avoid
#'   import, this slot can be set to an empty character value i.e., `""`.
#'
#' @slot sampleId `character(1)` A scalar specifying the sample identifier.
#'
#' @slot binSize `character(1)` An optional slot to store the image bin size
#'   when importing. The default slot value is an empty character i.e., `""`.
#'   When present, the value must be a character scalar, e.g., `'008'` for 8
#'   microns, and will match the directory name e.g., `square_008um`.
#'
#' @slot loadImage `logical(1)` Whether to load the images into memory as
#'   `SpatialImage` objects. If `FALSE`, the images are stored as file paths and
#'   loaded as `StoredSpatialImage` objects when the `TENxSpatialList` object is
#'   imported. The default is `FALSE` to avoid loading large images into memory.
#'   This functionality requires the `magick` package.
#'
#' @exportClass TENxSpatialList
.TENxSpatialList <- setClass(
    "TENxSpatialList",
    contains = "TENxFileList",
    slots = c(
        images = "character",
        scaleJSON = "character",
        tissuePos = "character",
        sampleId = "character",
        binSize = "character",
        loadImage = "logical"
    )
)

.check_file_pattern <- function(obj, pattern) {
    fname <- switch(
        pattern,
        "tissue_positions.*" = "tissue positions",
        "scalefactors.*\\.json$" = "scalefactor JSON"
    )
    if (!any(grepl(pattern, names(obj))))
        paste0("The '", fname, "' file was not found")
}

.check_file <- function(obj, filename) {
    if (!filename %in% names(obj))
        paste0("The '", filename, "' file was not found")
}

.validTENxSpatialList <- function(object) {
    c(
        if (nzchar(object@tissuePos))
            .check_file_pattern(object, "tissue_positions.*"),
        .check_file_pattern(object, "scalefactors.*\\.json$"),
        .check_file(object, object@scaleJSON)
    )
}

S4Vectors::setValidity2("TENxSpatialList", .validTENxSpatialList)

.SCALE_JSON_FILE <- "scalefactors_json.json"

#' @rdname TENxSpatialList-class
#'
#' @inheritParams TENxVisium
#'
#' @param bin_size `character(1)` An _optional_ scalar indicating the image bin
#'   size in microns, e.g., `'008'` for 8 microns. When provided, the function
#'   will look for a subfolder in the `binned_outputs` folder that corresponds
#'   to the specified bin size, e.g., `square_008um` for an input value of
#'   `'008'`. Bin sizes `'002'`, `'008'`, and `'016'` are typical of the space
#'   ranger pipeline but custom bin sizes may be used. If the `bin_size`
#'   argument is not provided, it is set to an empty character value i.e., `""`
#'   and the function will not look for any data in `binned_outputs` and will
#'   only import data from the main `spacerangerOut` directory.
#'
#' @importFrom BiocIO decompress
#'
#' @returns A `SpatialExperiment` object
#'
#' @examples
#' spatial_dir <- system.file(
#'     file.path("extdata", "10xVisium", "section1", "outs", "spatial"),
#'     package = "VisiumIO"
#' )
#'
#' TENxSpatialList(resources = spatial_dir, images = "lowres")
#'
#' TENxSpatialList(resources = spatial_dir, images = "lowres") |>
#'     metadata() |> lapply(names)
#'
#' TENxSpatialList(resources = spatial_dir, images = "lowres") |>
#'     import()
#'
#'
#' @export
TENxSpatialList <- function(
    resources,
    sample_id = "sample01",
    images = c(
        "lowres", "hires", "detected", "aligned",
        "aligned_fiducials", "cytassist"
    ),
    jsonFile = .SCALE_JSON_FILE,
    tissuePattern = "tissue_positions.*",
    bin_size = c("008", "016", "002"),
    loadImage = FALSE,
    ...
) {
    images <- match.arg(images, several.ok = TRUE)
    if (!missing(bin_size)) {
        if (!isScalarCharacter(bin_size, zchar = TRUE))
            stop("The 'bin_size' argument must be a single character value.")
        if (nzchar(bin_size) && !grepl("^\\d{3}$", bin_size))
            stop("The 'bin_size' argument must be a 3-digit string.")
    }

    if (!is(resources, "TENxFileList"))
        resources <- TENxFileList(resources, ...)
    if (resources@compressed)
        resources <- decompress(con = resources)
    tissuePos <- tissuePattern
    if (nzchar(tissuePattern)) {
        tissuePos <- grep(tissuePattern, names(resources), value = TRUE)
        if (!length(tissuePos))
            stop(
                "No tissue positions file found with pattern: ", tissuePattern
            )
    }

    if (missing(bin_size) && any(grepl("square_\\d{3}", path(resources))))
        bin_size <- unique(
            gsub(".*?square_(\\d{3}).*", "\\1", path(resources))
        )
    else if (missing(bin_size))
        bin_size <- ""

    if (nzchar(bin_size) && !bin_size %in% c("002", "008", "016"))
        message("Using custom 'bin_size': ", bin_size)

    .TENxSpatialList(
        resources, images = images, scaleJSON = jsonFile,
        tissuePos = tissuePos, sampleId = sample_id,
        binSize = bin_size, loadImage = loadImage
    )
}

#' @describeIn TENxSpatialList Import a `TENxSpatialList` object
#'
#' @param ... Parameters to pass to the format-specific method.
#'
#' @inheritParams BiocIO::import
#'
#' @importFrom BiocIO FileForFormat
#'
#' @exportMethod import
setMethod("import", "TENxSpatialList", function(con, format, text, ...) {
    jsonFile <- con@scaleJSON
    sampid <- con@sampleId
    sfs <- jsonlite::fromJSON(txt = path(con)[jsonFile])

    DFs <- lapply(
        X = con@images,
        FUN = .getImgRow,
        con = con,
        sampleId = sampid,
        scaleFx = sfs
    )
    res <- list(
        imgData = DataFrame(
            do.call(rbind, DFs)
        )
    )
    if (nzchar(con@tissuePos)) {
        fff <- FileForFormat(
            path(con)[con@tissuePos],
            prefix = "TENxSpatial", suffix = NULL
        )
        ffcolData <- import(fff)
        if (nzchar(con@binSize))
            ffcolData[["bin_size"]] <- con@binSize
        ffcolData <- as(ffcolData, "DataFrame")
        if (length(ffcolData[["barcode"]]))
            rownames(ffcolData) <- ffcolData[["barcode"]]
        res <- c(res, colData = ffcolData)
    }
    res
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
    img <- imgPath

    if (con@loadImage) {
        checkInstalled("magick")
        img <- magick::image_read(imgPath) |> grDevices::as.raster()
    }

    spi <- SpatialExperiment::SpatialImage(img)

    if (identical(image, "cytassist"))
        image <- "regist_target"
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
