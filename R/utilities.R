bcInFiltered <- function(
    raw_filt_resources, spacerangerOut, processing = c("raw", "filtered"), outs = "outs", barcodes
) {
    if (!identical(processing, c("raw", "filtered")))
        stop("'processing' types must be 'raw' and 'filtered' for comparison")
    res <- structure(vector("list", 2L), .Names = processing)
    if (!missing(spacerangerOut)) {
        if (isScalarCharacter(spacerangerOut))
            stopifnot(dir.exists(spacerangerOut))
        for (process in processin) {
            res[[process]] <-
                .find_convert_resources(spacerangerOut, processing, ...)
        }
    } else {
        if (
            !identical(length(raw_filt_resources), 2L) &&
            isCharacter(raw_filt_resources)
        )
            stop("'raw_filt_resources' must be a character vector of length 2")
        for (process in processing) {
            file <- grep(process, raw_filt_resources, value = TRUE)
            txfl <- TENxFileList(file)
            if (txfl@compressed)
                res[[process]] <- decompress(con = txfl)
            else
                res[[process]] <- txfl
        }
    }

    bcode_files <- lapply(
        res, function(x) {
            TENxFile(path(x)[startsWith(names(x), "barcode")])
        }
    )
    bcode_frames <- lapply(bcode_files, import)
    do.call(
        `%in%`, unname(lapply(bcode_frames, unlist))
    )
}
