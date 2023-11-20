testcsv <- tempfile(fileext = ".csv")

expect_error(
    TENxSpatialCSV(testcsv)
)

write.csv(
    x = data.frame(
        barcode = "A",
        in_tissue = TRUE,
        array_row = 1L,
        array_col = 1L,
        pxl_row_in_fullres = 1L,
        pxl_col_in_fullres = 1L
    ),
    file = testcsv,
    row.names = FALSE
)

tcsv <- TENxSpatialCSV(testcsv)

expect_identical(
    tcsv@extension, "csv"
)

expect_true(
    is(tcsv, "TENxSpatialCSV")
)

expect_identical(
    tcsv@colnames, VisiumIO:::.TISSUE_POS_COLS
)

expect_true(
    is(import(tcsv), "DataFrame")
)

expect_false(
    tcsv@isList
)

testcsv <- tempfile(fileext = "_list.csv")

write.csv(
    x = data.frame(
        barcode = "A",
        in_tissue = TRUE,
        array_row = 1L,
        array_col = 1L,
        pxl_row_in_fullres = 1L,
        pxl_col_in_fullres = 1L
    ),
    file = testcsv,
    row.names = FALSE
)

tcsv <- TENxSpatialCSV(testcsv)

expect_true(
    tcsv@isList
)

