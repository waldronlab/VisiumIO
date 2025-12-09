##  Prepare smaller example files from the
## Visium_HD_3prime_Human_Ovarian_Cancer_segmented_outputs
## https://www.10xgenomics.com/datasets/
## visium-hd-three-prime-ovarian-cancer-fresh-frozen

library(sf)

# cell segmentations ------------------------------------------------------
cseg <- file.path(
    "~/data",
    "Visium_HD_3prime_Human_Ovarian_Cancer",
    "Visium_HD_3prime_Human_Ovarian_Cancer_segmented_outputs",
    "segmented_outputs",
    "cell_segmentations.geojson"
)
geo_data <- st_read(dsn = cseg, quiet = TRUE, stringsAsFactors = FALSE)
st_crs(geo_data) <- NA
rownames(geo_data) <- geo_data[["cell_id"]]
geo_data <- geo_data[1:2, ]

st_write(
    geo_data,
    dsn <- file.path(
        "inst", "extdata", "segmented_outputs", "cell_segmentations.geojson"
    )
)

# TENxGeoJSON("inst/extdata/segmented_outputs/cell_segmentations.geojson") |>
#     import()

# nucleus segmentations ---------------------------------------------------
nseg <- file.path(
    "~/data",
    "Visium_HD_3prime_Human_Ovarian_Cancer",
    "Visium_HD_3prime_Human_Ovarian_Cancer_segmented_outputs",
    "segmented_outputs",
    "nucleus_segmentations.geojson"
)
geo_data <- st_read(dsn = nseg, quiet = TRUE, stringsAsFactors = FALSE)
st_crs(geo_data) <- NA
rownames(geo_data) <- geo_data[["cell_id"]]
geo_data <- geo_data[1:2, ]
st_write(
    geo_data,
    dsn <- file.path(
        "inst", "extdata", "segmented_outputs", "nucleus_segmentations.geojson"
    )
)

# TENxGeoJSON("inst/extdata/segmented_outputs/nucleus_segmentations.geojson") |>
#     import()

# Write smaller H5 file for examples --------------------------------------

library(rhdf5)
library(HDF5Array)

h5f <- file.path(
    "~/data",
    "Visium_HD_3prime_Human_Ovarian_Cancer",
    "Visium_HD_3prime_Human_Ovarian_Cancer_segmented_outputs",
    "segmented_outputs",
    "filtered_feature_cell_matrix.h5"
)

h5new <- file.path(
    "inst", "extdata", "segmented_outputs",
    "filtered_feature_cell_matrix.h5"
)

newmat <- TENxMatrix(h5f, "matrix")[1:10, 1:2]
HDF5Array::writeTENxMatrix(newmat, h5new, group = "matrix", verbose = TRUE)


tkeys <- h5read(h5f, "/matrix/features/_all_tag_keys") ## same
myidx <- list(1:10)
h1 <- h5read(h5f, "/matrix/features/feature_type", index = myidx)
h2 <- h5read(h5f, "/matrix/features/genome", index = myidx)
h3 <- h5read(h5f, "/matrix/features/id", index = myidx)
h4 <- h5read(h5f, "/matrix/features/name", index = myidx)

h5createGroup(h5new, "/matrix/features/")
h5write(tkeys, h5new, "/matrix/features/_all_tag_keys")
h5write(h1, h5new, "/matrix/features/feature_type")
h5write(h2, h5new, "/matrix/features/genome")
h5write(h3, h5new, "/matrix/features/id")
h5write(h4, h5new, "/matrix/features/name")

TENxH5(h5new)
file.info(h5new)$size

# Write smaller barcode_mappings.parquet ----------------------------------

park <- arrow::read_parquet(
    file.path(
        "~/data",
        "Visium_HD_3prime_Human_Ovarian_Cancer",
        "Visium_HD_3prime_Human_Ovarian_Cancer_barcode_mappings.parquet"
    )
)
ss <- TENxH5(
    "inst/extdata/segmented_outputs/filtered_feature_cell_matrix.h5",
    ranges = NA_character_
) |> import()
writepq <- park[park$cell_id %in% colnames(ss),]

arrow::write_parquet(
    writepq,
    sink = "inst/extdata/barcode_mappings.parquet"
)
