
# Introduction

The `VisiumIO` package provides a set of functions to import 10X
Genomics Visium experiment data into a `SpatialExperiment` object. The
package makes use of the `SpatialExperiment` data structure, which
provides a set of classes and methods to handle spatially resolved
transcriptomics data.

# TENxIO Supported Formats

| **Extension**       | **Class**     | **Imported as**                    |
|---------------------|---------------|------------------------------------|
| .h5                 | TENxH5        | SingleCellExperiment w/ TENxMatrix |
| .mtx / .mtx.gz      | TENxMTX       | SummarizedExperiment w/ dgCMatrix  |
| .tar.gz             | TENxFileList  | SingleCellExperiment w/ dgCMatrix  |
| peak_annotation.tsv | TENxPeaks     | GRanges                            |
| fragments.tsv.gz    | TENxFragments | RaggedExperiment                   |
| .tsv / .tsv.gz      | TENxTSV       | tibble                             |

# VisiumIO Supported Formats

| **Extension**  | **Class**          | **Imported as**   |
|----------------|--------------------|-------------------|
| spatial.tar.gz | TENxSpatialList    | DataFrame list \* |
| .parquet       | TENxSpatialParquet | tibble \*         |

**Note**. (\*) Intermediate format

# Installation

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("VisiumIO")
```

# Loading package

``` r
library(VisiumIO)
```

# TENxVisium

The `TENxVisium` class is used to import a **single** sample of 10X
Visium data. The `TENxVisium` constructor function takes the following
arguments:

``` r
TENxVisium(
    resources = "path/to/10x/visium/file.tar.gz",
    spatialResource = "path/to/10x/visium/spatial/file.spatial.tar.gz",
    spacerangerOut = "path/to/10x/visium/sample/folder",
    sample_id = "sample01",
    images = c("lowres", "hires", "detected", "aligned"),
    jsonFile = "scalefactors_json.json",
    tissuePattern = "tissue_positions.*\\.csv",
    spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres")
)
```

The `resource` argument is the path to the 10X Visium file. The
`spatialResource` argument is the path to the 10X Visium spatial file.
It usually ends in `spatial.tar.gz`.

## Example from SpatialExperiment

Note that we use the `images = "lowres"` and `processing = "raw"`
arguments based on the name of the `tissue_*_image.png` file and
`*_feature_bc_matrix` folder in the `spaceranger` output. The directory
structure for a **single** sample is shown below:

        section1
        └── outs
            ├── spatial
            │   ├── tissue_lowres_image.png
            │   └── tissue_positions_list.csv
            └── raw_feature_bc_matrix
                ├── barcodes.tsv
                ├── features.tsv
                └── matrix.mtx

### Creating a TENxVisium instance

Using the example data in `SpatialExperiment`, we can load the
`section1` sample using `TENxVisium`.

``` r
sample_dir <- system.file(
    file.path("extdata", "10xVisium", "section1"),
    package = "SpatialExperiment"
)

vis <- TENxVisium(
    spacerangerOut = sample_dir, processing = "raw", images = "lowres"
)
vis
#> An object of class "TENxVisium"
#> Slot "resources":
#> TENxFileList of length 3
#> names(3): barcodes.tsv features.tsv matrix.mtx
#> 
#> Slot "spatialList":
#> TENxSpatialList of length 3
#> names(3): scalefactors_json.json tissue_lowres_image.png tissue_positions_list.csv
#> 
#> Slot "coordNames":
#> [1] "pxl_col_in_fullres" "pxl_row_in_fullres"
#> 
#> Slot "sampleId":
#> [1] "sample01"
```

The show method of the `TENxVisium` class displays the object’s
metadata.

### Importing into SpatialExperiment

The `TEnxVisium` object can be imported into a `SpatialExperiment`
object using the `import` function.

``` r
import(vis)
#> class: SpatialExperiment 
#> dim: 50 50 
#> metadata(0):
#> assays(1): counts
#> rownames: NULL
#> rowData names(1): Symbol
#> colnames(50): AAACAACGAATAGTTC-1 AAACAAGTATCTCCCA-1 ...
#>   AAAGTCGACCCTCAGT-1 AAAGTGCCATCAATTA-1
#> colData names(4): in_tissue array_row array_col sample_id
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
#> imgData names(4): sample_id image_id data scaleFactor
```

# TENxVisiumList

The `TENxVisiumList` class is used to import multiple samples of 10X
Visium. The interface is a bit more simple in that you only need to
provide the `space ranger` output folder as input to the function.

``` r
TENxVisiumList(
    sampleFolders = "path/to/10x/visium/sample/folder",
    sample_ids = c("sample01", "sample02"),
    ...
)
```

The `sampleFolders` argument is a character vector of paths to the
`spaceranger` output folder. Note that each folder must contain an
`outs` directory. The `sample_ids` argument is a character vector of
sample ids.

## Example from SpatialExperiment

The directory structure for **multiple** samples (`section1` and
`section2`) is shown below:

        section1
        └── outs
        |   ├── spatial
        |   └── raw_feature_bc_matrix
        section2
        └── outs
            ├── spatial
            └── raw_feature_bc_matrix

### Creating a TENxVisiumList

The main inputs to `TENxVisiumList` are the `sampleFolders` and
`sample_ids`. These correspond to the `spaceranger` output sample
folders and a vector of sample identifiers, respectively.

``` r
sample_dirs <- list.dirs(
    system.file(
        file.path("extdata", "10xVisium"), package = "VisiumIO"
    ),
    recursive = FALSE, full.names = TRUE
)
    
vlist <- TENxVisiumList(
    sampleFolders = sample_dirs,
    sample_ids = basename(sample_dirs),
    processing = "raw",
    images = "lowres"
)
vlist
#> An object of class "TENxVisiumList"
#> Slot "VisiumList":
#> List of length 2
```

### Importing into SpatialExperiment

The `import` method combines both `SingleCellExperiment` objects along
with the spatial information into a single `SpatialExperiment` object.
The number of columns in the SpatialExperiment object is equal to the
number of cells across both samples (`section1` and `section2`).

``` r
import(vlist)
#> class: SpatialExperiment 
#> dim: 50 99 
#> metadata(0):
#> assays(1): counts
#> rownames: NULL
#> rowData names(1): Symbol
#> colnames(99): AAACAACGAATAGTTC-1 AAACAAGTATCTCCCA-1 ...
#>   AAAGTCGACCCTCAGT-1 AAAGTGCCATCAATTA-1
#> colData names(4): in_tissue array_row array_col sample_id
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
#> imgData names(4): sample_id image_id data scaleFactor
```

### Visium HD folder structure

The directory structure for a single bin size is shown below.

        Visium_HD
        └── binned_outputs
            └─── square_002um
            │   └── filtered_feature_bc_matrix
            │   │   └── barcodes.tsv.gz
            │   │   └── features.tsv.gz
            │   │   └── matrix.mtx.gz
            │   └── filtered_feature_bc_matrix.h5
            │   └── raw_feature_bc_matrix/
            │   └── raw_feature_bc_matrix.h5
            │   └── spatial
            │       └── [ ... ]
            │       └── tissue_positions.parquet
            └── square_*

### Import Visium HD into SpatialExperiment

``` r
TENxVisiumHD(
    spacerangerOut = "./Visium_HD/",
    sample_id = "sample01",
    processing = c("filtered", "raw"),
    images = c("lowres", "hires", "detected", "aligned_fiducials"),
    bin_size = c("002", "008", "016"),
    jsonFile = .SCALE_JSON_FILE,
    tissuePattern = "tissue_positions\\.parquet",
    spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    ...
)
```

### In-package example

By default, the MatrixMarket format is read in (`format = "mtx"`).

``` r
visfold <- system.file(
    package = "VisiumIO", "extdata", mustWork = TRUE
)
TENxVisiumHD(
    spacerangerOut = visfold, images = "lowres", bin_size = "002"
) |> import()
#> class: SpatialExperiment 
#> dim: 10 10 
#> metadata(2): resources spatialList
#> assays(1): counts
#> rownames(10): ENSMUSG00000051951 ENSMUSG00000025900 ... ENSMUSG00000033774 ENSMUSG00000025907
#> rowData names(3): ID Symbol Type
#> colnames(10): s_002um_02448_01644-1 s_002um_00700_02130-1 ... s_002um_01016_02194-1 s_002um_00775_02414-1
#> colData names(6): barcode in_tissue ... bin_size sample_id
#> reducedDimNames(0):
#> mainExpName: Gene Expression
#> altExpNames(0):
#> spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
#> imgData names(4): sample_id image_id data scaleFactor
```

H5 files are supported via the `format = "h5"` argument input.

``` r
TENxVisiumHD(
    spacerangerOut = visfold, images = "lowres", bin_size = "002",
    format = "h5"
) |> import()
#> class: SpatialExperiment 
#> dim: 10 10 
#> metadata(2): resources spatialList
#> assays(1): counts
#> rownames(10): ENSMUSG00000051951 ENSMUSG00000025900 ... ENSMUSG00000033774 ENSMUSG00000025907
#> rowData names(3): ID Symbol Type
#> colnames(10): s_002um_02448_01644-1 s_002um_00700_02130-1 ... s_002um_01016_02194-1 s_002um_00775_02414-1
#> colData names(6): barcode in_tissue ... bin_size sample_id
#> reducedDimNames(0):
#> mainExpName: Gene Expression
#> altExpNames(0):
#> spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
#> imgData names(4): sample_id image_id data scaleFactor
```

# Session Info

``` r
sessionInfo()
#> R version 4.5.0 Patched (2025-04-15 r88148)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.2 LTS
#> 
#> Matrix products: default
#> BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: America/New_York
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#>  [1] VisiumIO_1.5.1              TENxIO_1.11.1               SingleCellExperiment_1.31.0 SummarizedExperiment_1.39.0
#>  [5] Biobase_2.69.0              GenomicRanges_1.61.0        GenomeInfoDb_1.45.3         IRanges_2.43.0             
#>  [9] S4Vectors_0.47.0            BiocGenerics_0.55.0         generics_0.1.3              MatrixGenerics_1.21.0      
#> [13] matrixStats_1.5.0           colorout_1.3-2             
#> 
#> loaded via a namespace (and not attached):
#>  [1] rjson_0.2.23             xfun_0.52                rhdf5_2.53.0             lattice_0.22-7           tzdb_0.5.0              
#>  [6] rhdf5filters_1.21.0      vctrs_0.6.5              tools_4.5.0              parallel_4.5.0           tibble_3.2.1            
#> [11] pkgconfig_2.0.3          BiocBaseUtils_1.11.0     R.oo_1.27.0              Matrix_1.7-3             assertthat_0.2.1        
#> [16] lifecycle_1.0.4          compiler_4.5.0           codetools_0.2-20         htmltools_0.5.8.1        yaml_2.3.10             
#> [21] pillar_1.10.2            crayon_1.5.3             R.utils_2.13.0           rsconnect_1.3.4          DelayedArray_0.35.1     
#> [26] magick_2.8.6             abind_1.4-8              tidyselect_1.2.1         digest_0.6.37            purrr_1.0.4             
#> [31] arrow_19.0.1.1           fastmap_1.2.0            grid_4.5.0               archive_1.1.12           cli_3.6.5               
#> [36] SparseArray_1.9.0        magrittr_2.0.3           S4Arrays_1.9.0           h5mread_1.1.0            readr_2.1.5             
#> [41] UCSC.utils_1.5.0         bit64_4.6.0-1            rmarkdown_2.29           XVector_0.49.0           httr_1.4.7              
#> [46] bit_4.6.0                R.methodsS3_1.8.2        hms_1.1.3                SpatialExperiment_1.19.0 HDF5Array_1.37.0        
#> [51] evaluate_1.0.3           knitr_1.50               BiocIO_1.19.0            rlang_1.1.6              Rcpp_1.0.14             
#> [56] glue_1.8.0               rstudioapi_0.17.1        vroom_1.6.5              jsonlite_2.0.0           R6_2.6.1                
#> [61] Rhdf5lib_1.31.0
```
