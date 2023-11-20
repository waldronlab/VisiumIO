
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

| **Extension**  | **Class**       | **Imported as**       |
|----------------|-----------------|-----------------------|
| spatial.tar.gz | TENxSpatialList | inter. DataFrame list |

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

Note that we use the `image = "lowres"` and `processing = "raw"`
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
        file.path("extdata", "10xVisium"), package = "SpatialExperiment"
    ),
    recursive = FALSE, full.names = TRUE
)
    
vlist <- TENxVisiumList(
    sampleFolders = sample_dirs,
    sample_ids = basename(sample_dirs),
    processing = "raw",
    image = "lowres"
)
#> Warning in TENxVisiumList(sampleFolders = sample_dirs, sample_ids =
#> basename(sample_dirs), : partial argument match of 'image' to 'images'
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
