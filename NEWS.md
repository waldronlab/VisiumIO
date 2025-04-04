## Changes in version 1.4.0

### New features

* Added support for CytAssist images (@ZheFrench, #8)
* Enabled use of alternative readers via `VisiumIO.csvreader` option 
* Included metadata in `SpatialExperiment` outputs

### Bug fixes

* Resolved issue where `*_feature_bc_matrix` folder checks were incorrectly
  triggered for non-mtx formats (@estellad, #4)
* Fixed missing package anchor in `@param` documentation entries

## Changes in version 1.2.0

### New features

* Support VisiumHD file formats including `parquet`.
* Include `format` argument for `h5` file imports (default remains `mtx`).
* Add example data for Visium and VisiumHD imports.

## Changes in version 1.0.0

* Package released in Bioconductor!
