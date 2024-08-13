
<!-- README.md is generated from README.Rmd. Please edit that file -->

# geopipeline

<!-- badges: start -->
<!-- badges: end -->

The goal of geopipeline is to streamline the data access and processing
of geospatial datasets, particularly for modeling workflows. geopipeline
makes it easy to produce 2 types of outputs:

1.  **Data Frame** - Produced as an output when a user wants to extract
    the values associated with particular points in a region of
    interest. The data frame can return the pixel values at the exact
    point locations, or the pixel values of all pixels that fall within
    a specified radius of the point locations. The data frame format may
    be useful for training a model.

2.  **GeoTIFF File(s)** - Produced as an output when the user wants to
    extract the values associated with an entire region of interest.
    When a user defines an area of interest as a polygon, the output of
    the core geopipeline functions is a GeoTIFF file of the bounding box
    for the polygon.

There are currently 7 core functions in geopipeline that can return data
frame or GeoTIFF outputs for 7 unique datasets:

1.  `s2_process` - Sentinel-2 L2A
2.  `s1_process` - Sentinel-1 Radiometrically Terrain Corrected
3.  `dem_process` - Copernicus DEM GLO-30
4.  `esa_worldcover_process` - ESA WorldCover
5.  `alos_fnf_process` - ALOS Forest/Non-Forest
6.  `soil_process` - ISRIC Soil Grids
7.  `gedi_process` - GEDI L3

I would like to specifically acknowledge and cite the `rsi` package from
Permian-Global-Research
(<https://github.com/Permian-Global-Research/rsi/tree/main>) which is
the foundation for accessing, processing, and calculating indices for
the first 5 datasets listed above. That package is truly incredible and
has made working on this project much easier! Additional functionality I
have built on top of the `rsi` for the intended modeling use case of
geopipeline includes:

- Creating a data frame output for regions of interest based on point
  geometries
- Adding the ability to spatially mosaic neighboring image tiles when
  temporal compositing is not applied for both Sentinel-2 L2A data and
  Sentinel-1 RTC GRD data.

Still to come:

- A quick function for stacking and labeling non-composited images
- A function for summarizing pixel values within a buffer based on a
  user-defined formula.

## Installation

You can install the development version of geopipeline from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("laurenkwick/geopipeline")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(geopipeline)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
