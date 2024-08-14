
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
have built on top of `rsi` for the intended modeling use case of
geopipeline includes:

- Creating a data frame output for regions of interest based on point
  geometries
- Adding the ability to spatially mosaic neighboring image tiles when
  temporal compositing is not applied for both Sentinel-2 L2A data and
  Sentinel-1 RTC GRD data.

Still to come:

- A quick function for stacking and labeling non-composited images
  covering the same spatial extent.
- A function for summarizing pixel values within a buffer based on a
  user-defined formula.
- Core functions for accessing GEDI L2A and L2B gridded datasets.

## Installation

You can install the development version of geopipeline from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("laurenkwick/geopipeline")
```

## Example

Begin by loading example data. In this example, the region of interest
is the George W Mead Wildlife Area in Wisconsin. The point data includes
fires that recently occurred within the wildlife area.

``` r
library(geopipeline)

sf_point <- sf::st_read("~/ProcessingModule/national_forest_data/mead_fires.shp")
#> Reading layer `mead_fires' from data source 
#>   `/home/lwick/ProcessingModule/national_forest_data/mead_fires.shp' 
#>   using driver `ESRI Shapefile'
#> Simple feature collection with 77 features and 25 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 524429.1 ymin: 457321 xmax: 541478.1 ymax: 472452.7
#> Projected CRS: NAD83(HARN) / Wisconsin Transverse Mercator
sf_poly <- sf::st_read("~/ProcessingModule/national_forest_data/mead_boundary.shp", type = 3)
#> Reading layer `mead_boundary' from data source 
#>   `/home/lwick/ProcessingModule/national_forest_data/mead_boundary.shp' 
#>   using driver `ESRI Shapefile'
#> converted into: POLYGON
#> Simple feature collection with 1 feature and 5 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 515559.2 ymin: 454316.2 xmax: 542226.3 ymax: 474365.1
#> Projected CRS: NAD83(HARN) / Wisconsin Transverse Mercator
```

We will use s2_process to return a GeoTIFF for the region of interest
and a data frame containing values for our specific points of interest.

``` r
# Polygon data
s2_composite_im <- s2_process(sf_poly, start_dt = "2024-07-01", end_dt = "2024-07-31",
                              layers_sel = c("red", "nir"), 
                              composite_method = "median", idx_names = "NDVI", 
                              file_path = tempfile(pattern="s2_composite_im", tmpdir=tempdir()))
#> 1 indices were found based on criteria.
#> File saved to: /tmp/RtmpXBUTm0/s2_composite_im7134191a97c_2024-07-01_2024-07-31.tif

# Point data
s2_composite_df <- s2_process(sf_point, start_dt = "2024-07-01", end_dt = "2024-07-31",
                              uniqueID = "FIRE_ID",
                              layers_sel = c("red", "nir"), 
                              composite_method = "median", idx_names = "NDVI",
                              file_path = tempfile(pattern="s2_composite_df", tmpdir=tempdir()))
#> No radius defined. Will not apply a buffer to the point data.
#> 1 indices were found based on criteria.
#> File saved to: /tmp/RtmpXBUTm0/s2_composite_df713682d863c.fst
```

Plot the image:

``` r
terra::plot(terra::rast(s2_composite_im))
```

<img src="man/figures/README-plot_composites-1.png" width="100%" />

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
