#' Sentinel-2 L2A Processing
#'
#' @param aoi_layer 
#' @param radius 
#' @param start_dt 
#' @param end_dt 
#' @param uniqueID 
#' @param layers_sel 
#' @param composite_method 
#' @param resample_method 
#' @param keep_SCL 
#' @param apply_mask 
#' @param app_domains 
#' @param idx_names 
#' @param file_path 
#'
#' @return
#' @export
#'
#' @examples
s2_process <- function(aoi_layer, radius=NULL, start_dt, end_dt, uniqueID=NULL, layers_sel=NULL,
                       composite_method=c("median", "mean", "sum", "min", "max", NULL),
                       resample_method=c("bilinear", "average", "rms", "nearest", "gauss", "cube", "cubicspline", "lanczos", "average_magphase", "mode"),
                       keep_SCL=FALSE, apply_mask=TRUE, app_domains=NULL, idx_names=NULL, file_path) {
  
  # geometry & composite / resampling methods
  validate_geometry(aoi_layer)
  validate_method(composite_method, c("sum", "mean", "median", "min", "max"), "composite_method")
  validate_method(resample_method, c("nearest", "average", "rms", "bilinear", "gauss", "cube", "cubicspline",
                                     "lanczos", "average_magphase", "mode"), "resample_method")
  # uniqueID
  use_uniqueID <- NULL
  use_uniqueID <- valid_col_ID(aoi=aoi_layer, column_ID=uniqueID)
  
  # aoi_search for STAC query
  aoi_search <- prepare_aoi(aoi_layer, radius)
  
  # date range
  dates <- valid_dates(start_dt, end_dt)
  
  # Set rsi band mapping object
  s2_asset <- rsi::sentinel2_band_mapping$aws_v1
  s2_vals <- select_layers(layer_list=layers_sel, band_mapping=rsi::sentinel2_band_mapping$aws_v1)
  
  # indices matching s2_vals
  df_indices <- indices_df(band_names=s2_vals, index_application=app_domains, index_names=idx_names)
  
  # add SCL layer and/or masks if specified
  s2_vals <- add_scl_layer(keep_mask=keep_SCL, composite=composite_method[1], layers=s2_vals)
  mask_layer <- NULL
  mask_fxn <- NULL
  if (apply_mask == TRUE) {
    mask_layer <- attr(s2_asset, "mask_band")
    mask_fxn <- attr(s2_asset, "mask_function")
  }
  
  # Initialize image variable
  image <- process_image_data(stac_search_fxn = rsi::get_sentinel2_imagery,
                              aoi=aoi_search,
                              date_range=dates,
                              band_mapping=rsi::sentinel2_band_mapping$aws_v1,
                              search_layers=s2_vals,
                              mask_name=mask_layer,
                              mask_function=mask_fxn,
                              composite=composite_method[1],
                              resample=resample_method,
                              max_date_range=100,
                              date_interval=5)
  
  
  # Initialize the variable that will store the output files
  file_list <- NULL
  
  # Initialize the dataframe variable in the case that we extract values for multiple dates
  combined_df <- NULL
  
  # Initialize variable for keeping image dates 
  keep_dates <- FALSE
  
  # Preserve image dates if we are not temporally compositing images
  if (is.null(composite_method)) {
    keep_dates <- TRUE
  }
  
  for (i in 1:length(image)) {
    if (image[i] == "") next
    
    # Convert image file into raster
    DN_raster <- terra::rast(image[i])
    
    # Convert to SR and calculate indices if applicable
    SR_image <- add_indices(raster=DN_raster, file_name=file_path, count=names(image[i]), indices_df=df_indices)
    
    # Initialize data frame for extracting pixel values
    vals_df <- NULL
    
    if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
      # Convert to raster and extract pixel values
      SR_raster <- terra::rast(SR_image)
      vals_df <- extract_vals_df(SR_raster, aoi_search, use_uniqueID, 
                                 file_name=image[i], add_dates = keep_dates)
      
      if (length(image)==1) {
        # Write to fst
        file_nm <- paste0(file_path, ".fst")
        fst::write_fst(vals_df, path=file_nm)
        message("File saved to: ", file_nm)
        file_list <- file_nm
      }
      
      if (length(image) > 1) {
        if (i==1) {
          combined_df <- vals_df
        }
        else if (nrow(combined_df) == nrow(vals_df)){
          combined_df <- cbind(combined_df, vals_df[-1])
        }
      }
      
      if (i==length(image) && !is.null(combined_df)) {
        file_nm <- paste0(file_path, ".fst")
        fst::write_fst(combined_df, path=file_nm)
        message("File saved to: ", file_nm)
        file_list <- file_nm
      }
    }
    
    if ("POLYGON" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
      message("File saved to: ", SR_image)
      file_list <- append(file_list,SR_image)
    }
    
  }
  
  return(file_list)
}

s1_process <- function(aoi_layer, radius=NULL, start_dt, end_dt, uniqueID=NULL,
                     composite_method=c("median", "mean", "sum", "min", "max", NULL),
                     resample_method=c("bilinear", "average", "rms", "nearest", "gauss", "cube", "cubicspline", "lanczos", "average_magphase", "mode"),
                     calc_dB=FALSE, app_domains=NULL, idx_names=NULL, file_path) {
  
  # scoping for aoi_layer
  # Check for valid object type
  if (!inherits(aoi_layer, "sf")) {
    stop(deparse(substitute(aoi_layer)), "must be an sf object. \n")
  }
  
  # Only include polygon and point geometry types
  gtypes = c("POLYGON", "POINT")
  if (!any(gtypes %in% unique(as.character(sf::st_geometry_type(aoi_layer))))) {
    stop(deparse(substitute(aoi_layer)), " must be one of", paste(gtypes, collapse=", "),". \n")
  }
  
  # scoping for uniqueID
  use_uniqueID <- NULL
  use_uniqueID <- valid_col_ID(aoi=aoi_layer, column_ID=uniqueID)
  
  # scoping for start_dt & end_dt
  dates <- valid_dates(start_dt, end_dt)
  
  # scoping for composite_method.
  valid_composite <- c("sum", "mean", "median", "min", "max")
  if(!is.null(composite_method[1])){
    if (!composite_method[1] %in% valid_composite) {
      stop(deparse(substitute(composite_method[1])), " must be one of ", paste(valid_composite, collapse=", "),". \n")
    }
  }
  
  if (is.null(composite_method[1]) && ((as.Date(dates[2])-as.Date(dates[1]))) > 100) {
    warning("`composite_method` is `NULL` and date range is greater than 100 days. Several files will be created which may impact processing time and performance. \n",
            immediate. = TRUE)
  }
  
  # scoping for resample_method
  valid_resample <- c("nearest", "average", "rms", "bilinear", "gauss", "cube", "cubicspline",
                      "lanczos", "average_magphase", "mode")
  if (!resample_method[1] %in% valid_resample) {
    stop(deparse(substitute(resample_method[1])), " must be one of ", paste(valid_resample, collapse=", "),". \n")
  }
  
  # Set up aoi search parameter to use in STAC search
  aoi_search <- aoi_layer
  
  # Add a buffer to point geometries. Return new sf polygon layer.
  if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
    max_dist <- 1000
    aoi_search <- point_buffer(aoi_layer, radius, max_dist)
  }
  
  # Print a message if user give radius but sf object has polygon geometry.
  if ("POLYGON" %in%  unique(as.character(sf::st_geometry_type(aoi_layer))) && !is.null(radius)) {
    message("`radius` ignored for POLYGON geometry type. \n")
  }
  
  # Give warning message if aoi is very large
  area_limit <- 5000000000
  check_bbox_area(aoi_search, max_area = area_limit)
  
  # Initialize a variable to store indices
  df_indices <- NULL
  
  # Set the collection to Sentinel 1 RTC
  rtc_collection <- "sentinel-1-rtc"
  
  # Initialize image variables
  image_asc <- NULL
  image_des <- NULL
  
  # Try to get imagery from ascending orbit
  tryCatch({
    image_asc <- rsi::get_sentinel1_imagery(
      aoi = sf::st_buffer(aoi_search, dist=50),
      start_date = dates[1],
      end_date = dates[2],
      collection = rtc_collection,
      output_filename = tempfile(pattern="sentinel1_ascending", tmpdir=tempdir(), fileext=".tif"),
      query_function = s1_ascending_query,
      composite_function = composite_method[1],
      gdalwarp_options = c("-r", resample_method[1], "-multi", "-overwrite", "-co",
                           "COMPRESS=DEFLATE", "-co", "PREDICTOR=2", "-co", "NUM_THREADS=ALL_CPUS")
    )
  }, error = function(e) {
    message("No items were found for this combination of AOI and date range for Sentinel-1 ascending orbit. \n")
  })
  
  # Try to get imagery from descending orbit
  tryCatch({
    image_des <- rsi::get_sentinel1_imagery(
      aoi = sf::st_buffer(aoi_search, dist=50),
      start_date = dates[1],
      end_date = dates[2],
      collection = rtc_collection,
      output_filename = tempfile(pattern="sentinel1_descending", tmpdir=tempdir(), fileext=".tif"),
      query_function = s1_descending_query,
      composite_function = composite_method[1],
      gdalwarp_options = c("-r", resample_method[1], "-multi", "-overwrite", "-co",
                           "COMPRESS=DEFLATE", "-co", "PREDICTOR=2", "-co", "NUM_THREADS=ALL_CPUS")
    )
  }, error = function(e) {
    message("No items were found for this combination of AOI and date range for Sentinel-1 descending orbit. \n")
  })
  
  # Check if both images are NULL
  if(is.null(image_asc) && is.null(image_des)) {
    stop("No Sentinel-1 items were found for both ascending and descending orbit for this combination of AOI and date range. \n")
  }
  
  images_combined <- c(image_asc, image_des)
  file_list <- NULL
  
  for (i in 1:length(images_combined)) {
    
    # Set value for Sentinel-1 band names and return a dataframe of indices to use IF user provides
    # application domain or specific index names to search over. Otherwise, proceed with only radar data.
    out_raster <- terra::rast(images_combined[i])
    
    if (!is.null(app_domains) | !is.null(idx_names)) {
      s1_vals <- names(out_raster)
      df_indices <- indices_df(band_names=s1_vals, index_application=app_domains, index_names=idx_names)
    }
    
    if ((calc_dB==TRUE && is.null(df_indices)) | (calc_dB==TRUE && nrow(df_indices)==0)) {
      images_combined[i] <- s1_dB(out_raster, file_pth=file_path, counter=i)
    }
    
    if ((calc_dB==FALSE && is.null(df_indices)) | (calc_dB==FALSE && nrow(df_indices)==0)) {
      images_combined[i] <- paste0(file_path, "_", i, ".tif")
      w <- terra::writeRaster(out_raster, filename = images_combined[i],
                              datatype="FLT4S", filetype="GTiff",
                              gdal=c("COMPRESS=DEFLATE", "NUM_THREADS=ALL_CPUS", "PREDICTOR=2"),
                              tempdir=tempdir(), NAflag=NA, verbose=FALSE)
    }
    
    if ((calc_dB==TRUE && !is.null(df_indices)) && nrow(df_indices) > 0) {
      images_combined[i] <- s1_dB(out_raster)
    }
    
    # Add indices
    if (!is.null(df_indices)) {
      if (nrow(df_indices) > 0) {
        indices_image <- rsi::calculate_indices(raster = out_raster, indices = df_indices,
                                                output_filename = tempfile(pattern="sentinel1_indices", tmpdir=tempdir(), fileext=".tif"))
        
        images_combined[i] <- rsi::stack_rasters(c(images_combined[i],indices_image), output_filename=paste0(file_path,"_",i,".tif"))
      }
    }
    
    # Initialize dataframe for extracting pixel values
    vals_df <- NULL
    
    if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
      # Convert to raster and extract pixel values
      radar_raster <- terra::rast(images_combined[i])
      vals_df <- extract_vals_df(radar_raster, aoi_search, use_uniqueID)
      
      # Write to fst
      file_nm <- paste0(file_path,"_",i,".fst")
      fst::write_fst(vals_df, path=file_nm)
      message("File saved to: ", file_nm)
      file_list <- append(file_list, file_nm)
    }
    
    if ("POLYGON" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
      message("File saved to: ", images_combined[i])
      file_list <- append(file_list,images_combined[i])
    }
    
  }
  
  return(file_list)
}

dem_process <- function(aoi_layer, radius=NULL, uniqueID=NULL, file_path) {
  
  # scoping for aoi_layer
  # Check for valid object type
  if (!inherits(aoi_layer, "sf")) {
    stop(deparse(substitute(aoi_layer)), " must be an sf object. \n")
  }
  
  # Only include polygon and point geometry types
  gtypes = c("POLYGON", "POINT")
  if (!any(gtypes %in% unique(as.character(sf::st_geometry_type(aoi_layer))))) {
    stop(deparse(substitute(aoi_layer)), " must be one of ", paste(gtypes, collapse=", "),". \n")
  }
  
  # scoping for uniqueID
  use_uniqueID <- NULL
  use_uniqueID <- valid_col_ID(aoi=aoi_layer, column_ID=uniqueID)
  
  
  # Set up aoi search parameter to use in STAC search
  aoi_search <- aoi_layer
  
  # Add a buffer to point geometries. Return new sf polygon layer.
  if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
    max_dist <- 1000
    aoi_search <- point_buffer(aoi_layer, radius, max_dist)
  }
  
  if ("POLYGON" %in%  unique(as.character(sf::st_geometry_type(aoi_layer))) && !is.null(radius)) {
    message("`radius` ignored for POLYGON geometry type. \n")
  }
  
  # Give warning message if aoi is very large
  area_limit <- 5000000000
  check_bbox_area(aoi_search, max_area = area_limit)
  
  # Initialize image variable
  image <- NULL
  
  # Set default image file name: 
  image_file <- tempfile(pattern="dem_image", tmpdir=tempdir(), fileext=".tif")
  
  # Override image file name for polygon geometry
  if ("POLYGON" %in%  unique(as.character(sf::st_geometry_type(aoi_layer)))) {
    image_file <- paste0(file_path, ".tif")
  }
  
  tryCatch({
    # Retrieve raster data
    # Add a small buffer to aoi layer to get all pixels within the aoi boundary after reprojection (see rsi::get_stac_data help text)
    image <- rsi::get_dem(
      aoi = sf::st_buffer(aoi_search, dist=50),
      output_filename = image_file
    )
  }, error = function(e) {
    message("No DEM items were found.")
  })
  
  if (is.null(image)) {
    stop("No DEM items were found.")
  }
  
  file_list <- NULL
  
  for (i in 1:length(image)) {
    # Convert image file into raster
    DEM_raster <- terra::rast(image[i])
    
    vals_df <- NULL
    
    if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
      # Convert to raster and extract pixel values
      vals_df <- extract_vals_df(DEM_raster, aoi_search, use_uniqueID)
      
      # Write to fst
      file_nm <- paste0(file_path,"_",i,".fst")
      fst::write_fst(vals_df, path=file_nm)
      message("File saved to: ", file_nm)
      file_list <- append(file_list, file_nm)
    }
    
    if ("POLYGON" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
      message("File saved to: ", image[i])
      file_list <- append(file_list,image[i])
    }
    
  }
  
  return(file_list)
}


esa_worldcover_process <- function(aoi_layer, radius=NULL, start_dt, end_dt, uniqueID=NULL, file_path) {
  
  # scoping for aoi_layer
  # Check for valid object type
  if (!inherits(aoi_layer, "sf")) {
    stop(deparse(substitute(aoi_layer)), " must be an sf object. \n")
  }
  
  # Only include polygon and point geometry types
  gtypes = c("POLYGON", "POINT")
  if (!any(gtypes %in% unique(as.character(sf::st_geometry_type(aoi_layer))))) {
    stop(deparse(substitute(aoi_layer)), " must be one of ", paste(gtypes, collapse=", "),". \n")
  }
  
  # scoping for uniqueID
  use_uniqueID <- NULL
  use_uniqueID <- valid_col_ID(aoi=aoi_layer, column_ID=uniqueID)
  
  # scoping for start_dt & end_dt
  dates <- valid_dates(start_dt, end_dt)
  
  # Set up aoi search parameter to use in STAC search
  aoi_search <- aoi_layer
  
  # Add a buffer to point geometries. Return new sf polygon layer.
  if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
    max_dist <- 1000
    aoi_search <- point_buffer(aoi_layer, radius, max_dist)
  }
  
  if ("POLYGON" %in%  unique(as.character(sf::st_geometry_type(aoi_layer))) && !is.null(radius)) {
    message("`radius` ignored for POLYGON geometry type. \n")
  }
  
  # Give warning message if aoi is very large
  area_limit <- 5000000000
  check_bbox_area(aoi_search, max_area = area_limit)
  
  # Initialize image variable
  image <- NULL
  
  # Set up sign function to sign the item URLs for PC
  esa_coll <- "esa-worldcover"
  esa_assets <- "map"
  pc_stac_source <- "https://planetarycomputer.microsoft.com/api/stac/v1/"
  
  for(i in seq(2020, 2021, 1)) {
    # Check if query dates contain the ESA Worldcover Temporal extent
    search_dates <- check_date_in_range(dates[1], dates[2], as.character(i))
    if (is.null(search_dates)) {
      next
    }
    
    # Set default image file name: 
    image_file <- tempfile(pattern=paste0("esa_img", as.character(i)), tmpdir=tempdir(), fileext=".tif")
    
    # Override image file name for polygon geometry
    if ("POLYGON" %in%  unique(as.character(sf::st_geometry_type(aoi_layer)))) {
      image_file <- paste0(file_path, "_", as.character(i), ".tif")
    }
    
    tryCatch({
      # Retrieve raster data
      # Add a small buffer to aoi layer to get all pixels within the aoi boundary after reprojection (see rsi::get_stac_data help text)
      image <- c(image, rsi::get_stac_data(
        aoi = sf::st_buffer(aoi_search, dist=50),
        start_date = as.character(search_dates[1]),
        end_date = as.character(search_dates[2]),
        pixel_x_size = 10,
        pixel_y_size = 10,
        asset_names = esa_assets,
        stac_source = pc_stac_source,
        collection = esa_coll,
        sign_function = rsi::sign_planetary_computer,
        rescale_bands=FALSE,
        download_function = rsi::rsi_download_rasters,
        composite_function = "max",
        gdalwarp_options = c("-r", "near", "-multi", "-overwrite", "-co",
                             "COMPRESS=DEFLATE", "-co", "PREDICTOR=2", "-co", "NUM_THREADS=ALL_CPUS"),
        output_filename = image_file
      ))
    }, error = function(e) {
      message("No ESA WorldCover items were found for ", as.character(i))
    })
    
  }
  
  if (is.null(image)) {
    stop("No ESA WorldCover items were found.")
  }
  
  file_list <- NULL
  
  for (i in 1:length(image)) {
    # Convert image file into raster
    ESA_raster <- terra::rast(image[i])
    
    vals_df <- NULL
    
    if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
      # Convert to raster and extract pixel values
      vals_df <- extract_vals_df(ESA_raster, aoi_search, use_uniqueID)
      
      # Write to fst
      file_nm <- paste0(file_path,"_",i,".fst")
      fst::write_fst(vals_df, path=file_nm)
      message("File saved to: ", file_nm)
      file_list <- append(file_list, file_nm)
    }
    
    if ("POLYGON" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
      message("File saved to: ", image[i])
      file_list <- append(file_list,image[i])
    }
    
  }
  
  return(file_list)
}

alos_fnf_process <- function(aoi_layer, radius=NULL, start_dt, end_dt, uniqueID=NULL, file_path) {
  
  # scoping for aoi_layer
  # Check for valid object type
  if (!inherits(aoi_layer, "sf")) {
    stop(deparse(substitute(aoi_layer)), " must be an sf object. \n")
  }
  
  # Only include polygon and point geometry types
  gtypes = c("POLYGON", "POINT")
  if (!any(gtypes %in% unique(as.character(sf::st_geometry_type(aoi_layer))))) {
    stop(deparse(substitute(aoi_layer)), " must be one of ", paste(gtypes, collapse=", "),". \n")
  }
  
  # scoping for uniqueID
  use_uniqueID <- NULL
  use_uniqueID <- valid_col_ID(aoi=aoi_layer, column_ID=uniqueID)
  
  # scoping for start_dt & end_dt
  dates <- valid_dates(start_dt, end_dt)
  
  # Set up aoi search parameter to use in STAC search
  aoi_search <- aoi_layer
  
  # Add a buffer to point geometries. Return new sf polygon layer.
  if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
    max_dist <- 1000
    aoi_search <- point_buffer(aoi_layer, radius, max_dist)
  }
  
  if ("POLYGON" %in%  unique(as.character(sf::st_geometry_type(aoi_layer))) && !is.null(radius)) {
    message("`radius` ignored for POLYGON geometry type. \n")
  }
  
  # Give warning message if aoi is very large
  area_limit <- 5000000000
  check_bbox_area(aoi_search, max_area = area_limit)
  
  # Initialize image variable
  image <- NULL
  
  # Set up sign function to sign the item URLs for PC
  alos_coll <- "alos-fnf-mosaic"
  alos_assets <- "C"
  pc_stac_source <- "https://planetarycomputer.microsoft.com/api/stac/v1/"
  
  for (i in seq(2015, 2020, 1)) {
    # Check if query dates contain the ALOS FNF temporal extent
    search_dates <- check_date_in_range(dates[1], dates[2], as.character(i))
    if (is.null(search_dates)) {
      next
    }
    
    # Set default image file name
    image_file <- tempfile(pattern=paste0("alos_img", as.character(i)), tmpdir=tempdir(), fileext=".tif")
    
    # Override image file name for polygon geometry
    if ("POLYGON" %in%  unique(as.character(sf::st_geometry_type(aoi_layer)))) {
      image_file <- paste0(file_path, "_", as.character(i), ".tif")
    }
    
    tryCatch({
      # Retrieve raster data
      # Add a small buffer to aoi layer to get all pixels within the aoi boundary after reprojection (see rsi::get_stac_data help text)
      image <- c(image, rsi::get_stac_data(
        aoi = sf::st_buffer(aoi_search, dist=50),
        start_date = as.character(search_dates[1]),
        end_date = as.character(search_dates[2]),
        pixel_x_size = 25,
        pixel_y_size = 25,
        asset_names = alos_assets,
        stac_source = pc_stac_source,
        collection = alos_coll,
        sign_function = rsi::sign_planetary_computer,
        rescale_bands=FALSE,
        download_function = rsi::rsi_download_rasters,
        composite_function = "max",
        gdalwarp_options = c("-r", "near", "-multi", "-overwrite", "-co",
                             "COMPRESS=DEFLATE", "-co", "PREDICTOR=2", "-co", "NUM_THREADS=ALL_CPUS"),
        output_filename = image_file
      ))
    }, error = function(e) {
      message("No ALOS FNF items were found for ", as.character(i))
    })
  }
  
  if (is.null(image)) {
    stop("No ALOS FNF items were found.")
  }
  
  file_list <- NULL
  
  for (i in 1:length(image)) {
    # Convert image file into raster
    FNF_raster <- terra::rast(image[i])
    
    vals_df <- NULL
    
    if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
      # Convert to raster and extract pixel values
      vals_df <- extract_vals_df(FNF_raster, aoi_search, use_uniqueID)
      
      # Write to fst
      file_nm <- paste0(file_path,"_",i,".fst")
      fst::write_fst(vals_df, path=file_nm)
      message("File saved to: ", file_nm)
      file_list <- append(file_list, file_nm)
    }
    
    if ("POLYGON" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
      message("File saved to: ", image[i])
      file_list <- append(file_list,image[i])
    }
    
  }
  
  return(file_list)
}

soil_process <- function(aoi_layer, radius=NULL, uniqueID=NULL, soil_layers, file_path) {
  
  # scoping for aoi_layer
  # Check for valid object type
  if (!inherits(aoi_layer, "sf")) {
    stop(deparse(substitute(aoi_layer)), " must be an sf object. \n")
  }
  
  # Only include polygon and point geometry types
  gtypes = c("POLYGON", "POINT")
  if (!any(gtypes %in% unique(as.character(sf::st_geometry_type(aoi_layer))))) {
    stop(deparse(substitute(aoi_layer)), " must be one of ", paste(gtypes, collapse=", "),". \n")
  }
  
  # scoping for uniqueID
  use_uniqueID <- NULL
  use_uniqueID <- valid_col_ID(aoi=aoi_layer, column_ID=uniqueID)
  
  # Set up aoi search parameter to use in API call
  aoi_search <- aoi_layer
  
  # Add a buffer to point geometries. Return new sf polygon layer.
  if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
    max_dist <- 1000
    aoi_search <- point_buffer(aoi_layer, radius, max_dist)
  }
  
  if ("POLYGON" %in%  unique(as.character(sf::st_geometry_type(aoi_layer))) && !is.null(radius)) {
    message("`radius` ignored for POLYGON geometry type. \n")
  }
  
  # Get bounding box in Homolosine projection
  # As the ISRIC layers use the Homolosine projection, we need to reproject the SF object:
  igh <- "+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs"
  bbox_igh <- get_bbox_in_crs(sf_aoi=aoi_search, radius_src=radius, crs_val=igh)
  
  # Get valid soil layers to use in GDAL translate utility
  valid_soil_srch <- validate_soil_layers(soil_layers)
  
  # Get list of GeoTIFF soil files: 
  soil_files <- get_soil_data(valid_soil_layers=valid_soil_srch, bounding_box=bbox_igh)
  
  # Convert output GeoTIFFs into terra SpatRaster object
  soil_rast <- terra::rast(soil_files)
  
  # Reproject the raster to match the original sf object projection
  soil_rast <- rast_reproject(raster_layer=soil_rast, sf_layer=aoi_layer)
  
  # Extract values into dataframe for point geometry. 
  if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
    #Intialize dataframe variable: 
    vals_df <- NULL
    
    # Convert to raster and extract pixel values
    vals_df <- extract_vals_df(soil_rast, aoi_search, use_uniqueID)
    
    # Write to fst
    file_nm <- paste0(file_path, ".fst")
    fst::write_fst(vals_df, path=file_nm)
    message("File saved to: ", file_nm)
  }
  
  # Save as .tif for polygon geometry. 
  if ("POLYGON" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
    file_nm <- paste0(file_path, ".tif")
    writeRaster(soil_rast, filename=file_nm)
    message("File saved to: ", file_nm)
  }
  
  return(file_nm)
  
}

gedi_process <- function(aoi_layer, radius=NULL, uniqueID=NULL, gedi_product='GEDI03', netrc_path, file_path) {
  # scoping for aoi_layer
  # Check for valid object type
  if (!inherits(aoi_layer, "sf")) {
    stop(deparse(substitute(aoi_layer)), " must be an sf object. \n")
  }
  
  # Only include polygon and point geometry types
  gtypes = c("POLYGON", "POINT")
  if (!any(gtypes %in% unique(as.character(sf::st_geometry_type(aoi_layer))))) {
    stop(deparse(substitute(aoi_layer)), " must be one of ", paste(gtypes, collapse=", "),". \n")
  }
  
  # scoping for uniqueID
  use_uniqueID <- NULL
  use_uniqueID <- valid_col_ID(aoi=aoi_layer, column_ID=uniqueID)
  
  # scoping for GEDI product. Currently just GEDI03, but might incorporate more. 
  valid_product <- "GEDI03"
  if(!is.null(gedi_product[1])){
    if (!gedi_product[1] %in% valid_product) {
      stop(deparse(substitute(gedi_product[1])), " must be one of ", paste(valid_product, collapse=", "), ". \n")
    }
  }
  
  # Set up aoi search parameter to use in STAC search
  aoi_search <- aoi_layer
  
  # Add a buffer to point geometries. Return new sf polygon layer.
  if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
    max_dist <- 1000
    aoi_search <- point_buffer(aoi_layer, radius, max_dist)
  }
  
  if ("POLYGON" %in%  unique(as.character(sf::st_geometry_type(aoi_layer))) && !is.null(radius)) {
    message("`radius` ignored for POLYGON geometry type. \n")
  }
  
  # Give warning message if aoi is very large
  area_limit <- 5000000000
  check_bbox_area(aoi_search, max_area = area_limit)
  
  # Run the gedi_finder function
  granules <- gedi_finder(gedi_product)
  
  # Subset output to last 5 granules (most recent datasets)
  granules_5 <- gedi_get_latest_data(granules)
  
  # Download the data: 
  gedi_data <- gedi_download(granules_5, netrc_path)
  
  # Transform .tif files into 1 raster
  gedi_rast <- gedi_raster(gedi_data)
  
  # Reproject and crop
  gedi_crop <- gedi_reproject(gedi_rast, aoi_search)
  
  file_nm <- NULL
  
  if ("POINT" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
    # Extract pixel values
    vals_df <- extract_vals_df(gedi_crop, aoi_search, use_uniqueID)
    
    # Write to fst
    file_nm <- paste0(file_path, ".fst")
    fst::write_fst(vals_df, path=file_nm)
    message("File saved to: ", file_nm)
  }
  
  if ("POLYGON" %in% unique(as.character(sf::st_geometry_type(aoi_layer)))) {
    # Write to GeoTIFF
    file_nm <- paste0(file_path,".tif")
    terra::writeRaster(gedi_crop, filename = file_nm,
                       datatype="FLT4S", filetype="GTiff",
                       gdal=c("COMPRESS=DEFLATE", "NUM_THREADS=ALL_CPUS", "PREDICTOR=2"),
                       tempdir=tempdir(), NAflag=NA, verbose=FALSE)
    
    message("File saved to: ", file_nm)
  }
  
  on.exit(file.remove(unlist(gedi_data)), add=TRUE)
  return(file_nm)
  
}