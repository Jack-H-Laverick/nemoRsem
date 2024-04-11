
#### Utilities to help with manipulating data ####


#' Categorise files for processing with `NEMO_ERSEM`
#'
#' This function takes a directory, recursively lists the file names, and extracts meaningful meta-data
#'
#' @param dir A path string to the directory containing NEMO-ERSEM netcdf files.
#' @param recursive T/F whether to process files in all sub-directories.
#' @return A dataframe that can be passed to `NEMO_ERSEM` using pmap.
#' @family NEMO-ERSEM spatial tools
#' @export
categorise_files <- function(dir, recursive = TRUE) {


  ersem_files <- list.files(dir, recursive = recursive, full.names = TRUE, pattern = ".nc") %>%  # Retunr netcdf files in the directory
    as.data.frame() %>%                                                                     # Turn the vector into a dataframe
    tidyr::separate(".", into = c("Path", "File"), sep = "STRATH") %>%                      # Split file names from path for ease
    dplyr::mutate(File = paste0("STRATH", File)) %>%                                        # Reintroduce separation character string
    dplyr::mutate(date = stringr::str_sub(File, start = -9, end = -4),                      # Pull time information
           Month = stringr::str_sub(date, start = 5, end = 6),
           Year = stringr::str_sub(date, start = 1, end = 4)) %>%
    dplyr::mutate(String = dplyr::case_when(stringr::str_detect(File, "thetao_con") ~ "thetao_con-Temperature", # Categorise the file types
                            stringr::str_detect(File, "so_abs") ~ "so_abs-Salinity",
                            stringr::str_detect(File, "uo") ~ "uo-Zonal currents",
                            stringr::str_detect(File, "vo") ~ "vo-Meridional currents",
                            stringr::str_detect(File, "wo") ~ "wo-Vertical velocity",
                            stringr::str_detect(File, "difvho") ~ "difvho-Vertical diffusivity",
                            stringr::str_detect(File, "R1") ~ "R1-Dissolved organic nitrogen",
                            stringr::str_detect(File, "O2") ~ "O2-Oxygen",
                            stringr::str_detect(File, "N4") ~ "N4-Ammonium",
                            stringr::str_detect(File, "N3") ~ "N3-Nitrate",
                            stringr::str_detect(File, "RP") ~ "RP-Detritus (Particulate organic nitrogen)",
                            stringr::str_detect(File, "B1") ~ "B1-Bacterial nitrogen",
                            stringr::str_detect(File, "P1") ~ "P1-Diatom nitrogen",
                            stringr::str_detect(File, "P234_n") ~ "P234_n-Other phytoplankton nitrogen")) %>%
    tidyr::separate("String", into = c("Type", "Name"), sep = "-") %>%                     # cheat way to succinctly get both a "type" and named variable column
    dplyr::select(Path, File, date, Year, Month, Type, Name)                               # Limit to columns of use

}

#' Calculate Water Layer Thicknesses Within a Vector
#'
#' This function calculates the thickness of the water layer around each point in a vector. These thicknesses
#' are needed to calculate weighted averages across a depth window.
#'
#' The function calculates the midpoints between values in a vector and an optional set of boundary depths.
#' The differences between these midpoints are calculated to return the thickness of a depth layer centered
#' on each point in the original vector. In the absence of boundary depths, the maximum and minimum depths are used.
#'
#' If a depth falls outside the target window it gets a thickness of 0. If all depths fall outside the window
#' the function returns an all 0 vector. If this is passed to `weighted.mean` the result will be NaN.
#'
#'
#' @param depths A numeric vector of increasing depths.
#' @param min_depth The shallowest depth in the depth window. (defaults to minimum depth in the vector)
#' @param max_depth The deepest depth in the depth window. (defaults to minimum depth in the vector)
#' @return A vector of water layer thicknesses to match the length of the original vector.
#' @family NEMO-ERSEM spatial tools
#' @examples
#' # Get a vector of depths
#' depths <- seq(0, 100, by = 10)
#'
#' # Water layer thickness within the vector
#' calculate_depth_share(depths)
#'
#' # Water layer thickness using limits of a depth window
#' calculate_depth_share(depths, min_depth = 25, max_depth = 75)
#'
#' # Special case when the depth vector falls outside the target depth window
#' calculate_depth_share(depths, min_depth = 400, max_depth = 600)
#' @export
calculate_depth_share <- function(depths, min_depth = min(depths), max_depth = max(depths)) {

  contained <- depths <= max_depth & depths >= min_depth       # Which depths are within our window?

  if (all(!contained)) {

    depths <- rep(0, length(depths))                           # If none of the vector entries are within the depth window, return all 0s

  } else {

    weights <- diff(c(min_depth, RcppRoll::roll_mean(depths[contained], n = 2), max_depth)) # Calculate the midpoints and the difference between them (width either side of the original point)

    depths[contained] <- weights                               # Assign the weights to positions in the vector in our window
    depths[!contained] <- 0                                    # Positions outside the vector get 0 weight
  }
  return(depths)

}

#' Calculate the Weights for Linear Interpolation Between Two Depths
#'
#' This function calculates the distance between a shallower and deeper depth to a target depth. These values
#' are then swapped so they can be used in a weighted average to linearly interpolate to the target depth.
#'
#' @param depths A numeric vector of a shallower and deeper depth.
#' @param target An intermediate depth we would like to interpolate to.
#' @return A vector of weights to match the vector of depths.
#' @family NEMO-ERSEM spatial tools
#' @examples
#' calculate_proximity_weight(c(0, 100), target = 30)
#'
#' calculate_proximity_weight(c(0, 100), target = 50)
#'
#' calculate_proximity_weight(c(0, 100), target = 80)
#' @export
calculate_proximity_weight <- function(depths, target) rev(diff(c(depths[1], target, depths[2])))

#' Get Latitudes, Longitudes, & Depths From NEMO-ERSEM Model Outputs
#'
#' This function gets the latitudes, longitudes, and depths which define the spatial location of points in an array of NEMO-ERSEM outputs.
#'
#' Each variable of interest in the netcdf file is imported, and then collected into a list.
#'
#' @param file The full name of a netcdf file.
#' @param depthvar The name of the depth variable netcdf file (when present).
#' @param depthvar The name of the depth dimension in a netcdf file (when extracting values from dimension metadata).
#' @return A list of three elements:
#' \itemize{
#'  \item{\emph{nc_lat -}}{ A matrix of latitudes which maps onto the first and second dimension of a NEMO-ERSEM array.}
#'  \item{\emph{nc_lon -}}{ A matrix of longitudes which maps onto the first and second dimension of a NEMO-ERSEM array.}
#'  \item{\emph{nc_depth -}}{ A vector of depths which match the third dimension of a NEMO-ERSEM array.}
#'  }
#' @family NEMO-ERSEM variable extractors
#' @export
get_spatial <- function(file, depthvar = NULL, depthdim = NULL) {

  nc_raw <- ncdf4::nc_open(file)                       # Open up a netcdf file to see it's raw contents (var names)

  nc_lat <- ncdf4::ncvar_get(nc_raw, "nav_lat")        # Extract a matrix of all the latitudes
  nc_lon <- ncdf4::ncvar_get(nc_raw, "nav_lon")        # Extract a matrix of all the longitudes

  if(is.null(depthvar) == FALSE) nc_depth <- ncdf4::ncvar_get(nc_raw, depthvar)
  if(is.null(depthdim) == FALSE) eval(parse(text = paste0("nc_depth <- nc_raw$dim$",depthdim,"$vals")))

  ncdf4::nc_close(nc_raw)                              # You must close an open netcdf file when finished to avoid data loss

  all <- list("nc_lat" = nc_lat, "nc_lon" = nc_lon, "nc_depth" = nc_depth)
  return(all)
}

#' Convert XY coordinates to a value index for a matrix
#'
#' This function converts a double index for a matrix (xy) into a single index for values. This allows vectorised subsetting
#' of incomplete rows and columns.
#'
#' @param x vector of row numbers for target values.
#' @param y vector of column numbers for target values.
#' @param nrow The number of rows in the matrix
#' @return A vector of indices for values in a matrix.
#' @family NEMO-ERSEM spatial tools
#' @export
xyindex_to_nindex <- function(x, y, nrow) {x + ((y-1)*nrow)}

#' Calculate the Domain Area per Grid Point
#'
#' This function takes an array of a variable, and an array of water thicknesses to perform a weighted average across depth. The depth
#' window to be averaged can be specified, so this function can be used to create both shallow and deep layers (or more for that matter).
#'
#' @param points A Simple Feature object og the grid points within the model domain.
#' @param area A Simple Feature object containing the model domain.
#' @return the `points` object is returned, but instead of points, the geometry column now contains polygons representing the area closest to each point. A column for the size of this area is also gained.
#' @family NEMO-ERSEM spatial tools
#' @export
voronoi_grid <- function(points, area) {

  result <- purrr::map(1:nrow(area), ~{                            # For each polygon in area
    voronoi <- points %>%                                          # Take the grid points
      sf::st_geometry() %>%                                        # To get sfc from sf
      sf::st_union() %>%                                           # To get a sfc of MULTIPOINT type
      sf::st_voronoi(envelope = sf::st_geometry(area[.x,])) %>%    # Voronoi polygon for the area
      sf::st_collection_extract(type = "POLYGON") %>%              # A list of polygons
      sf::st_sf() %>%                                              # From list to sf object
      sf::st_join(points) %>%                                      # put names back
      sf::st_intersection(area[.x,]) %>%                           # Cut to shape of target area
      dplyr::mutate(Cell_area = units::drop_units(sf::st_area(.))) # Area of each polygon
  }) %>%
    dplyr::bind_rows() %>%                                         # Combine the results from each area
    sf::st_sf(geomc = .$geometry, crs = 4326)                      # Reinstate attributes of the geometry column

}

#' Summarise NEMO-ERSEM Output into Time Series Along Transects
#'
#' This function averages NEMO-ERSEM monthly summaries into time series for the target boundaries of StrathE2E.
#'
#' The function subsets the NEMO-ERSEM grid according to the transects object provided. Water exchanges between
#' model compartments are totaled. The boundary conditions of the model domain for variables needed by
#' StrathE2E are summarised as a flow-weighted mean, applying the flow rate at each transect.
#'
#' @param saved A dataframe containing a summarised month from NEMO-ERSEM model outputs.
#' @param transects A dataframe containing the labelled transects along the model domain boundaries.
#' @param vars A character vector containing the column names to be summarised for boundary conditions. Defaults to the targets for StrathE2E
#' @return A list containing two summaries.
#' \itemize{
#'  \item{Element 1}{A dataframe containing the total water exchanged between model compartments.}
#'  \item{Element 2}{A dataframe containing the flow-weighted boundary conditions around the model domain.}}
#' @family NEMO-MEDUSA averages
#' @export
NE_boundary_summary <- function(saved, transects, vars = c("NO3", "NH4", "Chlorophyll", "Temperature")) {

  Data <- readRDS(saved) %>%                                                  # Import a NM summary object
    dplyr::select(-c(Shore, weights))                                         # Drop duplicated columns which vonflict
  data.table::setDT(Data, key = c("x", "y", "slab_layer"))                    # Convert to a data.table keyed spatially for quick summaries.

  join <- Data[transects] %>%                                                 # 0.5% of the transects don't catch data (NA), it's because the deep offshore layer doesn't have a buffer on the shore side.
    dplyr::mutate(Flow = ifelse(current == "Zonal", Zonal, Meridional)) %>%   # Grab current perpendicular to the transect
    dplyr::mutate(Flow = ifelse(Flip == T, -1*Flow, Flow)) %>%                # Correct flow so that + always goes IN to a model box
    dplyr::mutate(Flow = Flow * weights,                                      # Weight by transect area to get the volume of water
                  Direction = ifelse(Flow > 0, "In", "Out"))                         # Label direction based on flow rate

  ## Summarise water exchanges

  water <- join[, .(Flow = sum(Flow, na.rm = T)),                             # Tally up water movements
                by = c("Shore", "slab_layer", "Direction",                    # By exchanges we want to keep track of
                       "Neighbour", "Month", "Year")] %>%
    tidyr::drop_na()                                                          # The NA transects introduce a dead group, remove.

  ## Summarise boundary conditions

  #*# How do we weight by flow? some are negative
  boundary <- join[perimeter == T & Direction == "In",                        # For transects which bound the perimeter of the model domain
                   lapply(.SD, weighted.mean, w = Flow, na.rm = T),           # Flow-weighted mean of target variables
                   by = c("Shore", "slab_layer", "Neighbour",                 # By groups we want to keep track of
                          "Month", "Year"),
                   .SDcols = vars] %>%                                        # Specify target variables
    tidyr::drop_na() %>%                                                      # The NA transects introduce a dead group, remove.
    dplyr::mutate(Date = as.Date(paste(15, Month, Year, sep = "/"), format = "%d/%m/%Y"),
                  Compartment = paste(Shore, slab_layer)) %>%
    tidyr::pivot_longer(eval(vars), names_to = "Variable", values_to = "Measured") # reshape

  result <- list(Flows = water, Boundarys = boundary)                         # Combine both summaries so they can be returned together

  return(result)
}

#' Summarise NEMO-ERSEM Output into Time Series Within Model Compartments
#' This function averages NEMO-ERSEM monthly summaries into time series for each model compartment.
#'
#' The function groups by model compartment (Depth and Shore zone) and time step (Month and Year).
#' The mean for every target variable is calculated within these groups.
#'
#' The ice-threshold parameter controls how high the concentration of ice must be at a pixel to count as ice-affected.
#' The default is set to 0, so any pixels with ice are counted, but if you decide that ice-concentrations below 0.05
#' are meaningless for your purposes, you can set the threshold to have these pixels classed as ice free instead. This will
#' also drop the pixels from calculations of Ice_concentration, Ice_Thickness, and Snow_Thickness.
#'
#' @param saved A dataframe containing a summarised month from NEMO-ERSEM model outputs.
#' @param ice_threshold A value between 0 and 1 defining ice free pixels.
#' @param ice A TRUE FALSE switch for whether the ice_mod files were included when extracting from NEMO-ERSEM.
#' @return A dataframe containing a mean monthly time series of all target variables in NEMO-ERSEM outputs.
#' @family NEMO-MEDUSA averages
#' @export
NE_volume_summary <- function(saved, ice_threshold = 0, ice = FALSE) {

  if (ice == TRUE) {

    Groups <- readRDS(file = saved) %>%                                          # Read in wide format data file
      tidyr::drop_na(Year, Shore) %>%                                            # Drop points outside of the polygons
      dplyr::group_by(Shore, Year, Month, slab_layer) %>%
      dplyr::mutate(Ice_pres = ifelse(Ice_conc < ice_threshold, 0, Ice_pres))    # Specify how much ice actually matters when labelling something as ice-affected

    Ice <- dplyr::filter(Groups, Ice_pres > 0) %>%                               # Remove ice free pixels before averaging
      dplyr::summarise(Ice_Thickness_avg = mean(Ice_Thickness, na.rm = TRUE),    # Get monthly mean sea ice thickness
                       Snow_Thickness_avg = mean(Snow_Thickness, na.rm = TRUE),         # Get monthly mean snow thickness
                       Ice_conc_avg = mean(Ice_conc, na.rm = TRUE))                     # Get monthly mean sea ice concentration

    Averaged <- Groups %>%
      dplyr::summarise(Salinity_avg = stats::weighted.mean(Salinity, weights, na.rm = TRUE), # Get monthly mean salinity
                       Temperature_avg = stats::weighted.mean(Temperature, weights, na.rm = TRUE),
                       DIN_avg = stats::weighted.mean(DIN, weights, na.rm = TRUE),
                       Detritus_avg = stats::weighted.mean(Detritus, weights, na.rm = TRUE),
                       Phytoplankton_avg = stats::weighted.mean(Phytoplankton, weights, na.rm = TRUE),
                       Ice_pres = mean(Ice_pres, na.rm = TRUE),                         # Proportion of pixels covered by ice
                       Meridional_avg = stats::weighted.mean(Meridional, weights, na.rm = TRUE),
                       Zonal_avg = stats::weighted.mean(Zonal, weights, na.rm = TRUE)) %>%
      dplyr::left_join(Ice) %>%                                                # Add in ice and snow thicknesses
      dplyr::ungroup()
  }

  if (ice == FALSE) {

    Groups <- readRDS(file = saved) %>%                                          # Read in wide format data file
      tidyr::drop_na(Year, Shore) %>%                                            # Drop points outside of the polygons
      dplyr::group_by(Shore, Year, Month, slab_layer)

    Averaged <- Groups %>%
      # dplyr::summarise(Salinity_avg = stats::weighted.mean(Salinity, weights, na.rm = TRUE), # Get monthly mean salinity
      #                  Temperature_avg = stats::weighted.mean(Temperature, weights, na.rm = TRUE),
      #                  NO3_avg = stats::weighted.mean(NO3, weights, na.rm = TRUE),
      #                  NH4_avg = stats::weighted.mean(NH4, weights, na.rm = TRUE),
      #                  #Detritus_avg = stats::weighted.mean(Detritus, weights, na.rm = TRUE),
      #                  Chlorophyll_avg = stats::weighted.mean(Chlorophyll, weights, na.rm = TRUE),
      #                  Meridional_avg = stats::weighted.mean(Meridional, weights, na.rm = TRUE),
      #                  Zonal_avg = stats::weighted.mean(Zonal, weights, na.rm = TRUE)) %>%
      dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(.x, na.rm = TRUE), .names = "{.col}_avg"))
      dplyr::ungroup()
  }


  return(Averaged) }

#' Summarise NEMO-ERSEM Output into Decadal Grids
#' This function averages cleaned NEMO-ERSEM monthly summaries into decadal grids.
#'
#' The function groups by all spatial variables (Longitude, Latitude, Depth, and Shore zone), and by decade and month.
#' The mean for every other variable is calculated within these groups.
#'
#' @param saved A dataframe containing a summarised month from NEMO-ERSEM model outputs. It must contain the columns:
#' Longitude, Latitude, Decade, Month, Shore, and Depth.
#' @param dt Switch for using either data.table or dplyr methods (TRUE/FALSE respectively)
#' @return A dataframe containing a summarised decade of spatialy resolved NEMO-ERSEM outputs.
#' @family NEMO-ERSEM averages
#' @export
NE_decadal_summary <- function(decade, dt) {

  if(dt == TRUE){                                                               # Run data.table method
    # data.table::setDT(decade)                                                 # set as a data.table, not needed if decade is already a data.table
    Averaged <- decade[, lapply(.SD, mean, na.rm = TRUE),                       # Average data columns which aren't groups
                       by = c("longitude", "latitude", "Decade", "Month", "Shore", "slab_layer")] # Group by pixel and decade
  } else{                                                                       # Run dplyr method
    Averaged <- decade %>%
      dplyr::group_by(longitude, latitude, Decade, Month, Shore, slab_layer) %>%# Group by pixel and decade
      dplyr::summarise_all(mean, na.rm = TRUE) %>%                              # Average data columns
      dplyr::ungroup()                                                          # Ungroup
  }
  return(Averaged)
}
