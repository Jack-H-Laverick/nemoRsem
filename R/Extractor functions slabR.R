
#### slabR ####

#' Extract summaries from NEMO-ERSEM arrays using slabR
#'
#' These functions read in target variables from NEMO-ERSEM model outputs and return weighted averages
#' according to a summary scheme.
#'
#' @details Each variable of interest in a netcdf file is imported, only reading within an x/y window specified
#' with `start` and `count`. The values are then passed to `array_w_mean()` to summarise according to a scheme.
#'
#' Each model output file type contains a different set of variables, extracted by the relevant function variant:
#'
#' | File type &nbsp; &nbsp; &nbsp; &nbsp;| Variables    |
#' |--------------|-------------|
#' | grid_T_      | Salinity, temperature.|
#' | grid_U_      | Zonal currents.|
#' | grid_V_      | Meridional currents.|
#' | grid_W_      | Vertical velocity.|
#' | diad_        | Surface irradiance.|
#' | ptrc_T_      | N03, NH4, POC, phytoplankton Chl-a.|
#'
#' Some function variants have different arguments:
#'
#' grid_W_ expects it's own `scheme_w` as depth levels are different between these and other files. it also expects it's
#' own `start_w` and `count_w`.
#'
#' @md
#' @param path the path to the NEMO-ERSEM model outputs.
#' @param file the name of a netcdf file containing the title variables.
#' @param scheme a summary scheme as expected by `array_w_mean()`.
#' @param start an optional vector of indices to start subsetting at. See ncdf4 documentation.
#' @param count an optional vector of steps to subset along. See ncdf4 documentation.
#' @param scheme_w as above but for grid_W files.
#' @param start_w as above but for grid_W files.
#' @param count_w as above but for grid_W files.
#' @param scheme2d a separate summary scheme for the 2D light variable.
#' @param ... soaks up unused function arguments passed by the wrapper functions handling file architecture.
#' @return A matrix with a column of group averages per title variable for a single day.
#' @family NEMO-ERSEM variable extractors
#' @name extractors_slabR
NULL

#' @rdname extractors_slabR
#' @export
get_grid_T_slabR <- function(path, file, scheme, start = c(1,1,1,1), count = c(-1,-1,-1,-1), ...) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_saline <- ncdf4::ncvar_get(nc_raw, "vosaline", start, count)      # Extract an array of salinities
  nc_temp <- ncdf4::ncvar_get(nc_raw, "votemper", start, count)        # Extract an array of temperatures
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                                # Bind as columns
    Salinity = array_w_mean(nc_saline, scheme),                                # Summarise salinity according to scheme
    Temperature = array_w_mean(nc_temp, scheme))                               # Summarise temperature according to scheme

    return(all)
}

#' @rdname extractors_slabR
#' @export
get_diad_slabR <- function(path, file, scheme, start = c(1,1,1,1), count = c(-1,-1,-1,-1), scheme2d, ...) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_light <- ncdf4::ncvar_get(nc_raw, "SurfaceE", start[-3], count[-3])       # Extract a matrix of surface irradiance
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                                # Bind as columns
    Light = c(nc_light[scheme2d],                                              # Subset light instead of average as it's 2D
              rep(NA, max(scheme$group)-length(scheme2d))))                    # Introduce NAs to fill deep layer summaries

  return(all)
}

#' @rdname extractors_slabR
#' @export
get_ptrc_T_slabR <- function(path, file, scheme, start = c(1,1,1,1), count = c(-1,-1,-1,-1), ...) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_NO3 <- ncdf4::ncvar_get(nc_raw, "N3n", start, count)                      # Extract an array for the variable
  nc_NH4 <- ncdf4::ncvar_get(nc_raw, "N4n", start, count)
#  nc_DET <- ncdf4::ncvar_get(nc_raw, "Q1n", start, count)
#  nc_PHD <- ncdf4::ncvar_get(nc_raw, "PHD", start, count)
#  nc_PHN <- ncdf4::ncvar_get(nc_raw, "PHN", start, count)
#  nc_Phyt <- nc_PHD + nc_PHN
  nc_POC1 <- ncdf4::ncvar_get(nc_raw, "R4c", start, count)
  nc_POC2 <- ncdf4::ncvar_get(nc_raw, "R8c", start, count)
  nc_Chl1 <- ncdf4::ncvar_get(nc_raw, "Chl1", start, count)
  nc_Chl2 <- ncdf4::ncvar_get(nc_raw, "Chl2", start, count)
  nc_Chl3 <- ncdf4::ncvar_get(nc_raw, "Chl3", start, count)
  nc_Chl4 <- ncdf4::ncvar_get(nc_raw, "Chl4", start, count)
  nc_Chl <- nc_Chl1 + nc_Chl2 + nc_Chl3 + nc_Chl4 ; rm(nc_Chl1, nc_Chl2, nc_Chl3, nc_Chl4)
  nc_POC <- nc_POC1 + nc_POC2 ; rm(nc_POC1, nc_POC2)
  ncdf4::nc_close(nc_raw)                                                          # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                                    # Bind as matrix
    NO3 = array_w_mean(nc_NO3, scheme),                                            # summarise dissolved nitrate according to scheme
    NH4 = array_w_mean(nc_NH4, scheme),                                            # summarise Dissolved ammonium according to scheme
    POC = array_w_mean(nc_POC, scheme),                                            # stand in for detritus
#    Detritus = array_w_mean(nc_DET, scheme),                                      # summarise Detritus according to scheme
    Chlorophyll = array_w_mean(nc_Chl, scheme))                                    # summarise Chlorophyll according to scheme
    return(all)
}

#' @rdname extractors_slabR
#' @export
get_grid_W_slabR   <- function(path, file, scheme_w, start_w = c(1,1,1,1), count_w = c(-1,-1,-1,-1), ...) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_vel <- ncdf4::ncvar_get(nc_raw, "vovecrtz", start_w, count_w) # Extract an array for the variable
  nc_dif <- ncdf4::ncvar_get(nc_raw, "votkeavt", start_w, count_w)
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                                # Bind as matrix
    Vertical_velocity = array_w_mean(nc_vel, scheme_w),                        # Summarise vertical velocity according to scheme
    Vertical_diffusivity = array_w_mean(nc_dif, scheme_w))                     # Summarise diffusivity according to scheme
  return(all)
}

#' @rdname extractors_slabR
#' @export
get_grid_V_slabR <- function(path, file, scheme, start = c(1,1,1,1), count = c(-1,-1,-1,-1), ...) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                               # Open up a netcdf file to see it's raw contents (var names)
  nc_merid <- ncdf4::ncvar_get(nc_raw, "vomecrty", start, count)         # Pull meridinal currents
  ncdf4::nc_close(nc_raw)                                                    # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(Meridional = array_w_mean(nc_merid, scheme))                  # summarise meridional currents according to scheme
  return(all)
}

#' @rdname extractors_slabR
#' @export
get_grid_U_slabR <- function(path, file, scheme, start = c(1,1,1,1), count = c(-1,-1,-1,-1), ...) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                               # Open up a netcdf file to see it's raw contents (var names)
  nc_zonal <- ncdf4::ncvar_get(nc_raw, "vozocrtx", start, count) # Pull zonal current
  ncdf4::nc_close(nc_raw)                                                    # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(Zonal = array_w_mean(nc_zonal, scheme))                       # summarise zonal currents according to scheme
  return(all)
}

