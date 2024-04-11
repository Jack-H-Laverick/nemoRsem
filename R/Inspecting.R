
#' Prepare for Averaging by Decade
#'
#' This function cleans the saved NEMO-ERSEM monthly summaries, for averaging into decades.
#'
#' @param saved A dataframe containing a summarised month from NEMO-ERSEM model outputs.
#' @return A dataframe containing a summarised month of NEMO-ERSEM output, gaining a decade column, and dropping columns
#' which aren't needed for spatial maps.
#' @family NEMO-ERSEM averages
#' @export
decadal <- function(saved) {

  import <- readRDS(file = saved) %>%                                   # Read in wide format data file
#    dplyr::select(-c(weights, Bathymetry)) %>%
    dplyr::select(-c(weights)) %>%
    dplyr::rename(Decade = Year)

  stringr::str_sub(import$Decade, -1, -1) <- "0"                        # Overwite the 4th digit with a 0 to get the decade

  return(import)
}

#' Plot Temporal Summaries
#'
#' This function builds a time series, split by model compartment.
#'
#' @param var A "quoted" column name denoting how to colour the map.
#' @return The function creates a time series by model compartment of a summarised variable.
#' The plot is saved into the NEMO-ERSEM folder.
#' @family NEMO-ERSEM plots
#' @export
ts_plot <- function(var) {

  ts <- ggplot2::ggplot(TS, aes(x=date, y= get(var), colour = Compartment)) +
    ggplot2::geom_line(size = 0.2) +
    ggplot2::geom_smooth(span = 0.08, size = 0.2, se = FALSE) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top") +
    ggplot2::labs(caption = paste("NE", var, "time series by compartment"), y = var) +
    NULL
  ggplot2::ggsave(paste0("./Figures/NEMO-ERSEM/TS_", var, ".png"), plot = ts, width = 16, height = 10,
                  units = "cm", dpi = 500, bg = "white")

}

#' Map Spatial Summaries
#'
#' This function builds a map of the passed variable, facetted by month.
#'
#' @param data A dataframe containing the columns "Zonal" and "Meridional" for currents.
#' @param var A "quoted" column name denoting how to colour the map.
#' @return The function creates a facetted map of a summarised variable. Each facet shows a month.
#' The plot is saved into the NEMO-ERSEM/grids folder.
#' @family NEMO-ERSEM plots
#' @export
point_plot <- function(data, var) {

  decade <- data$Decade[1]; depth <- data$slab_layer[1]                             # Find out what the data is

  if(depth == "D" & var %in% c("Ice", "Ice_conc", "Ice_Thickness", "Snow_Thickness")) {
    print("Skipped deep ice plot") } else {

      print(paste("plotting", var, "for", decade, depth))                          # Show things are working

      map <- ggplot2::ggplot() +                                                            # Create the base
        ggplot2::theme_minimal() +
        ggplot2::labs(title = paste("Decade:", decade),
                      subtitle = paste("Water layer:", depth), x = NULL, y = NULL) +
        ggplot2::geom_raster(data = data, aes(x=x, y=y, fill = get(var))) +
        viridis::scale_fill_viridis(option = "viridis", name = var, na.value = "red") +
        ggplot2::facet_wrap(vars(Month)) +
        NULL
      ggplot2::ggsave(paste0("./Figures/NEMO-ERSEM/grids/map ", var, " ", depth, " ", decade, ".png"),
                      plot = map, scale = 1, width = 32, height = 20, units = "cm", dpi = 500, bg = "white")
    }
}

