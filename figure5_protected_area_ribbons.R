# Figure 5: Protected area coverage along ungulate movement trajectories
# GIUM manuscript -- illustrative analytical and visualization workflow.
# Movement and protected area datasets are required and are not included.
# Set data_path below to your local data directory, then adapt the relative
# filenames and folders in the input calls if your data are arranged differently.
# Expected structure under data_path:
#   protected_areas/GIUM_PAs_May2026v2.shp (and shapefile companion files)
#   movement_data/proj_Burchell_zebra_Botswana/sequenceShapefiles/
#   movement_data/proj_Eland_Serengeti/sequenceShapefiles/
#   movement_data/proj_mongoliangobi_Mongolia/pts_Khulan_SGobi_MNG.shp
#   plots/ (create this folder before enabling optional output commands)
# This script documents the supplied workflow; it is not a standalone analysis.
# Original calculations, colours, and plot layout are retained. The final layout
# retains indices for 34 migratory and 3 nomadic populations, while the example
# inputs below include only zebra, eland, and khulan. Adapt the population inputs
# and layout together when using this example.

######################################
## GIUM - Lines for ribbons         ##
## Nomads and migratory individuals ##
## Updated: 05 November 2026        ##
######################################

rm(list = ls())

# Set the path to your local data directory.
data_path <- "your/path/to/data"

library(sf)
library(tidyverse)
library(stringr)
library(mapview)
library(ggstream)
library(cowplot)
library(grid)
library(gridExtra)
library(ggspatial)
library(patchwork)
library(ggplotify)

# load funtions
#Helper to load and process shapefiles
read_points <- function(path, layers) {
  points <- lapply(layers, function(layer) {
    st_read(dsn = file.path(path, "sequenceShapefiles"), layer = layer, quiet = TRUE) %>%
      mutate(date = as.POSIXct(date, format = "%Y-%m-%d %H:%M:%S"))
  })
  do.call(rbind, points)
}

read_points2 <- function(path, layers) {
  points <- lapply(layers, function(layer) {
    readRDS(
      file.path(path, "sequences", layer, paste0(layer, ".rds"))
    )
  })
  do.call(rbind, points)
}

utm_crs_from_data <- function(sf_obj) {
  # Calculate centroid of the combined geometry
  centroid <- st_centroid(st_union(sf_obj))
  coords <- st_coordinates(centroid)
  lon <- coords[, 1]
  lat <- coords[, 2]
  
  # Compute UTM zone
  zone <- floor((lon + 180) / 6) + 1
  
  # EPSG: 32600 = Northern Hemisphere, 32700 = Southern Hemisphere
  if (lat >= 0) {
    epsg <- 32600 + zone
  } else {
    epsg <- 32700 + zone
  }
  
  return(epsg)
}


process_species <- function(proj_path, out_path, country, species, layers, population, movement_type) {
  message("Processing: ", species, " (", country, ")")
  # ---- Read & clean points ----
  dat <- read_points(proj_path, layers) %>%
    st_transform(4326) %>%  # ensure consistent CRS
    arrange(id, date) %>%
    st_make_valid() %>%
    group_by(id) %>%
    filter(n() > 1) %>%
    ungroup() %>%
    mutate(id_year = paste0(id, "_", year(date)),
           season = str_extract(id_year, "[^_]+$")) %>% 
    filter(is.na(date) == FALSE)
  
  # ---- Project to UTM for accurate distance ----
  utm_crs <- utm_crs_from_data(dat)
  
  dat_proj <- dat %>%
    st_transform(crs = utm_crs) %>%
    mutate(
      x = st_coordinates(.)[,1],
      y = st_coordinates(.)[,2]
    ) %>%
    arrange(id, date) %>%
    group_by(id_year) %>%
    mutate(
      step_dist = sqrt((x - lag(x))^2 + (y - lag(y))^2)
    ) %>%
    filter(is.na(step_dist) | step_dist > 200) %>%
    ungroup()
  
  # ---- Build movement lines ----
  dat_line <- dat_proj %>%
    group_by(id_year) %>%
    summarise(
      geometry = st_combine(geometry) %>% st_cast("LINESTRING"),
      firstdate = min(date, na.rm = TRUE),
      lastdate = max(date, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      id = str_extract(id_year, "^[^_]+"),
      year = str_extract(id_year, "_(\\d{2,4})") %>% str_remove("_"),
      season = str_extract(id_year, "[^_]+$"),
      species = species
    )
  
  # ---- Sample points every 1 km along filtered lines ----
  pts <- dat_line %>%
    st_transform(crs = utm_crs) %>%
    mutate(geometry = st_line_sample(., density = 1 / 1000)) %>%
    st_cast("POINT", warn = FALSE) %>%
    group_by(id_year) %>%
    mutate(position = row_number() - 1) %>%
    ungroup()
  
  pts <- pts %>%
    arrange(id_year, position) %>%
    st_transform(crs = 4326)
  
  # ---- PA overlay ----
  pts$PA_status <- ifelse(lengths(st_intersects(pts, PA)) > 0, "PA", "NoPA")
  
  # ---- Plot PA + trajectories ----
  xlimits <- st_bbox(pts)[c("xmin", "xmax")]
  ylimits <- st_bbox(pts)[c("ymin", "ymax")]
  
  ggplot() +
    geom_sf(data = PA, fill = "lightgreen", color = "darkgreen", alpha = 0.5) +
    geom_sf(data = pts, aes(color = PA_status), size = 1) +
    geom_path(
      data = as.data.frame(st_coordinates(pts)) %>%
        mutate(id_year = pts$id_year, PA_status = pts$PA_status),
      aes(x = X, y = Y, group = id_year, color = PA_status),
      linewidth = 0.5, alpha = 0.5
    ) +
    scale_color_manual(values = c("PA" = "#5ab4ac", "NoPA" = "#d8b365"),breaks = c("PA", "NoPA"),labels = c("PA", "NoPA")) +
    theme_minimal() +
    labs(
      title = paste(species, country, "Sampled Points Every 1 km with PA Status"),
      color = "PA Status"
    ) +
    theme(legend.position = "bottom") +
    annotation_scale(location = "bl", width_hint = 0.3) +
    coord_sf(xlim = xlimits, ylim = ylimits)
  
  # ggsave(d, filename = file.path(out_path, paste0("PA_map_", species, "_", country, "V2.png")),
  #        width = 8, height = 6, units = "in", dpi = 300)
  
  # ---- Add lat/long and factor levels ----
  pts <- pts %>%
    ungroup() %>%
    mutate(
      PA_status = factor(PA_status, levels = pa_levels),
      country = country,
      species = species,
      lat = st_coordinates(.)[,2],
      long = st_coordinates(.)[,1]
    ) %>%
    st_drop_geometry()
  
  # ---- Annual PA percentages per individual ----
  individual_PA_year <- pts %>%
    filter(!is.na(PA_status)) %>%
    group_by(id, id_year) %>%
    summarise(
      population = first(population),
      movement_type = first(movement_type),
      species = first(species),
      country = first(country),
      n_PA = sum(PA_status == "PA"),
      n_NoPA = sum(PA_status == "NoPA"),
      n_points = n(),
      percent_in_PA = 100 * n_PA / n_points,
      .groups = "drop"
    )
  
  # ---- Average annual PA percentages per individual ----
  individual_PA <- individual_PA_year %>%
    group_by(id) %>%
    summarise(
      population = first(population),
      movement_type = first(movement_type),
      species = first(species),
      country = first(country),
      n_years = n(),
      percent_in_PA = mean(percent_in_PA, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Save the individual-level results for this population
  # write_csv(
  #   individual_PA,
  #   file.path(
  #     out_path,
  #     paste0("individual_PA_", make.names(population), ".csv")
  #   )
  # )
  
  # ---- Proportional summaries ----
  PA_prop <- pts %>%
    summarise(.by = c(position, PA_status),
              count = n(),
              species = unique(species),
              country = unique(country))
  
  PA_summary <- pts %>%
    count(PA_status, name = "count") %>%
    mutate(
      percent = round(count / sum(count) * 100),
      percent = if_else(percent == 0, "<1", as.character(percent)),
      label = PA_status
    ) %>%
    select(-count)
  
  PA_result <- left_join(PA_prop, PA_summary, by = "PA_status")
  
  # --- Overall % in PA: mean of individual-level percentages ---
  overall_pa_summary <- individual_PA %>%
    summarise(
      overall_percent_in_PA = mean(percent_in_PA, na.rm = TRUE)
    ) %>%
    pull(overall_percent_in_PA)
  
  # Add it to IUCN_result as a new column
  PA_result <- PA_result %>%
    mutate(overall_percent_in_PA = overall_pa_summary)
  
  return(PA_result)  
}

process_nomad_species <- function(proj_path, out_path, country, species, population, movement_type) {
  message("Processing: ", species, " (", country, ")")
  # ---- Read & clean points ----
  dat <- st_read(proj_path) %>%
    mutate(
      timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S"),
      date = as.Date(timestamp),
      hour = hour(timestamp)
    ) %>%
    #filter(hour == 12) %>%
    group_by(tag_id) %>%
    arrange(date) %>%
    # get first date per tag_id
    mutate(start_date = first(date),
           end_date = start_date + days(364)) %>%
    # filter to keep only dates within that 365-day window
    filter(date >= start_date & date <= end_date) %>%
    ungroup() %>% 
    st_transform(crs = 4326)
  
  dat_line <- dat %>%
    group_by(tag_id) %>%
    summarise(
      geometry = st_combine(geometry) %>% st_cast("LINESTRING"),
      firstdate = min(date, na.rm = TRUE),
      lastdate = max(date, na.rm = TRUE)
    ) %>%
    mutate(
      id = tag_id,
      year = year(firstdate),
      species = species
    )
  
  # Determine UTM CRS for transformation
  utm_crs <- utm_crs_from_data(dat_line)
  
  pts <- dat_line %>%
    st_transform(crs = utm_crs) %>%
    mutate(geometry = st_line_sample(., density = 1 / 1000)) %>%
    st_cast("POINT", warn = FALSE) %>%
    group_by(id, year) %>%
    mutate(position = row_number()-1)
  
  pts <- pts %>% 
    st_transform(crs = 4326)
  
  # ---- Plot PA + trajectories ----
  xlimits <- st_bbox(pts)[c("xmin", "xmax")]
  ylimits <- st_bbox(pts)[c("ymin", "ymax")]
  
  ggplot() +
    geom_sf(data = PA, fill = "lightgreen", color = "darkgreen", alpha = 0.5) +
    geom_sf(data = pts, aes(color = PA_status), size = 1) +
    geom_path(
      data = as.data.frame(st_coordinates(pts)) %>%
        mutate(id_year = pts$id_year, PA_status = pts$PA_status),
      aes(x = X, y = Y, group = id_year, color = PA_status),
      linewidth = 0.5, alpha = 0.5
    ) +
    scale_color_manual(values = c("PA" = "#5ab4ac", "NoPA" = "#d8b365"),breaks = c("PA", "NoPA"),labels = c("PA", "NoPA")) +
    theme_minimal() +
    labs(
      title = paste(species, country, "Sampled Points Every 1 km with PA Status"),
      color = "PA Status"
    ) +
    theme(legend.position = "bottom") +
    annotation_scale(location = "bl", width_hint = 0.3) +
    coord_sf(xlim = xlimits, ylim = ylimits)
  
  # ggsave(d, filename = file.path(out_path, paste0("PA_map_", species, "_", country, "V2.png")),
  #        width = 8, height = 6, units = "in", dpi = 300)
  
  # Join with PA data
  pts <- pts %>%
    st_join(PA["NAME"], left = TRUE) %>%
    mutate(PA_status = if_else(is.na(NAME), "NoPA", "PA"),
           PA_status = factor(PA_status, levels = pa_levels),
           #IUCN = if_else(is.na(IUCN_CAT), "NoPA", IUCN_CAT),
           #IUCN = factor(IUCN, levels = iucn_levels),
           country = country,
           species = species,
           lat = st_coordinates(.)[, 2],
           long = st_coordinates(.)[, 1]    
    ) %>%
    st_drop_geometry()
  
  # ---- Individual-level PA summaries ----
  individual_PA <- pts %>%
    filter(!is.na(PA_status)) %>%
    group_by(id) %>%
    summarise(
      population = population,
      movement_type = movement_type,
      species = first(species),
      country = first(country),
      n_PA = sum(PA_status == "PA"),
      n_NoPA = sum(PA_status == "NoPA"),
      n_points = n(),
      percent_in_PA = 100 * n_PA / n_points,
      .groups = "drop"
    )
  
  # Save the individual-level results for this population
  # write_csv(
  #   individual_PA,
  #   file.path(
  #     out_path,
  #     paste0("individual_PA_", make.names(population), ".csv")
  #   )
  # )
  # 
  # Proportional calculation
  IUCN_prop <- pts %>%
    summarise(.by = c(position, PA_status), count = n(),
              #summarise(.by = c(day_of_year, PA_status), count = n(),
              species = unique(species), country = unique(country))
  
  IUCN_summary <- pts %>%
    #count(IUCN, name = "count") %>%
    count(PA_status, name = "count") %>%
    mutate(
      percent = round(count / sum(count) * 100),
      percent = if_else(percent == 0, "<1", as.character(percent)),
      #label = paste0(IUCN, " (", percent, "%)"),
      label = PA_status
    ) %>%
    select(-count)
  
  IUCN_result <- left_join(IUCN_prop, IUCN_summary, by = "PA_status")
  
  # --- Overall % in PA: mean of individual-level percentages ---
  overall_pa_summary <- individual_PA %>%
    summarise(
      overall_percent_in_PA = mean(percent_in_PA, na.rm = TRUE)
    ) %>%
    pull(overall_percent_in_PA)
  
  # Add it to IUCN_result as a new column
  IUCN_result <- IUCN_result %>%
    mutate(overall_percent_in_PA = overall_pa_summary)
  return(IUCN_result)  
}

make_stream_plot <- function(df, species_name, n_grid = 10000, show_legend = FALSE) {
  df <- df %>%
    mutate(
      position_scaled = 100 * (position - min(position, na.rm = TRUE)) /
        (max(position, na.rm = TRUE) - min(position, na.rm = TRUE))
    )
  ggplot(df) +
    geom_stream(aes(x = position_scaled, y = count, fill = PA_status),
                type = "proportional", n_grid = n_grid, sorting = 'none') +
    annotate("rect",
             xmin = 0, xmax = 100, ymin = 0, ymax = 1,
             color = "#b5b6b6", fill = NA, linetype = "solid") +
    scale_x_continuous(
      name = "Protection status along a selected seasonal migration direction",
      breaks = scales::pretty_breaks(n = 5),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_fill_manual("Protection status", values = pa_colors, labels = pa_labels) +
    theme_classic(base_size = 10) +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      plot.title = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )
}

combine_for_legend <- function() {
  # Create a dummy dataframe just to pull out the IUCN categories for the legend
  combined_data <- bind_rows(
    IUCN_zebra
  ) %>%
    select(PA_status, position,count) %>%
    #group_by(IUCN) %>% 
    distinct(PA_status,.keep_all=T) %>%
    #slice(n=1) %>% 
    mutate(species = "Legend")  # Add a dummy species label for legend only
  
  return(combined_data)
}

get_legend2 <- function(plot, legend = NULL) {
  if (is_ggplot(plot)) {
    gt <- ggplotGrob(plot)
  } else {
    if (is.grob(plot)) {
      gt <- plot
    } else {
      stop("Plot object is neither a ggplot nor a grob.")
    }
  }
  pattern <- "guide-box"
  if (!is.null(legend)) {
    pattern <- paste0(pattern, "-", legend)
  }
  indices <- grep(pattern, gt$layout$name)
  not_empty <- !vapply(
    gt$grobs[indices], 
    inherits, what = "zeroGrob", 
    FUN.VALUE = logical(1)
  )
  indices <- indices[not_empty]
  if (length(indices) > 0) {
    return(gt$grobs[[indices[1]]])
  }
  return(NULL)
}

#Load Protected Areas once
PA <- st_read(file.path(data_path, "protected_areas", "GIUM_PAs_May2026v2.shp")) %>%
  st_make_valid() %>%
  st_transform(crs = 4326) %>% 
  mutate(IUCN_CAT = if_else(is.na(IUCN_CAT), "NoPA", IUCN_CAT))

# PA status and color palette
pa_levels <- c("PA","NoPA")

pa_colors <- c(
  "PA" = "#5ab4ac", "NoPA" = "#d8b365")

pa_labels <- c(
  "PA" = "Within PA",
  "NoPA" = "Outside PA"
)

# Process each species
IUCN_zebra <- process_species(
  proj_path = file.path(data_path, "movement_data", "proj_Burchell_zebra_Botswana"),
  out_path = file.path(data_path, "plots"),
  country = "Botswana",
  species = "Zebra",
  population = "Zebra_Botswana",
  movement_type = "Migratory",
  layers = c("preshift1_points", "shift1_points")
)

make_stream_plot(IUCN_zebra, show_legend = FALSE)

IUCN_eland <- process_species(
  proj_path = file.path(data_path, "movement_data", "proj_Eland_Serengeti"),
  out_path = file.path(data_path, "plots"),
  country = "Tanzania",
  species = "Eland",
  layers = c("preshift1_points", "shift1_points"),
  population = "Eland_Serengeti",
  movement_type = "Migratory"
)
make_stream_plot(IUCN_eland, show_legend = FALSE)

# List of dataframes with names

IUCN_mig_list <- list(
  "Plains zebra, Makgadikgadi (BWA)" = IUCN_zebra,
  "Eland, Serengeti NP (TZA/KEN)" = IUCN_eland
)

species_by_PA_percent_mig <- IUCN_mig_list %>%
  imap_dfr(~ tibble(species = .y, max_pos = max(.x$overall_percent_in_PA))) %>%
  arrange(desc(max_pos)) %>%
  pull(species)

# Create plots without legend
dist_plots <- imap(IUCN_mig_list, ~ make_stream_plot(.x, .y, show_legend = FALSE))

# Reorder plots list
dist_plots <- dist_plots[species_by_PA_percent_mig]

# Create a plot for the legend
legend_data <- combine_for_legend()

legend_plot <- make_stream_plot(legend_data, "Legend", n_grid=1000,show_legend = TRUE)

shared_legend <- get_legend2(legend_plot)

# Nomads

IUCN_khulan <- process_nomad_species(
  proj_path = file.path(data_path, "movement_data", "proj_mongoliangobi_Mongolia", "pts_Khulan_SGobi_MNG.shp"),
  out_path = file.path(data_path, "plots"),
  country = "Mongolia",
  species = "Khulan",
  population = "SGobi",
  movement_type = "Nomad"
)

# List of dataframes with names
IUCN_nomad_list <- list(
  "Khulan, Mongolian Gobi (MNG)" = IUCN_khulan)

species_by_PA_percent_nomad <- IUCN_nomad_list %>%
  imap_dfr(~ tibble(species = .y, max_pos = max(.x$overall_percent_in_PA))) %>%
  arrange(desc(max_pos)) %>%
  pull(species)

# Create plots without legend
doy_plots <- imap(IUCN_nomad_list, ~ make_stream_plot(.x, .y, show_legend = FALSE))

# Reorder plots list
doy_plots <- doy_plots[species_by_PA_percent_nomad]

### SECTION HEADINGS
heading_grob <- function(text) {
  ggdraw(textGrob(text, gp = gpar(fontface = "bold", fontsize = 10), x = 0.02, just = "left"))
}

### LEFT LABELS
species_labels_dist <- names(dist_plots)

label_grobs_dist <- lapply(species_labels_dist, function(label) {
  ggdraw(textGrob(label, x = 0.05, just = "left", gp = gpar(fontface = "plain", fontsize = 10)))
})

species_labels_doy <- names(doy_plots)

label_grobs_doy <- lapply(species_labels_doy, function(label) {
  ggdraw(textGrob(label, x = 0.05, just = "left", gp = gpar(fontface = "plain", fontsize = 10)))
})


## RIGHT LABELS
IUCN_mig_list <- IUCN_mig_list[species_by_PA_percent_mig]
IUCN_nomad_list <- IUCN_nomad_list[species_by_PA_percent_nomad]
IUCN_list <- c(IUCN_mig_list,IUCN_nomad_list)

species_summary <- lapply(IUCN_list, function(df) {
  df %>%
    slice(1) %>%  # Get the first row only, since n_individuals and max_distance_km are the same in all rows
    transmute(
      `PA (%)` = overall_percent_in_PA,
      )
}) %>%
  bind_rows(.id = "species_label")  # Use list names as labels

species_summary_table <- lapply(IUCN_list, function(df) {
  df %>%
    slice(1) %>%  # Get the first row only, since n_individuals and max_distance_km are the same in all rows
    transmute(
      `PA (%)` = overall_percent_in_PA
    )
}) %>%
  bind_rows(.id = "species_label")  # Use list names as labels

right_labels <- species_summary %>%
  mutate(
    label = ifelse(
      `PA (%)` < 1,
      sprintf("%.1f%%", `PA (%)`),
      sprintf("%.0f%%", `PA (%)`)
    )
  ) %>%
  pull(label)

right_label_grobs <- lapply(right_labels, function(label) {
  ggdraw(textGrob(label, x = 0.95, just = "right", gp = gpar(fontface = "italic", fontsize = 10)))
})

# --- Create tableGrob for right side ---
right_table_grob <- tableGrob(
  species_summary,
  rows = NULL,                 # remove row numbers
  theme = ttheme_default(
    core = list(fg_params = list(cex = 0.8)),
    colhead = list(fg_params = list(fontface = "bold", cex = 0.85))
  )
)

# Convert tableGrob to ggplot object
right_table_gg <- as.ggplot(right_table_grob)

# Separate labels for DOY and DIST
left_labels_doy  <- plot_grid(plotlist = label_grobs_doy, ncol = 1,rel_heights = rep(1, 1))
left_labels_dist <- plot_grid(plotlist = label_grobs_dist, ncol = 1, rel_heights = rep(1, 2))  # adjust if needed

right_labels_doy  <- plot_grid(plotlist = right_label_grobs[35:37], ncol = 1,rel_heights = rep(1, 1))
right_labels_dist <- plot_grid(plotlist = right_label_grobs[1:34], ncol = 1,rel_heights = rep(1, 2))

### STACKED PLOT GRIDS
doy_stack  <- plot_grid(plotlist = doy_plots, ncol = 1, align = "v", rel_heights = rep(1, 1))
dist_stack <- plot_grid(plotlist = dist_plots, ncol = 1, align = "v", rel_heights = rep(1, 2))

# Align each block with left labels
doy_grid <- plot_grid(
  left_labels_doy,
  doy_stack,
  right_labels_doy,
  ncol = 3,
  rel_widths = c(0.47, 1, 0.06)
)

#print(doy_grid)

doy_grid_with_xlab <- plot_grid(
  doy_grid,
  ggdraw() + draw_label("Protection status during a full annual movement cycle", x = 0.6, fontface = "bold", size = 11),
  ncol = 1,
  rel_heights = c(1, 0.2)
)

#print(doy_grid_with_xlab)

dist_grid <- plot_grid(
  left_labels_dist,
  dist_stack,
  right_labels_dist,
  ncol = 3,
  rel_widths = c(0.47, 1, 0.06)
)

dist_grid_with_xlab <- plot_grid(
  dist_grid,
  ggdraw() + draw_label("Protection status along a selected seasonal migration direction", x = 0.6,fontface = "bold", size = 11),
  ncol = 1,
  rel_heights = c(1, 0.015)
)

### FINAL LAYOUT: FIX PLOT HEIGHTS
final_main <- plot_grid(
  heading_grob("A) Seasonal migrants"),
  dist_grid_with_xlab,
  heading_grob("B) Nomadic populations"),
  doy_grid_with_xlab,
  ncol = 1,
  rel_heights = c(1, 37, 1, 4.2)  # Adjust the bottom slightly for the added label
)

### COMBINE WITH LEGEND
final_plot <- plot_grid(
  final_main,
  shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.1)
)

print(final_plot)

# ggsave(
#   filename = "ribbon_plots.jpeg",
#   plot = final_plot,
#   width = 10,
#   height = 15,
#   dpi = 600
# )



