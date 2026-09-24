library(terra)
library(sf)
library(dplyr)

extracted_pts <- data.frame(
  Plot_ID   = sprintf("P%02d", 1:30),
  Land_Type = c(rep("Forest", 10), rep("Peri-Urban", 10), rep("Grassland", 10)),
  Longitude = c(85.35727, 85.35819, 85.35912, 85.36004, 85.36097, 85.35727, 85.36004, 85.36097, 85.36375, 85.36467,
                85.36004, 85.36097, 85.36375, 85.36364, 85.35819, 85.35912, 85.36004, 85.36097, 85.36282, 85.36375,
                85.35541, 85.35727, 85.36282, 85.36467, 85.36560, 85.35541, 85.36467, 85.36560, 85.36375, 85.36467),
  Latitude  = c(27.75346, 27.75346, 27.75346, 27.75346, 27.75346, 27.75256, 27.75256, 27.75256, 27.75166, 27.75166,
                27.75077, 27.75077, 27.75077, 27.74987, 27.74987, 27.74987, 27.74987, 27.74987, 27.74987, 27.74987,
                27.75436, 27.75436, 27.75436, 27.75436, 27.75436, 27.75346, 27.75346, 27.75346, 27.75256, 27.75256),
  pH_table        = c(5.59, 5.58, 5.65, 5.71, 5.66, 5.58, 5.68, 5.68, 5.64, 5.63,
                      5.65, 5.62, 5.65, 5.62, 5.63, 5.67, 5.60, 5.56, 5.61, 5.65,
                      5.56, 5.58, 5.62, 5.56, 5.59, 5.66, 5.46, 5.61, 5.61, 5.46),
  Clay_table      = c(18.55, 19.53, 18.11, 18.68, 18.12, 19.05, 18.13, 18.47, 17.28, 17.34,
                      17.44, 17.87, 17.72, 18.33, 18.24, 17.63, 18.75, 18.03, 17.14, 19.05,
                      19.54, 19.11, 17.94, 18.71, 18.14, 18.22, 17.20, 17.72, 17.83, 17.79),
  OM_table        = c(3.15, 3.11, 3.29, 3.37, 3.17, 3.17, 3.47, 3.60, 3.02, 3.17,
                      3.31, 3.42, 3.04, 3.59, 3.53, 3.40, 2.97, 3.15, 3.47, 3.13,
                      2.78, 2.79, 2.74, 2.91, 2.86, 2.84, 2.80, 2.68, 2.81, 2.99),
  N_table         = c(0.154, 0.149, 0.162, 0.173, 0.162, 0.158, 0.181, 0.175, 0.146, 0.151,
                      0.171, 0.174, 0.157, 0.183, 0.172, 0.181, 0.157, 0.150, 0.178, 0.151,
                      0.150, 0.148, 0.131, 0.135, 0.140, 0.148, 0.136, 0.139, 0.134, 0.145),
  P2O5_table      = c(127.89, 132.05, 144.06, 146.80, 147.97, 123.28, 162.25, 157.01, 140.16, 140.37,
                      153.38, 152.18, 136.28, 161.74, 170.62, 165.85, 146.02, 160.70, 174.55, 122.66,
                      131.04, 143.55, 157.32, 147.94, 142.24, 135.58, 130.04, 152.22, 130.56, 137.06),
  K_table         = c(291.93, 281.04, 279.85, 322.50, 301.26, 270.51, 307.25, 317.40, 261.33, 281.90,
                      293.30, 289.26, 284.56, 322.43, 325.50, 311.35, 287.45, 269.00, 301.61, 253.59,
                      272.24, 303.49, 247.44, 245.89, 261.83, 258.53, 250.54, 266.63, 257.71, 254.07)
)

# Convert coordinates to spatial points vector (EPSG:4326)
pts_vect <- vect(extracted_pts, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")

# Path setup for rasters
desktop_dir <- "C:/Users/user/OneDrive/Desktop"

raster_files <- c(
  pH   = file.path(desktop_dir, "ph.tif"),
  Clay = file.path(desktop_dir, "clay.tif"),
  K    = file.path(desktop_dir, "k.tif"),
  N    = file.path(desktop_dir, "nitrogen.tif"),
  OM   = file.path(desktop_dir, "organic.tif"),
  P2O5 = file.path(desktop_dir, "p2o5.tif")
)

# Extract raster values at plot locations and compare
results <- extracted_pts %>% select(Plot_ID, Land_Type)

for (var_name in names(raster_files)) {
  tif_path <- raster_files[[var_name]]
  
  if (file.exists(tif_path)) {
    r <- rast(tif_path)
    
    # Project spatial points if raster uses a different CRS (e.g., UTM Zone 45N)
    if (crs(r) != crs(pts_vect)) {
      pts_proj <- project(pts_vect, crs(r))
    } else {
      pts_proj <- pts_vect
    }
    
    extracted_val <- terra::extract(r, pts_proj)[[2]]
    
    tbl_col  <- paste0(var_name, "_table")
    tbl_vals <- extracted_pts[[tbl_col]]
    
    diff_val <- abs(extracted_val - tbl_vals)
    
    results[[paste0(var_name, "_table")]]  <- tbl_vals
    results[[paste0(var_name, "_raster")]] <- round(extracted_val, 2)
    results[[paste0(var_name, "_diff")]]   <- round(diff_val, 4)
  } else {
    warning(paste("File not found:", tif_path))
  }
}

# Check for any non-zero differences
cat("\n=== RASTER VERIFICATION SUMMARY ===\n")
diff_cols <- grep("_diff$", names(results), value = TRUE)

mismatches <- results %>%
  filter(if_any(all_of(diff_cols), ~ . > 0.01))

if (nrow(mismatches) == 0) {
  cat("SUCCESS: All extracted values match the raster pixel values at the given coordinates!\n")
} else {
  cat("DISCREPANCIES FOUND:\n")
  print(mismatches)
}
