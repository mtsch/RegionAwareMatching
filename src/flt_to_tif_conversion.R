library(raster)
raster_folder = "/dhps_mapping_data/540_rasters" # where rasters are saved

for (year in seq(1990, 2020)){ # loop over years
  print(year)
  map <- raster(paste(getwd(),raster_folder,sprintf("/Resistance_median_%s.flt", year), sep = ""))
  map
  jpeg(file = paste(getwd(),raster_folder,sprintf("/Resistance_median_%s.jpeg", year), sep = ""))
  plot(map)
  dev.off()
  writeRaster(map, paste(getwd(),raster_folder,sprintf("/Resistance_median_%s.tif", year), sep = ""), overwrite = TRUE)
}
