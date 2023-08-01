library(rasterSp)
library(rgdal)
library(tidyverse)
library(magrittr)

# # Read and merge temporary files converted to shapefile format by QGIS (selects only presence = 1)
# botw_1 <- sf::st_read(dsn="C:\\bird\\spp_shp\\all_spp_pt1.shp", layer="all_spp_pt1")
# botw_2 <- sf::st_read(dsn="C:\\bird\\spp_shp\\all_spp_pt2.shp", layer="all_spp_pt2")
# botw_3 <- sf::st_read(dsn="C:\\bird\\spp_shp\\all_spp_pt3.shp", layer="all_spp_pt3")
# 
# sf::st_write(bind_rows(botw_1, botw_2, botw_3), dsn="C:\\bird\\spp_shp\\all_spp.shp")

rasterizeRange(dsn="C:\\bird\\spp_shp\\all_spp_pt1.shp",
                       id="sci_name", 
                       resolution=0.5,
                       save=TRUE,
                       path="C:\\bird\\spp_shp\\raster\\")

rasterizeRange(dsn="C:\\bird\\spp_shp\\all_spp_pt2.shp",
                       id="sci_name", 
                       resolution=0.5,
                       save=TRUE,
                       path="C:\\bird\\spp_shp\\raster\\")

rasterizeRange(dsn="C:\\bird\\spp_shp\\all_spp_pt3.shp",
                       id="sci_name", 
                       resolution=0.5,
                       save=TRUE,
                       path="C:\\bird\\spp_shp\\raster\\")

# This may also exist within rasterSp::ter_birds_dist, takes at least 3h
start_time <- Sys.time()
all_rich <- calcSR(species_names=c(list.files("C:/bird/spp_shp/raster") %>% gsub("_0.5.tif", "", .) %>% gsub("_", " ", .)), 
                                    path="C:/bird/spp_shp/raster/")
end_time <- Sys.time()


# rast <- rasterizeRange(dsn=list.files("C:\\bird\\spp_shp\\", pattern=".shp", full.names=TRUE), 
#                        id="sci_name", 
#                        resolution=0.5,
#                        save=TRUE,
#                        path="C:\\bird\\spp_shp\\raster\\")

