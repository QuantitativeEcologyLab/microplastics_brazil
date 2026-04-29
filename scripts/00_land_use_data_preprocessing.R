#This script converts the discrete land-use data on development available from MapBiomas
#To a continuous 'distance to development' rasters.
#It also projects the MapBiomas to EPSG:5641

#Note: This script only needs to be run once. The individual rasters are saved
#in the environmental_data subdirectory of the project.
#Due to file size limitations, these rasters were not uploaded to the github repository
#however, the base land classification data are openly from: https://brasil.mapbiomas.org/en/

# Written by Michael Noonan

#load packages
library(terra)
library(raster)
library(sf)



#-------------------------------------------------------------
# Import and process the land class information
#-------------------------------------------------------------

# Import the MapBiomas 2022 land classification raster for Brazil
# Land classification data are from here: https://brasil.mapbiomas.org/en/
land_types <- rast("~/Dropbox/UBC/Projects/microplastics_brazil/data/environmental_data/brasil_coverage_2022_collection10.tif")


# Reproject land use data from lat/long to SIRGAS 2000 / South America Albers Equal Area Conic (EPSG:5641)
# Note: This step is slow, but it speeds up subsequent calculations.
land_types <- project(land_types, "EPSG:5641", method = "near")
writeRaster(land_types,
            "~/Dropbox/UBC/Projects/microplastics_brazil/data/environmental_data/brasil_coverage_2022_projected.tif",
            #filetype = "COG",
            overwrite = TRUE)


#Create a raster for development
#Legend codes come from here: https://brasil.mapbiomas.org/wp-content/uploads/sites/4/2025/08/Legenda-Colecao-10-Legend-Code.pdf
development <- land_types
development <- ifel(development %in% c(24, 25, 30, 75), 1, NA)
writeRaster(development,
            "~/Dropbox/UBC/Projects/microplastics_brazil/data/environmental_data/development.tif",
            filetype = "COG",
            overwrite = TRUE)
