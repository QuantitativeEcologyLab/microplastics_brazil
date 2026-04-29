# This script calculates the change in HFI between 2001 and 2020
# Generates some summary statistics

# Written by Michael Noonan


#Load in any necessary packages
library(terra)
library(sf)
library(fields)


#---------------------------------------------------------------------
# Data import and processing
#---------------------------------------------------------------------

#Import the 2001 and 2020 HII rasters
#Note: Files are too large to store on github, but are available here: https://wcshumanfootprint.org/data-access
hii_2001 <- rast("~/Dropbox/UBC/Projects/microplastics_brazil/data/environmental_data/hii_2001-01-01.tif")/100
hii_2020 <- rast("~/Dropbox/UBC/Projects/microplastics_brazil/data/environmental_data/hii_2020-01-01.tif")/100
#Create a difference raster (2020 minus 2001)
hii_diff <- hii_2020 - hii_2001

#Get the boundaries of Brasil
world_sf <- st_as_sf(rworldmap::getMap(resolution = "low"))
brasil_sf <- subset(world_sf, SOVEREIGNT == "Brazil")
brasil <- vect(lwgeom::st_transform_proj(brasil_sf, crs = crs(hii_2001)))


#Import the three biome boundaries obtained from: https://terrabrasilis.dpi.inpe.br/en/download-files/
biomes <- st_read("~/Dropbox/UBC/Projects/microplastics_brazil/data/environmental_data/biomes",
                  crs = 4326,
                  quiet = TRUE)
biomes <- st_transform(biomes, crs = crs(hii_diff))
biomes <- subset(biomes, NOM_BIOMA %in% c("Amaz\xf4nia", "Cerrado","Pantanal"))
amazon <- subset(biomes, NOM_BIOMA == "Amaz\xf4nia")
cerrado <- subset(biomes, NOM_BIOMA == "Cerrado")
pantanal <- subset(biomes, NOM_BIOMA == "Pantanal")

#---------------------------------------------------------------------
# Global means and change
#---------------------------------------------------------------------

#Calculate mean HFI in 2001
global_mean_2001 <- global(hii_2001, fun = "mean", na.rm = TRUE)[[1]]
#Calculate mean HFI in 2020
global_mean_2020 <- global(hii_2020, fun =  "mean", na.rm = TRUE)[[1]]
#Calculate difference in mean HFI between 2001 and 2020
global_difference <- global(hii_2020, fun =  "mean", na.rm = TRUE)[[1]]

#Percent change
(global_mean_2020 - global_mean_2001)/global_mean_2001

#---------------------------------------------------------------------
# Brasil HFI change
#---------------------------------------------------------------------

# Crop and mask the rasters to Brasil
hii_2001_br <- crop(hii_2001, brasil)
hii_2001_br <- mask(hii_2001_br, brasil)
hii_2020_br <- crop(hii_2020, brasil)
hii_2020_br <- mask(hii_2020_br, brasil)
diff_br <- crop(hii_diff, brasil)
diff_br <- mask(diff_br, brasil)


#Calculate Brasil-wide means
brasil_mean_2001 <- global(hii_2001_br, fun = "mean", na.rm = TRUE)[[1]]
brasil_mean_2020 <- global(hii_2020_br, fun = "mean", na.rm = TRUE)[[1]]
brasil_difference <- global(diff_br, fun =  "mean", na.rm = TRUE)[[1]]

#Percent change
(brasil_mean_2020 - brasil_mean_2001)/brasil_mean_2001

#---------------------------------------------------------------------
# Amazon HFI change
#---------------------------------------------------------------------

# Crop and mask the rasters to the Amazon biome
hii_2001_amazon <- crop(hii_2001, amazon)
hii_2001_amazon <- mask(hii_2001_amazon, amazon)
hii_2020_amazon <- crop(hii_2020, amazon)
hii_2020_amazon <- mask(hii_2020_amazon, amazon)
diff_amazon <- crop(hii_diff, amazon)
diff_amazon <- mask(diff_amazon, amazon)


#Calculate Amazon-wide means
amazon_mean_2001 <- global(hii_2001_amazon, fun = "mean", na.rm = TRUE)[[1]]
amazon_mean_2020 <- global(hii_2020_amazon, fun = "mean", na.rm = TRUE)[[1]]
amazon_difference <- global(diff_amazon, fun =  "mean", na.rm = TRUE)[[1]]

#Percent change
(amazon_mean_2020 - amazon_mean_2001)/amazon_mean_2001

#---------------------------------------------------------------------
# Cerrado HFI change
#---------------------------------------------------------------------

# Crop and mask the rasters to the Cerrado biome
hii_2001_cerrado <- crop(hii_2001, cerrado)
hii_2001_cerrado <- mask(hii_2001_cerrado, cerrado)
hii_2020_cerrado <- crop(hii_2020, cerrado)
hii_2020_cerrado <- mask(hii_2020_cerrado, cerrado)
diff_cerrado <- crop(hii_diff, cerrado)
diff_cerrado <- mask(diff_cerrado, cerrado)


#Calculate Cerrado-wide means
cerrado_mean_2001 <- global(hii_2001_cerrado, fun = "mean", na.rm = TRUE)[[1]]
cerrado_mean_2020 <- global(hii_2020_cerrado, fun = "mean", na.rm = TRUE)[[1]]
cerrado_difference <- global(diff_cerrado, fun =  "mean", na.rm = TRUE)[[1]]

#Percent change
(cerrado_mean_2020 - cerrado_mean_2001)/cerrado_mean_2001

#---------------------------------------------------------------------
# Pantanal HFI change
#---------------------------------------------------------------------

# Crop and mask the rasters to the Pantanal biome
hii_2001_pantanal <- crop(hii_2001, pantanal)
hii_2001_pantanal <- mask(hii_2001_pantanal, pantanal)
hii_2020_pantanal <- crop(hii_2020, pantanal)
hii_2020_pantanal <- mask(hii_2020_pantanal, pantanal)
diff_pantanal <- crop(hii_diff, pantanal)
diff_pantanal <- mask(diff_pantanal, pantanal)


#Calculate Pantanal-wide means
pantanal_mean_2001 <- global(hii_2001_pantanal, fun = "mean", na.rm = TRUE)[[1]]
pantanal_mean_2020 <- global(hii_2020_pantanal, fun = "mean", na.rm = TRUE)[[1]]
pantanal_difference <- global(diff_pantanal, fun =  "mean", na.rm = TRUE)[[1]]

#Percent change
(pantanal_mean_2020 - pantanal_mean_2001)/pantanal_mean_2001