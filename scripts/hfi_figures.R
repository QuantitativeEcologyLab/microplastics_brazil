# This script generates figure 1 in the main text that
# visualises differences in the study sites and general patterns
# the plastic concentrations

# Written by Michael Noonan


#Import the packages
library(ggplot2)
library(ctmm)
library(terra)
library(sf)
library(tidyterra)
library(gridExtra)
library(rphylopic)
library(ggspatial)


#---------------------------------------------------------------------
# Raster import and processing
#---------------------------------------------------------------------


#Import the human footprint index raster from: https://www.frontiersin.org/articles/10.3389/frsen.2023.1130896/full
#Note: File is too large to store on github, but are available here: https://source.coop/repositories/vizzuality/hfp-100/description
HFI <- rast("~/Dropbox/UBC/Projects/microplastics_brazil/data/environmental_data/hfp_2021_100m_v1-2_cog.tif")

#Reproject and process the HFI data
world_sf <- st_as_sf(rworldmap::getMap(resolution = "low"))
brasil_sf <- subset(world_sf, SOVEREIGNT == "Brazil")
brasil <- vect(lwgeom::st_transform_proj(brasil_sf, crs = crs(HFI)))
HFI <- crop(HFI, brasil) #Note: can mask directly, but cropping first keeps it memory friendly
HFI <- mask(HFI, brasil) #Note: this step is slow
HFI <- HFI/1000 #rescale; also slow


#Import the land classification data for Brazil
#Data are from here: https://brasil.mapbiomas.org/en/
land_types <- rast("~/Dropbox/UBC/Projects/microplastics_brazil/data/environmental_data/brasil_coverage_2022.tif")


#---------------------------------------------------------------------
# Individual figure generation
#---------------------------------------------------------------------

#Identify the file paths to the gps data
gps_data_paths <- list.files("~/Dropbox/UBC/Projects/microplastics_brazil/data/gps_data",
                             full.names = TRUE)


for(i in 1:length(gps_data_paths)){
  
  DATA <- read.csv(gps_data_paths[i])
  DATA[which(DATA$ID == "Alvinho/Alvaro"),"ID"] <- "Alvinho_Alvaro"
  DATA <- as.telemetry(DATA, mark.rm = TRUE)
  
  #Identify the correct file path for importing the home range
  if(grepl("anteater", gps_data_paths[i])){HR_path <- paste("~/Dropbox/UBC/Projects/microplastics_brazil/results/anteater_movement/akdes", sep = "")}
  if(grepl("armadillo", gps_data_paths[i])){HR_path <- paste("~/Dropbox/UBC/Projects/microplastics_brazil/results/armadillo_movement/akdes", sep = "")}
  if(grepl("tapir", gps_data_paths[i])){HR_path <- paste("~/Dropbox/UBC/Projects/microplastics_brazil/results/tapir_movement/akdes", sep = "")}
  
  #Identify the current species for saving and figure colours
  if(grepl("anteater", gps_data_paths[i])){species <- "anteater"; COL <- "#619b8a"}
  if(grepl("armadillo", gps_data_paths[i])){species <- "armadillo"; COL <- "#bb3e03"}
  if(grepl("tapir", gps_data_paths[i])){species <- "tapir"; COL <- "#005f73"}
  
  for(j in 1:length(DATA)){
    
    #grab the jth individual
    cilla <- DATA[[j]]
    
    #Convert tracking data to sf format
    cilla_sf <- as.sf(cilla)
    cilla_sf <- st_transform(cilla_sf, crs = crs(HFI))
    
    #----------------------------------------
    #prep the HFI raster to match the extent of the data
    
    #Some calculations to buffer the data for nicer figures
    EXT <- ext(cilla_sf)
    
    #Calculate width and height
    width <- EXT[2] - EXT[1]  # xmax - xmin
    height <- EXT[4] - EXT[3] # ymax - ymin
    
    #Determine the larger dimension
    size <- max(width*1.2, height*1.2)
    
    #Center the extent to keep the same midpoint
    x_center <- (EXT[1] + EXT[2])/2
    y_center <- (EXT[3] + EXT[4])/2
    
    #Create a new square extent
    square_ext <- ext(x_center - size/2, x_center + size/2,
                      y_center - size/2, y_center + size/2)
    
    # cilla_sf_buff <- st_buffer(cilla_sf,
    #                            dist = mmax(width*.2, height*.2))
    
    #Crop HFI based on the buffered extent
    HFI_cropped <- crop(HFI,
                        square_ext,
                        snap = "out")
    
    # #Crop HFI based on the buffered extent
    # HFI_cropped <- crop(HFI,
    #                     cilla_sf_buff,
    #                     snap = "out")
    
    
    #----------------------------------------
    #prep the land class raster to match the extent of the data
    #Some calculations to buffer the data for nicer figures
    cilla_sf_2 <- st_transform(cilla_sf, crs = crs(land_types))
    #Some calculations to buffer the data for nicer figures
    EXT <- ext(cilla_sf_2)
    
    #Calculate width and height
    width <- EXT[2] - EXT[1]  # xmax - xmin
    height <- EXT[4] - EXT[3] # ymax - ymin
    
    #Determine the larger dimension
    size <- max(width*1.2, height*1.2)
    
    #Center the extent to keep the same midpoint
    x_center <- (EXT[1] + EXT[2])/2
    y_center <- (EXT[3] + EXT[4])/2
    
    #Create a new square extent
    square_ext <- ext(x_center - size/2, x_center + size/2,
                      y_center - size/2, y_center + size/2)
    
    # EXT <- ext(cilla_sf_2)
    # cilla_sf_buff <- st_buffer(cilla_sf_2,
    #                            dist = max(c((EXT[2] - EXT[1])*30000, (EXT[4] - EXT[3])*30000)))
    land_cropped <- crop(land_types,
                         square_ext,
                         snap = "out")
    
    #Process the land class raster for easy plotting
    land_cropped <- ifel(land_cropped %in% c(1,3,4,5,6,49,29), 99, land_cropped) #99 = Native_forest
    land_cropped <- ifel(land_cropped %in% c(12,15), 98, land_cropped) #98 = Pasture
    land_cropped <- ifel(land_cropped %in% c(18,19,20,21,39,40,41,62,36,46,47,35,48), 97, land_cropped) # 97 = Agriculture
    land_cropped <- ifel(land_cropped %in% c(24,25,50), 96, land_cropped) # 96 = Development
    land_cropped <- ifel(land_cropped %in% c(11,26,33), 95, land_cropped) # 95 = water
    land_cropped <- ifel(land_cropped %in% 0:94, 94, land_cropped) # 94 = other
    land_cropped <- as.factor(land_cropped)
    
    
    #----------------------------------------
    #Generate the HFI figure
    
    hfi_figure <- 
      ggplot() +
      geom_spatraster(data = HFI_cropped, maxcell = 5e+07,
                      alpha = 0.7) +
      geom_sf(data = cilla_sf, size = 0.2, alpha = 0.9, col = COL, shape = 16) +
      scale_fill_gradient(name = "Human Footprint Index",
                          low = "#eeeeee",
                          high = "black",
                          na.value = NA,
                          limits = c(0,50)) +
      guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA, barwidth = 4,
                                   barheight = 0.2, direction = "horizontal")) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_rect(fill = "transparent"),
            plot.background = element_rect(fill = "transparent", color = NA),
            legend.position = "top",
            legend.title = element_text(size=3, family = "sans", face = "bold", vjust = -4, hjust = 0.5),
            legend.text = element_text(size=2, family = "sans", face = "bold", vjust = 8),
            legend.background = element_rect(fill = "transparent"),
            legend.key.size = unit(0.2, 'cm'),
            legend.spacing.y = unit(0.1, 'cm'),
            plot.title = element_text(hjust = .01, vjust = -4, size = 6, family = "sans", face = "bold"),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x  = element_blank(),
            axis.ticks = element_blank(),
            strip.background=element_blank(),
            plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
      annotation_scale(height = unit(0.010, "npc"),
                       width_hint = 0.4,
                       line_width = 0.2,
                       pad_x = unit(0.07, "npc"),
                       pad_y = unit(0.07, "npc"),
                       text_pad = unit(0.01, "npc"),
                       text_cex = .2,
                       text_family = "sans",
                       text_face = "bold")
    
    
    
    #----------------------------------------
    #Generate the land class figure
    
    land_use_figure <- 
    ggplot() +
      geom_spatraster(data = land_cropped, maxcell = 5e+07,
                      alpha = 0.7, aes(fill = brasil_coverage_2022)) +
      geom_sf(data = cilla_sf_2, size = 0.2, alpha = 0.9, col = COL, shape = 16) +
      
      scale_fill_manual(breaks = c(94,95,96,97,98,99),
                          labels = c("Other", "Water", "Development", "Agriculture", "Pasture", "Native Forest"),
                          values = c("grey70", "#168aad", "#343a40", "#e9c46a","#99d98c","#31572c"),
                        name = "Land Class") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_rect(fill = "transparent"),
            plot.background = element_rect(fill = "transparent", color = NA),
            legend.position = "top",
            #legend.position.inside = c(0.28,0.12),
            legend.title = element_text(size=3, family = "sans", face = "bold", hjust = 0.5),
            legend.text = element_text(size=2, family = "sans", face = "bold"),
            legend.background = element_rect(fill = "transparent"),
            legend.key.size = unit(0.2, 'cm'),
            legend.spacing.y = unit(0.1, 'cm'),
            plot.title = element_text(hjust = .01, vjust = -4, size = 6, family = "sans", face = "bold"),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x  = element_blank(),
            axis.ticks = element_blank(),
            strip.background=element_blank(),
            plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
      annotation_scale(height = unit(0.010, "npc"),
                       width_hint = 0.4,
                       line_width = 0.2,
                       pad_x = unit(0.07, "npc"),
                       pad_y = unit(0.07, "npc"),
                       text_pad = unit(0.01, "npc"),
                       text_cex = .2,
                       text_family = "sans",
                       text_face = "bold")
    
    
    #Specify the file path for saving the hfi figure
    hfi_fig_path <- paste("~/Dropbox/UBC/Projects/microplastics_brazil/figures/tracking_data_hfi/",
                          species,
                          "/",
                          cilla@info$identity,
                          ".png",
                          sep = "")
    
    #Specify the file path for saving the land class figure
    land_use_fig_path <- paste("~/Dropbox/UBC/Projects/microplastics_brazil/figures/tracking_data_land_use/",
                               species,
                               "/",
                               cilla@info$identity,
                               ".png",
                               sep = "")
    
    #Save the figures
    ggsave(hfi_figure,
           width = 2.375, height = 2.5, units = "in",
           dpi = 600,
           bg = "transparent",
           file=hfi_fig_path)
    
    ggsave(land_use_figure,
           width = 2.375, height = 2.5, units = "in",
           dpi = 600,
           bg = "transparent",
           file=land_use_fig_path)
    
  } # closes the loop over the tracked individuals
} #closes the loop over the gps data file paths
