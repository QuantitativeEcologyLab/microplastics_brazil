# This script generates figure 1 in the main text that
# visualises differences in the study sites and general patterns
# the plastic concentrations
# Note: the indiviudal panels are generated in this script
# but need to be assembled outside of R

# Written by Michael Noonan


#Import the packages
library(ggplot2)
library(terra)
library(sf)
library(tidyterra)
library(gridExtra)
library(rphylopic)


#---------------------------------------------------------------------
# Data import and processing
#---------------------------------------------------------------------

#Import the MP datasets
source("scripts/data_import.R")


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

#Import the three biome boundaries obtained from: https://terrabrasilis.dpi.inpe.br/en/download-files/
biomes <- st_read("~/Dropbox/UBC/Projects/microplastics_brazil/data/environmental_data/biomes",
                  crs = 4326,
                  quiet = TRUE)
biomes <- st_transform(biomes, crs = crs(HFI))
biomes <- subset(biomes, NOM_BIOMA %in% c("Amaz\xf4nia", "Cerrado","Pantanal"))
amazon <- subset(biomes, NOM_BIOMA == "Amaz\xf4nia")
cerrado <- subset(biomes, NOM_BIOMA == "Cerrado")
pantanal <- subset(biomes, NOM_BIOMA == "Pantanal")


#Drop any individuals without movement data
mp_data_sub <- mp_data[!is.na(mp_data$hr),]

#Get each individuals' mean location
anteater <- mp_data_sub[mp_data_sub$species == "Myrmecophaga_tridactyla",c("name","biome", "Long", "Lat")]
anteater <- st_as_sf(anteater,
                     coords = c("Long",
                                "Lat"),
                     crs = st_crs(4326))
anteater <- st_transform(anteater, crs = crs(HFI))


tapir <- mp_data_sub[mp_data_sub$species == "Tapirus_terrestris",c("name","biome", "Long", "Lat")]
tapir <- st_as_sf(tapir,
                  coords = c("Long",
                             "Lat"),
                  crs = st_crs(4326))
tapir <- st_transform(tapir, crs = crs(HFI))

armadillo <- mp_data_sub[mp_data_sub$species == "Priodontes_maximus",c("name","biome", "Long", "Lat")]
armadillo <- st_as_sf(armadillo,
                      coords = c("Long",
                                 "Lat"),
                      crs = st_crs(4326))
armadillo <- st_transform(armadillo, crs = crs(HFI))


#Get the animal silhouettes 
tapir_pic <- get_phylopic("7950e979-6738-45b3-a7c6-c573ef5559d1")
anteater_pic <- get_phylopic("d52b48fc-be52-46a1-94b7-ac7790b4730c")
armadillo_pic <- get_phylopic("5d59b5ce-c1dd-40f6-b295-8d2629b9775e")



#---------------------------------------------------------------------
# Figure 1A - map of tracking data with underlying human footprint
#---------------------------------------------------------------------


#Generate the figure
a <-
  ggplot() +
  ggtitle("A") +
  geom_spatraster(data = HFI, maxcell = 5e+07,
                  alpha = 1) +
  scale_fill_gradient(name = "Human Footprint Index",
                      low = "#eeeeee",
                      high = "black",
                      na.value = NA) +
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA, barwidth = 4,
                               barheight = 0.2, direction = "horizontal")) +
  geom_sf(data = brasil, size = 0.2, fill = "transparent", linewidth = 0.1, col = "black") +
  #geom_sf(data = biomes, aes(fill = NOM_BIOMA), size = 0.2, linewidth = 0.1, col = "transparent", alpha = 0.3) +
  geom_sf(data = amazon, size = 0.2, fill = "#007200", linewidth = 0.1, col = "transparent", alpha = 0.3) +
  geom_sf(data = cerrado, size = 0.2, fill = "#e9c46a", linewidth = 0.1, col = "transparent", alpha = 0.3) +
  geom_sf(data = pantanal, size = 0.2, fill = "#0a9396", linewidth = 0.1, col = "transparent", alpha = 0.3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.position = "inside",
        legend.position.inside = c(0.28,0.12),
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
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  
  #Add locations of the animals
  add_phylopic(armadillo_pic,
               x = mean(st_coordinates(armadillo)[,1]),
               y = mean(st_coordinates(armadillo)[,2])+100000,
               ysize = 150000, alpha = 1, fill = "#bb3e03") +
  add_phylopic(tapir_pic,
               x = mean(st_coordinates(tapir[tapir$biome == "Pantanal",])[,1]),
               y = mean(st_coordinates(tapir[tapir$biome == "Pantanal",])[,2]),
               ysize = 150000, alpha = 1, fill = "#005f73",
               horizontal = TRUE) +
  add_phylopic(tapir_pic,
               x = mean(st_coordinates(tapir[tapir$biome == "Cerrado",])[,1]),
               y = mean(st_coordinates(tapir[tapir$biome == "Cerrado",])[,2]),
               ysize = 150000, alpha = 1, fill = "#005f73",
               horizontal = TRUE) +
  add_phylopic(tapir_pic,
               x = median(st_coordinates(tapir[tapir$biome == "Amazon",])[,1]),
               y = median(st_coordinates(tapir[tapir$biome == "Amazon",])[,2]),
               ysize = 150000, alpha = 1, fill = "#005f73",
               horizontal = TRUE) +
  add_phylopic(tapir_pic,
               x = median(st_coordinates(tapir[tapir$name == "Krishna",])[,1]),
               y = median(st_coordinates(tapir[tapir$name == "Krishna",])[,2]),
               ysize = 150000, alpha = 1, fill = "#005f73",
               horizontal = TRUE) +
  add_phylopic(anteater_pic,
               x = mean(st_coordinates(anteater)[,1]),
               y = mean(st_coordinates(anteater)[,2]),
               ysize = 150000, alpha = 1, fill = "#619b8a")


#Save the figure
ggsave(a,
       width = 2.375, height = 5, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/figure_1a.png")



#---------------------------------------------------------------------
# Figure 1B - Heatmap of polymer abundances



# First, reorder the 'Sample' column by its numeric component
ORDER <- metadata[order(metadata$species), "sample"]
data_long$sample <- factor(data_long$sample, levels = ORDER, ordered = TRUE)


#Adjust the polymer names
data_long$Microplastic <- stringr::str_to_title(data_long$Microplastic)
data_long$Microplastic[data_long$Microplastic == "Abs"] <- "ABS"
data_long$Microplastic[data_long$Microplastic == "Pet"] <- "PET"
data_long$Microplastic[data_long$Microplastic == "Polyvinyl Chloride"] <- "PVC"

# Heatmap
b <- 
  ggplot(data_long, aes(Microplastic, sample, fill= log(Adjusted+1))) + 
  ggtitle("B") +
  geom_tile(alpha = 0.95) +
  scico::scale_fill_scico(palette = "lipari",
                          name = "Particles/mL",
                          breaks = c(0,log(10),log(50),log(400)),
                          labels = c(0,10,50,400)) +
  geom_hline(yintercept = 21.5, linetype = "dashed", col = "grey70", linewidth = 0.3) +
  geom_hline(yintercept = 26.5, linetype = "dashed", col = "grey70", linewidth = 0.3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x  = element_text(size=4,
                                    family = "sans",
                                    angle = 90,
                                    face = "bold",
                                    color = "black",
                                    hjust = 1, vjust = 0.5),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold"),
        legend.position = "top",
        legend.title = element_text(size=5, family = "sans", face = "bold", vjust = -2, hjust = 0.5),
        legend.text = element_text(size=4, family = "sans", face = "bold", vjust = 4),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA, barwidth = 10,
                               barheight = 0.3, direction = "horizontal"))


#Save the figure
ggsave(b,
       width = 2.375, height = 5, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/figure_1b.png")


#---------------------------------------------------------------------
# Figure 1C - E - Radar plots of polymer abundances
#---------------------------------------------------------------------


mean_polymers <- aggregate(Adjusted ~ species + Microplastic, FUN = "mean", data = data_long)
mean_polymers <- reshape(data = mean_polymers,
                         idvar = "species",
                         timevar = "Microplastic",
                         v.names = "Adjusted",
                         direction = "wide")
#Adjust the polymer names
names(mean_polymers)[-1] <- sub("Adjusted.","",names(mean_polymers)[-1])
names(mean_polymers)[11] <- "PVC"


#Define the min and max (specific requirement of the fmsb package)
mean_polymers[4,] <- 40 #Max
mean_polymers[5,] <- 0 #Min

row.names(mean_polymers) <- c("Myrmecophaga_tridactyla",
                              "Priodontes_maximus",
                              "Tapirus_terrestris",
                              "Max",
                              "Min")

mean_polymers <- mean_polymers[c(4:5,1:3),]

# Generate and save the figures
png(filename = "figures/figure_1cde.png",
    width = 6, height = 2, units = "in",
    bg = "transparent",
    res = 600)

par(mfrow = c(1,3),
    mar = c(0.2,0.1,0.2,0.1),
    font.axis = 2,
    font = 2)
# Create the radar charts for the anteaters
radarchart(mean_polymers[c(1:2,3),-1],
           axistype = 1,
           pty = NA,
           pcol = "#619b8a",
           pfcol = adjustcolor("#619b8a", alpha.f = 0.3),
           plwd = 1,
           cglcol = "grey",
           cglty = 1,
           axislabcol = "grey",
           caxislabels = c(0, 10, 20, 30, 40),
           calcex = 0.5,
           cglwd = 0.4,
           vlcex = 0.48)
title(main = "C",
      cex.main = 0.8, font.main= 2, col.main= "black", adj = 0, line = -1)

# Create the radar chart for the armadillos
radarchart(mean_polymers[c(1:2,4),-1],
           axistype = 1,
           pty = NA,
           pcol = "#bb3e03",
           pfcol = adjustcolor("#bb3e03", alpha.f = 0.3),
           plwd = 1,
           cglcol = "grey",
           cglty = 1,
           axislabcol = "grey",
           caxislabels = c(0, 10, 20, 30, 40),
           calcex = 0.5,
           cglwd = 0.4,
           vlcex = 0.48)
title(main = "D", cex.main = 0.8, font.main= 2, col.main= "black", adj = 0, line = -1)

# Create the radar chart for the tapirs
radarchart(mean_polymers[c(1:2,5),-1],
           axistype = 1,
           pty = NA,
           pcol = "#005f73",
           pfcol = adjustcolor("#005f73", alpha.f = 0.3),
           plwd = 1,
           cglcol = "grey",
           cglty = 1,
           axislabcol = "grey",
           caxislabels = c(0, 10, 20, 30, 40),
           calcex = 0.5,
           cglwd = 0.4,
           vlcex = 0.48)
title(main = "E", cex.main = 0.8, font.main= 2, col.main= "black", adj = 0, line = -1)

dev.off()

