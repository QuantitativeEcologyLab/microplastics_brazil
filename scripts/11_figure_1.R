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
source("scripts/04_microplastics_data_import.R")


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
                      na.value = NA,
                      limits = c(0,50)) +
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA, barwidth = 5,
                               barheight = 0.3, direction = "horizontal")) +
  geom_sf(data = brasil, size = 0.2, fill = "transparent", linewidth = 0.1, col = "black") +
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
        legend.title = element_text(size=5, family = "sans", face = "bold", vjust = -2, hjust = 0.5),
        legend.text = element_text(size=4, family = "sans", face = "bold", vjust = 6),
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
       file="figures/main_text_panels/figure_1a.png")



#---------------------------------------------------------------------
# Figure 1B - Heatmap of polymer abundances



# First, reorder the 'Sample' column by its numeric component
ORDER <- mp_data[order(mp_data$species), "sample"]
mp_data_long$sample <- factor(mp_data_long$sample, levels = ORDER, ordered = TRUE)


# Heatmap
b <- 
  ggplot(mp_data_long, aes(polymer, sample, fill= log(concentration+1))) + 
  ggtitle("B") +
  geom_tile(alpha = 0.95) +
  scico::scale_fill_scico(palette = "lipari",
                          name = "Particles/mL",
                          breaks = c(0,log(10),log(50),log(400)),
                          labels = c(0,10,50,400)) +
  geom_hline(yintercept = 21.5, linetype = "dashed", col = "grey70", linewidth = 0.3) +
  geom_hline(yintercept = 24.5, linetype = "dashed", col = "grey70", linewidth = 0.3) +
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
       file="figures/main_text_panels/figure_1b.png")


#---------------------------------------------------------------------
# Figure 1C - E - Radar barplots of polymer abundances across species
#---------------------------------------------------------------------


mean_polymers <- aggregate(concentration ~ species + polymer,
                           FUN = "mean", data = mp_data_long)


C <- 
  ggplot(mean_polymers[mean_polymers$species == "Tapirus_terrestris",],
         aes(x = polymer, y = concentration)) +
  ggtitle("C") +
    geom_segment(data = data.frame(polymer = unique(mean_polymers$polymer)),
                 aes(x = polymer, xend = polymer, y = -8, yend = 40),
                 color = "grey85", linewidth = 0.2, inherit.aes = FALSE) +
  geom_col(width = 0.85, show.legend = FALSE, fill = "#005f73") +
  scale_y_continuous(limits = c(-8, 43),
                     breaks = c(0, 10, 20, 30, 40),
                     expand = c(0, 0)) +
  coord_polar(theta = "x", start = pi / 1) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = c("grey85", "grey85", "grey85","grey85","grey85", NA), linewidth = 0.2),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_text(size=5, family = "sans", face = "bold"),
        axis.title.x = element_text(size=3, family = "sans", face = "bold"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=3, family = "sans", face = "bold", color = "black"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold", vjust = -6),
        axis.ticks = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0,0.1,0,0.1), "cm")) +
  ggplot2::annotate("text",
                    x = pi*2.2,
                    y = c(10, 20, 30, 40),
                    label = c("10", "20", "30", "40"),
                    color = "black",
                    family = "sans",
                    size = 1)


D <- 
  ggplot(mean_polymers[mean_polymers$species == "Priodontes_maximus",],
         aes(x = polymer, y = concentration)) +
  ggtitle("D") +
  geom_segment(data = data.frame(polymer = unique(mean_polymers$polymer)),
               aes(x = polymer, xend = polymer, y = -8, yend = 40),
               color = "grey85", linewidth = 0.2, inherit.aes = FALSE) +
  geom_col(width = 0.85, show.legend = FALSE, fill = "#bb3e03") +
  scale_y_continuous(limits = c(-8, 43),
                     breaks = c(0, 10, 20, 30, 40),
                     expand = c(0, 0)) +
  coord_polar(theta = "x", start = pi / 1) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = c("grey85", "grey85", "grey85","grey85","grey85", NA), linewidth = 0.2),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_text(size=5, family = "sans", face = "bold"),
        axis.title.x = element_text(size=3, family = "sans", face = "bold"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=3, family = "sans", face = "bold", color = "black"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold", vjust = -6),
        axis.ticks = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0,0.1,0,0.1), "cm")) +
  ggplot2::annotate("text",
                    x = pi*2.2,
                    y = c(10, 20, 30, 40),
                    label = c("10", "20", "30", "40"),
                    color = "black",
                    family = "sans",
                    size = 1)

E <- 
  ggplot(mean_polymers[mean_polymers$species == "Myrmecophaga_tridactyla",],
         aes(x = polymer, y = concentration)) +
  ggtitle("E") +
  geom_segment(data = data.frame(polymer = unique(mean_polymers$polymer)),
               aes(x = polymer, xend = polymer, y = -8, yend = 40),
               color = "grey85", linewidth = 0.2, inherit.aes = FALSE) +
  geom_col(width = 0.85, show.legend = FALSE, fill = "#619b8a") +
  scale_y_continuous(limits = c(-8, 43),
                     breaks = c(0, 10, 20, 30, 40),
                     expand = c(0, 0)) +
  coord_polar(theta = "x", start = pi / 1) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = c("grey85", "grey85", "grey85","grey85","grey85", NA), linewidth = 0.2),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_text(size=5, family = "sans", face = "bold"),
        axis.title.x = element_text(size=3, family = "sans", face = "bold"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=3, family = "sans", face = "bold", color = "black"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold", vjust = -6),
        axis.ticks = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0,0.1,0,0.1), "cm")) +
  geom_col(width = 0.85, show.legend = FALSE, fill = "#619b8a") +
  ggplot2::annotate("text",
                    x = pi*2.2,
                    y = c(10, 20, 30, 40),
                    label = c("10", "20", "30", "40"),
                    color = "black",
                    family = "sans",
                    size = 1)

#Combine
FIG <-
  grid.arrange(C,D,E,
               ncol=3,
               nrow=1)


#Save the figures
ggsave(FIG,
       width = 4.75, height = 1.5, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/main_text_panels/figure_1cde.png")

