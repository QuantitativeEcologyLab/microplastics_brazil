# This script generates figure 1 in the main text that
# visualises differences in the study sites and general patterns
# the plastic concentrations

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

# Import the biome boundaries
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
                  alpha = 0.7) +
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
        plot.margin = unit(c(0.2,3,-0.3,-1.5), "cm")) +
  
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



#---------------------------------------------------------------------
# Figure 1B - Barplot of polymer abundances
#---------------------------------------------------------------------

# Define the colors for each of the polymers
polymer_colors <- c(
  "ABS" = "#94D2BD",
  "Cellulose" = "#E9D8A6",
  "PET" = "grey35",
  "Polyamide" = "#CA6702",
  "Polyester" = "#9B2226",
  "Polyethylene" = "#6d597a",
  "Polypropylene" = "#005F73",
  "Polystyrene" = "#0A9396",
  "PVC" = "#8d96a3",
  "Rubber" = "green4",
  "Silicon" = "#355070",
  "Teflon" = "coral2",
  "Other" = "#EE9B00"
)

mp_data_long$polymer <- factor(mp_data_long$polymer, levels = names(polymer_colors))

# First, reorder the 'Sample' column by its numeric component
ORDER <- mp_data[order(mp_data$species, mp_data$mp_ml), "sample"]
mp_data_long$sample <- factor(mp_data_long$sample, levels = ORDER, ordered = TRUE)

# Create the bar plot with your custom colors
b <- 
  ggplot(data = mp_data_long, 
         aes(x = sample,
             y = concentration,
             fill = polymer)) +
  ggtitle("B") +
  geom_vline(xintercept = 21.5, linetype = "dashed", col = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 26.5, linetype = "dashed", col = "grey70", linewidth = 0.3) +
  geom_bar(stat = "identity", alpha = 0.9) +
  scale_fill_manual(values = polymer_colors) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(size = 0.2),
        axis.line.y = element_line(size = 0.2),
        axis.ticks.x = element_line(size = 0.2),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 5, family = "sans", face = "bold"),
        axis.title.y = element_blank(),#element_text(size = 5, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 4, family = "sans"),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold"),
        legend.position = c(0.8, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 3, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(0.08, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.spacing.y = unit(-0.1, "cm"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.4,0.1,-0,-5.5), "cm")
  ) +
  scale_x_discrete(breaks = ORDER) +
  scale_y_continuous(expand = c(0, 0)) +
  # scale_y_log10(breaks = c(0.01,0.1,1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000,10000000000,100000000000),
  #               labels = c(0.01,0.1,1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000,10000000000,100000000000),
  #               expand = c(0, 0)) +
  # annotation_logticks(sides="b",
  #                     outside = TRUE,
  #                     size = 0.2,
  #                     short = unit(0.05, "cm"),
  #                     mid = unit(0.05, "cm"),
  #                     long = unit(0.1, "cm")) +
  xlab(expression(bold(Sample))) +
  ylab(expression(bold(Microplastic~concentration~(particles/mL)))) +
  coord_flip() +
  add_phylopic(anteater_pic,
               x = 10,
               y = 300,
               ysize = 2, alpha = 1, fill = "#619b8a") +
  add_phylopic(armadillo_pic,
               x = 24,
               y = 300,
               ysize = 2, alpha = 1, fill = "#bb3e03") +
  add_phylopic(tapir_pic,
               x = 38,
               y = 300,
               ysize = 2, alpha = 1, fill = "#005f73",
               horizontal = TRUE)


top <-
  grid.arrange(a,b,
               ncol=2,
               nrow=1,
               widths=c(7,0.25))


#Save the figures
ggsave(top,
       width = 4.75, height = 2.5, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/figure_1_top.png")


#---------------------------------------------------------------------
# Figure 1C & D - Histograms of size distributions
#---------------------------------------------------------------------



#Load in the data
metadata <- read.csv("data/mp_data/All_general.csv")
sizes <- read.csv("data/mp_data/All_size.csv")
sizes <- merge(x = metadata, y = sizes, by.x = "sample", by.y = "Sample")


c <- 
  ggplot() +
  ggtitle("C") +
  geom_histogram(data = sizes, aes(Length, fill = species),
                 alpha = 0.8,
                 bins = 60,
                 col = "black",
                 linewidth = 0.05) +
  scale_fill_manual(values = c("#619b8a", "#bb3e03", "#005f73")) +
  scale_x_log10(expand = c(0,0.1)) +
  scale_y_continuous(limits = c(0,230), expand = c(0,.1)) +
  ylab("Number of particles") +
  xlab(bquote(bold('Particle length '(µm)))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=5, family = "sans", face = "bold"),
        axis.title.x = element_text(size=5, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=5, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  guides(fill = guide_legend(byrow = TRUE)) +
    add_phylopic(anteater_pic,
                 x = 150,
                 y = 200,
                 ysize = 14, alpha = 1, fill = "#619b8a") +
    add_phylopic(armadillo_pic,
                 x = 150,
                 y = 170,
                 ysize = 14, alpha = 1, fill = "#bb3e03") +
    add_phylopic(tapir_pic,
                 x = 150,
                 y = 140,
                 ysize = 14, alpha = 1, fill = "#005f73",
                 horizontal = TRUE)



d <- 
  ggplot() +
  ggtitle("D") +
  geom_histogram(data = sizes, aes(Width, fill = species),
                 alpha = 0.8,
                 bins = 50,
                 col = "black",
                 linewidth = 0.05) +
  scale_fill_manual(values = c("#619b8a", "#bb3e03", "#005f73")) +
  scale_x_log10(expand = c(0,0.01)) +
  scale_y_continuous(limits = c(0,280), expand = c(0,.1)) +
  ylab("Number of particles") +
  xlab(bquote(bold('Particle width '(µm)))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=5, family = "sans", face = "bold"),
        axis.title.x = element_text(size=5, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=5, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  guides(fill = guide_legend(byrow = TRUE))




bot <-
  grid.arrange(c,d,
               ncol=2,
               nrow=1)


#Save the figures
ggsave(bot,
       width = 4.75, height = 1.5, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/figure_1_bottom.png")

