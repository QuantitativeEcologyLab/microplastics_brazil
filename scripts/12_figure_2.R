# This script generates figure 2 in the main text that
# visualises the correlation with HFI and the difference in the human footprint index
#Note: The insets in figure 2a are generated in the script "habitat_figures.R"

# Written by Michael Noonan

#Load in any necessary packages
library(mgcv)
library(ggplot2)
library(fmsb)
library(rphylopic)
library(tidyterra)
library(gridExtra)

#Import the MP datasets
source("scripts/04_microplastics_data_import.R")

#Exclude the giant armadillos for the main text analyses
mp_data <- mp_data[mp_data$species != "Priodontes_maximus",]

#Get the animal silhouettes from the phylopic package
tapir_pic <- get_phylopic("7950e979-6738-45b3-a7c6-c573ef5559d1")
anteater_pic <- get_phylopic("d52b48fc-be52-46a1-94b7-ac7790b4730c")
armadillo_pic <- get_phylopic("5d59b5ce-c1dd-40f6-b295-8d2629b9775e")


#-------------------------------------------------------------
# Figure 2A - Map of difference in human footprint index
#-------------------------------------------------------------

#Import the 2001 and 2020 HII rasters
#Note: Files are too large to store on github, but are available here: https://wcshumanfootprint.org/data-access
hii_2001 <- rast("~/Dropbox/UBC/Projects/microplastics_brazil/data/environmental_data/hii_2001-01-01.tif")/100
hii_2020 <- rast("~/Dropbox/UBC/Projects/microplastics_brazil/data/environmental_data/hii_2020-01-01.tif")/100
#Create a difference raster (2020 minus 2001)
hii_diff <- hii_2020 - hii_2001

#Get the 0.01 and 0.99 quantiles
hii_diffs <- na.omit(values(diff_br))
quants <-  quantile(hii_diffs, c(0.01, 0.99))

#Clamp the difference raster to the 0.01 and 0.99 quantiles (to help the visualisation)
hii_diff[hii_diff < quants[1]] <- quants[1]
hii_diff[hii_diff > quants[2]] <- quants[2]

#Get the boundaries of Brasil
world_sf <- st_as_sf(rworldmap::getMap(resolution = "low"))
brasil_sf <- subset(world_sf, SOVEREIGNT == "Brazil")
brasil <- vect(lwgeom::st_transform_proj(brasil_sf, crs = crs(hii_2001)))


#Generate the figure
a <-
  ggplot() +
  geom_spatraster(data = hii_diff, maxcell = 5e+07,
                  alpha = 1) +
  scale_fill_scico(name = "HFI Difference",
                   palette = "vik",
                   midpoint = 0,
                   breaks = c(-4,0,4,8,12),
                   na.value = NA) +
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA, barwidth = 4,
                               barheight = 0.2, direction = "horizontal")) +
  #geom_sf(data = brasil, size = 0.2, fill = "transparent", linewidth = 0.1, col = "black") +
  coord_sf(ylim = c(-50, 30), xlim = c(-95, -10)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.position = "inside",
        legend.position.inside = c(0.78,0.12),
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
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) 


#Save the figures
ggsave(a,
       width = 4.75, height = 3, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/figure_2a.png")


#-------------------------------------------------------------
# Figure 2b - Correlation with maximum human footprint index
#-------------------------------------------------------------


b <- 
  ggplot(data = mp_data, aes(x = max_HFI, y = mp_ml)) +
  ggtitle("b") +
  geom_smooth(method = "gam",
              formula = y ~ x,
              method.args = list(family = tw(link = "log")),
              col = "black",
              fill = "grey80",
              linewidth = 0.2,
              linetype = "solid") +
  geom_point(aes(col = species),size = 0.4) +
  scale_colour_manual(values = c("#619b8a", "#005f73"), name = "") +
  ylab("Microplastics (particles/mL)") +
  xlab("Maximum HFI exposure") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=5, family = "sans", face = "bold"),
        axis.title.x = element_text(size=5, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold", vjust = -2),
        axis.ticks.length=unit(0.08, "cm"),
        axis.ticks = element_line(size = 0.3),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0,0.1,0,0.1), "cm")) +
  coord_cartesian(ylim = c(5,170))

#-------------------------------------------------------------
# Figure 2c - Correlation with distance to development
#-------------------------------------------------------------

c <- 
  ggplot(data = mp_data, aes(x = mean_dist_development, y = mp_ml)) +
  ggtitle("c") +
  geom_smooth(method = "gam", formula = y ~ x, method.args = list(family = tw(link = "log")), col = "black", fill = "grey80", linewidth = 0.2, linetype = "solid") +
  geom_point(aes(col = species),size = 0.4) +
  scale_colour_manual(values = c("#619b8a", "#005f73"), name = "") +
  ylab("Microplastics (particles/mL)") +
  xlab("Mean distance to development (km)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=5, family = "sans", face = "bold"),
        axis.title.x = element_text(size=5, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold", vjust = -2),
        axis.ticks.length=unit(0.08, "cm"),
        axis.ticks = element_line(size = 0.3),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0,0.1,0,0.1), "cm")) +
  coord_cartesian(ylim = c(5,170)) +
  # add_phylopic(armadillo_pic,
  #              x = 50,
  #              y = 160,
  #              height = 16, alpha = 1, fill = "#bb3e03") +
  add_phylopic(anteater_pic,
               x = 55,
               y = 160,
               height = 16, alpha = 1, fill = "#619b8a") +
  add_phylopic(tapir_pic,
               x = 65,
               y = 160,
               height = 16, alpha = 1, fill = "#005f73",
               horizontal = TRUE)




#Combine
FIG <-
  grid.arrange(b,c,
               ncol=2,
               nrow=1)


#Save the figures
ggsave(FIG,
       width = 4.75, height = 1.5, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/main_text_panels/figure_2bc.png")


#Save the underlying data
write.csv(na.omit(mp_data[,c("species","mp_ml","max_HFI")]),
          file = "data/figures/figure_2b.csv",
          row.names = F)

write.csv(na.omit(mp_data[,c("species","mp_ml","mean_dist_development")]),
          file = "data/figures/figure_2c.csv",
          row.names = F)
