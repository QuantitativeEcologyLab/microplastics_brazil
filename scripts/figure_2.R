# This script generates figure 2 in the main text that
# visualises differences in polymer concentrations between the three species

# Written by Michael Noonan

#Load in any necessary packages
library(mgcv)
library(ggplot2)
library(fmsb)

#Import the MP datasets
source("scripts/data_import.R")


#---------------------------------------------------------------------
# Figure 2A - heatmap of polymer abundances
#---------------------------------------------------------------------


#Reorder the sample ID column by species for cleaner plotting
ORDER <- mp_data[order(mp_data$species), "sample"]
mp_data_long$sample <- factor(mp_data_long$sample, levels = ORDER, ordered = TRUE)


# Heatmap of the polymer abundances
a <- 
  ggplot(mp_data_long, aes(polymer, sample, fill= log(concentration+1))) + 
  ggtitle("A") +
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


ggsave(a,
       width = 2.375, height = 5, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/figure_2_left.png")


#---------------------------------------------------------------------
# Figure 2B -D - Radar plots of polymer abundances in each species
#---------------------------------------------------------------------

#Mean polymer concentrations in each species
mean_polymers <- aggregate(concentration ~ species + polymer,
                           FUN = "mean",
                           data = mp_data_long)

#Convert to wide format
mean_polymers <- reshape(data = mean_polymers,
                         idvar = "species",
                         timevar = "polymer",
                         v.names = "concentration",
                         direction = "wide")

#Adjust the names
names(mean_polymers)[-1] <- sub("concentration.","",names(mean_polymers)[-1])


#Define the min and max (specific requirement of the fmsb package)
mean_polymers[4,] <- 40 #Max
mean_polymers[5,] <- 0 #Min

row.names(mean_polymers) <- c("Myrmecophaga_tridactyla",
                              "Priodontes_maximus",
                              "Tapirus_terrestris",
                              "Max",
                              "Min")

#Reorder the rows to match the requirement of the fmsb package
mean_polymers <- mean_polymers[c(4:5,1:3),]


# Generate and save the figures
png(filename = "figures/figure_2_right.png",
    width = 2.375, height = 5, units = "in",
    bg = "transparent",
    res = 600)

par(mfrow = c(3,1),
    mar = c(0.2,0.1,0.2,0.1),
    font.axis = 2,
    font = 2)

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
title(main = "B", cex.main = 0.8, font.main= 2, col.main= "black", adj = 0, line = -1)

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
title(main = "C", cex.main = 0.8, font.main= 2, col.main= "black", adj = 0, line = -1)


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
title(main = "D",
      cex.main = 0.8, font.main= 2, col.main= "black", adj = 0, line = -1)

dev.off()