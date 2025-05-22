# This script generates figure 3 in the main text that
# visualises differences in polymer concentrations between the three biomes
# Note: the indiviudal panels are generated in this script
# but need to be assembled outside of R

# Written by Michael Noonan

#Load in any necessary packages
library(ggplot2)
library(fmsb)

#Import the MP datasets
source("scripts/data_import.R")


#---------------------------------------------------------------------
# Figure 3A - heatmap of polymer abundances
#---------------------------------------------------------------------


#Reorder the sample ID column by biome for cleaner plotting
ORDER <- mp_data[order(mp_data$biome), "sample"]
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
  geom_hline(yintercept = 14.5, linetype = "dashed", col = "grey70", linewidth = 0.3) +
  geom_hline(yintercept = 36.5, linetype = "dashed", col = "grey70", linewidth = 0.3) +
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
       file="figures/figure_3a.png")


#---------------------------------------------------------------------
# Figure 3B -D - Radar plots of polymer abundances in each biome
#---------------------------------------------------------------------

#Mean polymer concentrations in each biome
mean_polymers <- aggregate(concentration ~ biome + polymer,
                           FUN = "mean",
                           data = mp_data_long)

#Convert to wide format
mean_polymers <- reshape(data = mean_polymers,
                         idvar = "biome",
                         timevar = "polymer",
                         v.names = "concentration",
                         direction = "wide")

#Adjust the names
names(mean_polymers)[-1] <- sub("concentration.","",names(mean_polymers)[-1])
names(mean_polymers)[11] <- "PVC"

#Define the min and max (specific requirement of the fmsb package)
mean_polymers[4,] <- 60 #Max
mean_polymers[5,] <- 0 #Min

row.names(mean_polymers) <- c("Amazon",
                              "Cerrado",
                              "Pantanal",
                              "Max",
                              "Min")

#Reorder the rows to match the requirement of the fmsb package
mean_polymers <- mean_polymers[c(4:5,1:3),]


# Generate and save the figures
png(filename = "figures/figure_3bcd.png",
    width = 2.375, height = 5, units = "in",
    bg = "transparent",
    res = 600)

par(mfrow = c(3,1),
    mar = c(0.2,0.1,0.2,0.1),
    font.axis = 2,
    font = 2)

# Create the radar chart for the amazon
radarchart(mean_polymers[c(1:2,3),-1],
           axistype = 1,
           pty = NA,
           pcol = "#007200",
           pfcol = adjustcolor("#007200", alpha.f = 0.3),
           plwd = 1,
           cglcol = "grey",
           cglty = 1,
           axislabcol = "grey",
           caxislabels = c(0, 15, 30, 45, 60),
           calcex = 0.5,
           cglwd = 0.4,
           vlcex = 0.48)
title(main = "B", cex.main = 0.8, font.main= 2, col.main= "black", adj = 0, line = -1)

# Create the radar chart for the armadillos
radarchart(mean_polymers[c(1:2,4),-1],
           axistype = 1,
           pty = NA,
           pcol = "#e9c46a",
           pfcol = adjustcolor("#e9c46a", alpha.f = 0.3),
           plwd = 1,
           cglcol = "grey",
           cglty = 1,
           axislabcol = "grey",
           caxislabels = c(0, 15, 30, 45, 60),
           calcex = 0.5,
           cglwd = 0.4,
           vlcex = 0.48)
title(main = "C", cex.main = 0.8, font.main= 2, col.main= "black", adj = 0, line = -1)


# Create the radar charts for the anteaters
radarchart(mean_polymers[c(1:2,5),-1],
           axistype = 1,
           pty = NA,
           pcol = "#0a9396",
           pfcol = adjustcolor("#0a9396", alpha.f = 0.3),
           plwd = 1,
           cglcol = "grey",
           cglty = 1,
           axislabcol = "grey",
           caxislabels = c(0, 15, 30, 45, 60),
           calcex = 0.5,
           cglwd = 0.4,
           vlcex = 0.48)
title(main = "D",
      cex.main = 0.8, font.main= 2, col.main= "black", adj = 0, line = -1)

dev.off()