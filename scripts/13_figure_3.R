# This script generates figure 3 in the main text that
# visualises differences in polymer concentrations between the three biomes
# Note: the indiviudal panels are generated in this script
# but need to be assembled outside of R

# Written by Michael Noonan

#Load in any necessary packages
library(ggplot2)
library(gridExtra)

#Import the MP datasets
source("scripts/04_microplastics_data_import.R")


#---------------------------------------------------------------------
# Figure 3A - heatmap of polymer abundances
#---------------------------------------------------------------------


#Reorder the sample ID column by biome for cleaner plotting
ORDER <- mp_data[order(mp_data$biome), "sample"]
mp_data_long$sample <- factor(mp_data_long$sample, levels = ORDER, ordered = TRUE)


# Heatmap of the polymer abundances
A <- 
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


#---------------------------------------------------------------------
# Figure 3B -D - Radar barplots of polymer abundances in each study site
#---------------------------------------------------------------------


mean_polymers <- aggregate(concentration ~ biome + polymer,
                           FUN = "mean", data = mp_data_long)


B <- 
  ggplot(mean_polymers[mean_polymers$biome == "Amazon",],
         aes(x = polymer, y = concentration)) +
  ggtitle("B") +
  geom_segment(data = data.frame(polymer = unique(mean_polymers$polymer)),
               aes(x = polymer, xend = polymer, y = -8, yend = 30),
               color = "grey85", linewidth = 0.2, inherit.aes = FALSE) +
  geom_col(width = 0.85, show.legend = FALSE, fill = "#007200") +
  scale_y_continuous(limits = c(-8, 33),
                     breaks = c(0, 10, 20, 30, 40),
                     expand = c(0, 0)) +
  coord_polar(theta = "x", start = pi / 1) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = c("grey85", "grey85", "grey85","grey85", NA), linewidth = 0.2),
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


C <- 
  ggplot(mean_polymers[mean_polymers$biome == "Cerrado",],
         aes(x = polymer, y = concentration)) +
  ggtitle("C") +
  geom_segment(data = data.frame(polymer = unique(mean_polymers$polymer)),
               aes(x = polymer, xend = polymer, y = -8, yend = 30),
               color = "grey85", linewidth = 0.2, inherit.aes = FALSE) +
  geom_col(width = 0.85, show.legend = FALSE, fill = "#e9c46a") +
  scale_y_continuous(limits = c(-8, 33),
                     breaks = c(0, 10, 20, 30, 40),
                     expand = c(0, 0)) +
  coord_polar(theta = "x", start = pi / 1) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = c("grey85", "grey85", "grey85","grey85", NA), linewidth = 0.2),
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
  ggplot(mean_polymers[mean_polymers$biome == "Pantanal",],
         aes(x = polymer, y = concentration)) +
  ggtitle("D") +
  geom_segment(data = data.frame(polymer = unique(mean_polymers$polymer)),
               aes(x = polymer, xend = polymer, y = -8, yend = 30),
               color = "grey85", linewidth = 0.2, inherit.aes = FALSE) +
  geom_col(width = 0.85, show.legend = FALSE, fill = "#0a9396") +
  scale_y_continuous(limits = c(-8, 33),
                     breaks = c(0, 10, 20, 30, 40),
                     expand = c(0, 0)) +
  coord_polar(theta = "x", start = pi / 1) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = c("grey85", "grey85", "grey85","grey85", NA), linewidth = 0.2),
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

#Combine
FIG <-
  grid.arrange(B,C,D,
               ncol=1,
               nrow=3)


#Combine
FIG <-
  grid.arrange(A, FIG,
               ncol=2,
               nrow=1)

#Save the figures
ggsave(FIG,
       width = 2.375*2, height = 5, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/supplementary/extended_data_figure_1.png")
