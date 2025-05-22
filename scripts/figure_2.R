# This script generates figure 2 in the main text that
# visualises the correlation with HFI
#Note: The insets in figure 2a are generated in the script "habitat_figures.R"

# Written by Michael Noonan

#Load in any necessary packages
library(mgcv)
library(ggplot2)
library(fmsb)
library(rphylopic)

#Import the MP datasets
source("scripts/data_import.R")


#Get the animal silhouettes 
tapir_pic <- get_phylopic("7950e979-6738-45b3-a7c6-c573ef5559d1")
anteater_pic <- get_phylopic("d52b48fc-be52-46a1-94b7-ac7790b4730c")
armadillo_pic <- get_phylopic("5d59b5ce-c1dd-40f6-b295-8d2629b9775e")



#-------------------------------------------------------------
# Figure 2B - Correlation with mean human footprint index



b <- 
  ggplot(data = mp_data, aes(x = mean_HFI, y = mp_ml)) +
  ggtitle("B") +
  geom_smooth(method = "gam",
              formula = y ~ x,
              method.args = list(family = tw(link = "log")),
              col = "black",
              fill = "grey80",
              size = 0.2,
              linetype = "solid") +
  geom_point(aes(col = species),size = 0.4) +
  scale_colour_manual(values = c("#619b8a", "#bb3e03", "#005f73"), name = "") +
  ylab("Microplastics (particles/mL)") +
  xlab("Mean human footprint index in home range") +
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
  #coord_cartesian(ylim = c(5,170)) +
  add_phylopic(anteater_pic,
               x = 26,
               y = 160,
               ysize = 12, alpha = 1, fill = "#619b8a") +
  add_phylopic(armadillo_pic,
               x = 22,
               y = 160,
               ysize = 12, alpha = 1, fill = "#bb3e03") +
  add_phylopic(tapir_pic,
               x = 18,
               y = 160,
               ysize = 12, alpha = 1, fill = "#005f73",
               horizontal = TRUE)


#-------------------------------------------------------------
# Figure 2C - Correlation with maximum human footprint


c <- 
  ggplot(data = mp_data, aes(x = max_HFI, y = mp_ml)) +
  ggtitle("C") +
  geom_smooth(method = "gam", formula = y ~ x, method.args = list(family = tw(link = "log")), col = "black", fill = "grey80", size = 0.2, linetype = "solid") +
  geom_point(aes(col = species),size = 0.4) +
  scale_colour_manual(values = c("#619b8a", "#bb3e03", "#005f73"), name = "") +
  ylab("Microplastics (particles/mL)") +
  xlab("Maximum human footprint index in home range") +
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





FIG <-
  grid.arrange(b,c,
               ncol=2,
               nrow=1)


#Save the figures
ggsave(FIG,
       width = 4.75, height = 1.5, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/figure_2bc.png")
