# This script generates figure 5 in the main text that
# visualises relationships in polymer concentrations and movement patterns

# Written by Michael Noonan

#Load in any necessary packages
library(ggplot2)
library(gridExtra)
library(rphylopic)

#Import the MP datasets
source("scripts/data_import.R")


#Get the animal silhouettes 
tapir_pic <- get_phylopic("7950e979-6738-45b3-a7c6-c573ef5559d1")
anteater_pic <- get_phylopic("d52b48fc-be52-46a1-94b7-ac7790b4730c")
armadillo_pic <- get_phylopic("5d59b5ce-c1dd-40f6-b295-8d2629b9775e")

#-------------------------------------------------------------
# Figure 5A - Relationship with home-range size
#-------------------------------------------------------------


a <- 
  ggplot(data = mp_data[mp_data$name != "Juliana",], aes(x = hr, y = mp_ml)) +
  ggtitle("A") +
  geom_smooth(method = "gam",
              formula = y ~ x,
              method.args = list(family = tw(link = "log")),
              col = "black",
              fill = "grey80",
              linewidth = 0.2,
              linetype = "dashed") +
  geom_point(aes(col = species), size = 0.4) +
  scale_colour_manual(values = c("#619b8a", "#bb3e03", "#005f73"),
                      name = "",
                      labels = c("Myrmecophaga tridactyla", "Priodontes maximus", "Tapirus terrestris")) +
  ylab("Microplastics (particles/mL)") +
  xlab(expression(bold(Home~range~size~(km^2))))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=5, family = "sans", face = "bold"),
        axis.title.x = element_text(size=5, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold", vjust = -2),
        legend.text  = element_text(size=5, family = "sans", face = "bold"),
        axis.ticks.length=unit(0.08, "cm"),
        axis.ticks = element_line(size = 0.3),
        legend.position = "none",
        legend.position.inside  = c(0.7,0.95),
        legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0,0.1,0,0.1), "cm")) +
  coord_cartesian(ylim = c(5,170)) +
  add_phylopic(anteater_pic,
               x = 28,
               y = 160,
               ysize = 12, alpha = 1, fill = "#619b8a") +
  add_phylopic(armadillo_pic,
               x = 24,
               y = 160,
               ysize = 12, alpha = 1, fill = "#bb3e03") +
  add_phylopic(tapir_pic,
               x = 20,
               y = 160,
               ysize = 12, alpha = 1, fill = "#005f73",
               horizontal = TRUE)






#-------------------------------------------------------------
# Figre 5B - Relationship with diffusion rate
#-------------------------------------------------------------


b <- 
  ggplot(data = mp_data, aes(x = diffusion, y = mp_ml)) +
  ggtitle("B") +
  geom_smooth(method = "gam",
              formula = y ~ x,
              method.args = list(family = tw(link = "log")),
              col = "black",
              fill = "grey80",
              size = 0.2,
              linetype = "dashed") +
  geom_point(aes(col = species),size = 0.4) +
  scale_colour_manual(values = c("#619b8a", "#bb3e03", "#005f73"), name = "") +
  ylab("Microplastics (particles/mL)") +
  xlab(expression(bold(Diffusion~rate~(m^2~sec^-1))))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=5, family = "sans", face = "bold"),
        axis.title.x = element_text(size=5, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold", vjust = -2),
        legend.text  = element_text(size=5, family = "sans", face = "bold"),
        axis.ticks.length=unit(0.08, "cm"),
        axis.ticks = element_line(size = 0.3),
        legend.position = "none",
        legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0,0.1,0,0.1), "cm")) +
  coord_cartesian(ylim = c(5,170))



#-------------------------------------------------------------
# Figure 5C - Correlation with mean human footprint index
#-------------------------------------------------------------



c <- 
  ggplot(data = mp_data, aes(x = mean_HFI, y = mp_ml)) +
  ggtitle("C") +
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
  coord_cartesian(ylim = c(5,170))


#-------------------------------------------------------------
# Figure 5D - Correlation with max HFI
#-------------------------------------------------------------


d <- 
  ggplot(data = mp_data, aes(x = max_HFI, y = mp_ml)) +
  ggtitle("D") +
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



#-------------------------------------------------------------
# Figure 5E - Correlation with agricultural land
#-------------------------------------------------------------


e <- 
  ggplot(data = mp_data, aes(x = Agriculture, y = mp_ml)) +
  ggtitle("E") +
  geom_smooth(method = "gam", formula = y ~ x, method.args = list(family = tw(link = "log")), col = "black", fill = "grey80", size = 0.2, linetype = "dashed") +
  geom_point(aes(col = species),size = 0.4) +
  scale_colour_manual(values = c("#619b8a", "#bb3e03", "#005f73"), name = "") +
  ylab("Microplastics (particles/mL)") +
  xlab("Proportion of home range in agricultural land") +
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
# Figure 5F - Correlation with water and wetlands
#-------------------------------------------------------------



f <- 
  ggplot(data = mp_data, aes(x = Water, y = mp_ml)) +
  ggtitle("F") +
  geom_smooth(method = "gam", formula = y ~ x, method.args = list(family = tw(link = "log")), col = "black", fill = "grey80", size = 0.2, linetype = "solid") +
  geom_point(aes(col = species),size = 0.4) +
  scale_colour_manual(values = c("#619b8a", "#bb3e03", "#005f73"), name = "") +
  ylab("Microplastics (particles/mL)") +
  xlab("Proportion of water and wetlands in home range") +
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






g <- 
  ggplot(data = mp_data, aes(x = Native_forest, y = mp_ml)) +
  ggtitle("G") +
  geom_smooth(method = "gam", formula = y ~ x, method.args = list(family = tw(link = "log")), col = "black", fill = "grey80", size = 0.2, linetype = "dashed") +
  geom_point(aes(col = species),size = 0.4) +
  scale_colour_manual(values = c("#619b8a", "#bb3e03", "#005f73"), name = "") +
  ylab("Microplastics (particles/mL)") +
  xlab("Proportion of natural forest in home range") +
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


figure_5 <-
  grid.arrange(a,b,c,d,e,f,
               ncol=2,
               nrow=3)


#Save the figures
ggsave(figure_5,
       width = 4.75, height = 4, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/figure_5.png")
