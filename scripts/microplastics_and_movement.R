# This script processes runs the analyses investigating the relationship
# between the observed concentrations of microplastics in blood
# and the animals' patterns of movement and habitat use.
# For each covariate, a statistical anaysis is run, and a figure is generated

# Written by Michael Noonan

# Last updated: October 30 2024

library(mgcv)
library(ggplot2)
library(gridExtra)

#Import the MP concentration dataset
mp_data <- read.csv("data/mp_data/All_general.csv")
mp_data$sex <- as.factor(mp_data$sex)
mp_data$species <- as.factor(mp_data$species)
mp_data[which(mp_data$name == "Alvinho/Alvaro"),"name"] <- "Alvinho_Alvaro"

#Import the estimated movement metrics for the 3 species
tapir_data <- read.csv("data/movement_data/tapir_movement.csv")
anteater_data <- read.csv("data/movement_data/anteater_movement.csv", skip = 2)
armadillo_data <- read.csv("data/movement_data/armadillo_movement.csv")
move_data <- dplyr::bind_rows(tapir_data, armadillo_data)

#Import the estimated land usefor the 3 species
tapir_data <- read.csv("data/movement_data/tapir_land_use.csv")
#anteater_data <- read.csv("data/movement_data/anteater_land_use.csv")
armadillo_data <- read.csv("data/movement_data/armadillo_land_use.csv")
land_use <- dplyr::bind_rows(tapir_data, armadillo_data)
land_use[is.na(land_use)] <- 0

#Merge with the microplastics information
mp_data <- merge(x = mp_data, y = move_data, by.x = c("name", "species"), by.y = c("ID", "binomial"), all.x = TRUE)
mp_data <- merge(x = mp_data, y = land_use, by.x = c("name", "species"), by.y = c("ID", "binomial"), all.x = TRUE)



#-------------------------------------------------------------
# Correlation with home-range size
#-------------------------------------------------------------

#Test for a correlation between home range size and blood MP concentration
fit <- gam(mp_ml ~ hr + s(species, bs = "re") + s(hr, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)



a <- 
  ggplot(data = mp_data, aes(x = hr*1e-6, y = mp_ml)) +
  ggtitle("A)") +
  geom_smooth(method = "gam", formula = y ~ x, method.args = list(family = tw(link = "log")), col = "black", fill = "grey80", size = 0.2, linetype = "dashed") +
  geom_point(aes(col = species),size = 0.4) +
  scale_colour_manual(values = c("#619b8a", "#bb3e03", "#005f73"),
                      name = "",
                      labels = c("Myrmecophaga tridactyla", "Priodontes maximus", "Tapirus terrestris")) +
  ylab("MP concentration (particles/mL)") +
  xlab(expression(bold(Home-range~size~(km^2))))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=5, family = "sans", face = "bold"),
        axis.title.x = element_text(size=5, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 8, family = "sans", face = "bold"),
        legend.text  = element_text(size=5, family = "sans", face = "bold"),
        axis.ticks.length=unit(0.08, "cm"),
        axis.ticks = element_line(size = 0.3),
        legend.position = c(0.7,0.95),
        legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) 




#-------------------------------------------------------------
# Correlation with diffusion rate
#-------------------------------------------------------------

#Test for a correlation between home range size and blood MP concentration
fit <- gam(mp_ml ~ diffusion + s(species, bs = "re") + s(diffusion, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)



b <- 
  ggplot(data = mp_data, aes(x = diffusion, y = mp_ml)) +
  ggtitle("B)") +
  geom_smooth(method = "gam", formula = y ~ x, method.args = list(family = tw(link = "log")), col = "black", fill = "grey80", size = 0.2, linetype = "dashed") +
  geom_point(aes(col = species),size = 0.4) +
  scale_colour_manual(values = c("#619b8a", "#bb3e03", "#005f73"), name = "") +
  ylab("MP concentration (particles/mL)") +
  xlab(expression(bold(Diffusion~rate~(m^2~sec^-1))))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=5, family = "sans", face = "bold"),
        axis.title.x = element_text(size=5, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 8, family = "sans", face = "bold"),
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
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) 




#-------------------------------------------------------------
# Correlation with human footprint index
#-------------------------------------------------------------

#Test for a correlation between mean HFI and blood MP concentration
fit <- gam(mp_ml ~ mean_HFI + s(species, bs = "re") + s(mean_HFI, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)


c <- 
  ggplot(data = mp_data, aes(x = mean_HFI, y = mp_ml)) +
  ggtitle("C)") +
  geom_smooth(method = "gam", formula = y ~ x, method.args = list(family = tw(link = "log")), col = "black", fill = "grey80", size = 0.2, linetype = "dashed") +
  geom_point(aes(col = species),size = 0.4) +
  scale_colour_manual(values = c("#619b8a", "#bb3e03", "#005f73"), name = "") +
  ylab("MP concentration (particles/mL)") +
  xlab("Mean Human Footprint Index") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=5, family = "sans", face = "bold"),
        axis.title.x = element_text(size=5, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 8, family = "sans", face = "bold"),
        axis.ticks.length=unit(0.08, "cm"),
        axis.ticks = element_line(size = 0.3),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm"))




#-------------------------------------------------------------
# Correlation with agricultural land
#-------------------------------------------------------------


#Test for a correlation between the ammount of agg land in the HR and blood MP concentration
fit <- gam(mp_ml ~ Agriculture + s(species, bs = "re") + s(Agriculture, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)

d <- 
  ggplot(data = mp_data, aes(x = Agriculture, y = mp_ml)) +
  ggtitle("D)") +
  geom_smooth(method = "gam", formula = y ~ x, method.args = list(family = tw(link = "log")), col = "black", fill = "grey80", size = 0.2, linetype = "dashed") +
  geom_point(aes(col = species),size = 0.4) +
  scale_colour_manual(values = c("#619b8a", "#bb3e03", "#005f73"), name = "") +
  ylab("MP concentration (particles/mL)") +
  xlab("Proportion of home range in agricultural land") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=5, family = "sans", face = "bold"),
        axis.title.x = element_text(size=5, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 8, family = "sans", face = "bold"),
        axis.ticks.length=unit(0.08, "cm"),
        axis.ticks = element_line(size = 0.3),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm"))




#-------------------------------------------------------------
# Correlation with water and wetlands
#-------------------------------------------------------------


#Test for a correlation between mean HFI and blood MP concentration
fit <- gam(mp_ml ~ Water + s(species, bs = "re") + s(Water, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)


e <- 
ggplot(data = mp_data, aes(x = Water, y = mp_ml)) +
  ggtitle("E)") +
  geom_smooth(method = "gam", formula = y ~ x, method.args = list(family = tw(link = "log")), col = "black", fill = "grey80", size = 0.2, linetype = "solid") +
  geom_point(aes(col = species),size = 0.4) +
  scale_colour_manual(values = c("#619b8a", "#bb3e03", "#005f73"), name = "") +
  ylab("MP concentration (particles/mL)") +
  xlab("Proportion of home range in water and wetlands") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=5, family = "sans", face = "bold"),
        axis.title.x = element_text(size=5, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 8, family = "sans", face = "bold"),
        axis.ticks.length=unit(0.08, "cm"),
        axis.ticks = element_line(size = 0.3),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm"))



#-------------------------------------------------------------
# Correlation with water and wetlands
#-------------------------------------------------------------


#Test for a correlation between mean HFI and blood MP concentration
fit <- gam(mp_ml ~ Native_forest + s(species, bs = "re") + s(Native_forest, species, bs = "re"),
           family = tw(link = "log"),
           data = mp_data,
           method = "REML")

summary(fit)


f <- 
  ggplot(data = mp_data, aes(x = Native_forest, y = mp_ml)) +
  ggtitle("F)") +
  geom_smooth(method = "gam", formula = y ~ x, method.args = list(family = tw(link = "log")), col = "black", fill = "grey80", size = 0.2, linetype = "dashed") +
  geom_point(aes(col = species),size = 0.4) +
  scale_colour_manual(values = c("#619b8a", "#bb3e03", "#005f73"), name = "") +
  ylab("MP concentration (particles/mL)") +
  xlab("Proportion of home range in natural forest") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=5, family = "sans", face = "bold"),
        axis.title.x = element_text(size=5, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 8, family = "sans", face = "bold"),
        axis.ticks.length=unit(0.08, "cm"),
        axis.ticks = element_line(size = 0.3),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm"))


FIG <-
  grid.arrange(a,b,c,d,e,f,
               ncol=2,
               nrow=3)


#Save the figures
ggsave(FIG,
       width = 4.75, height = 5.1, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/Movement_Trends.png")

