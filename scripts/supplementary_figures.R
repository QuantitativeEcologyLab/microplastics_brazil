# This script generates the 5 supplementary figures associated with 
# the paper entitled "Human disturbance drives microplastic abundance in terrestrial wildlife"

#Note: A small amount of post processing was done to fig S3 to include animal silhouttes 

# Written by Michael Noonan

library(ggplot2)
library(gridExtra)
library(mgcv)

#Import the MP datasets
source("scripts/data_import.R")



#---------------------------------------------------------------------
# Figure S1 - Overview of the dataset with microplastic characterization across species
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# Figure S1a - Barplot of polymer abundances


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
top <- 
  ggplot(data = mp_data_long, 
         aes(x = sample,
             y = concentration,
             fill = polymer)) +
  ggtitle("A") +
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
        legend.text = element_text(size = 5, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(0.08, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.spacing.y = unit(-0.1, "cm"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")
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



#---------------------------------------------------------------------
# Figure S1B & C - Histograms of size distributions



#Load in the data
metadata <- read.csv("data/mp_data/All_general.csv")
sizes <- read.csv("data/mp_data/All_size.csv")
sizes <- merge(x = metadata, y = sizes, by.x = "sample", by.y = "Sample")


b <- 
  ggplot() +
  ggtitle("B") +
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



c <- 
  ggplot() +
  ggtitle("C") +
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
  grid.arrange(b,c,
               ncol=2,
               nrow=1)


FIG <-
  grid.arrange(top,bot,
               ncol=1,
               nrow=2,
               heights = c(1.4,.6))


#Save the figures
ggsave(FIG,
       width = 4.75, height = 4.5, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/figure_s1.png")



#-------------------------------------------------------------
# Figure S2 - Relationships between habitat use and microplastic concentrations 
#-------------------------------------------------------------

#-------------------------------------------------------------
# Figure S2a - Relationship with home-range size


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
# Figure S2B - Relationship with diffusion rate


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
# Figure S2C - Correlation with agricultural land
#-------------------------------------------------------------


c <- 
  ggplot(data = mp_data, aes(x = Agriculture, y = mp_ml)) +
  ggtitle("C") +
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
# Figure S2D - Correlation with water and wetlands
#-------------------------------------------------------------



d <- 
  ggplot(data = mp_data, aes(x = Water, y = mp_ml)) +
  ggtitle("D") +
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



#Compile the panels into a single figure
figure_s2 <-
  grid.arrange(a,b,c,d,
               ncol=2,
               nrow=2)


#Save the figure
ggsave(figure_s2,
       width = 4.75, height = 3, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/figure_s2.png")


#---------------------------------------------------------------------
# Figure S3 - Relationships between mp concentrations and biological variables
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# Figure S3a - correlation with age in tapirs

a <-
  ggplot(data = mp_data[which(mp_data$species == "Tapirus_terrestris"),],
         aes(x = age,
             y = mp_ml+1,
             col = sex,
             fill = sex)) +
  ggtitle("A") +
  geom_smooth(method = "gam",
              formula = y ~ x,
              method.args = list(family = tw(link = "log")),
              linewidth = 0.2,
              linetype = "solid",
              show.legend = FALSE) +
  geom_point(size = 0.4) +
  scale_colour_manual(values = c("#fca311", "#14213d"), labels = c("Female", "Male"), name = "") +
  scale_fill_manual(values = c("#fca311", "#14213d"), labels = c("Female", "Male"), name = "") +
  xlab("Age (years)") +
  ylab("Microplastic concentration (particles/mL)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=5, family = "sans", face = "bold"),
        axis.title.x = element_text(size=5, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold"),
        strip.text.x = element_text(size=6, family = "sans", face = "bold", color = "black"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=5, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_y_log10(breaks = c(1,10+1,100+1,1000+1,10000,100000,1000000,10000000,100000000,1000000000,10000000000,100000000000,1000000000000,10000000000000),
                labels = c(0,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000,10000000000,100000000000,1000000000000,10000000000000)) +
  annotation_logticks(#outside = TRUE,
    sides = "l",
    linewidth = 0.1,
    short = unit(0.05, "cm"),
    mid = unit(0.05, "cm"),
    long = unit(0.08, "cm"))



#---------------------------------------------------------------------
# Figure S3B-D; differences between sexes


#Some processing to get the names in the right order/format
mp_data$species_common <- as.character(mp_data$species)
mp_data[mp_data$species_common == "Tapirus_terrestris","species_common"] <- "Lowland tapir"
mp_data[mp_data$species_common == "Myrmecophaga_tridactyla","species_common"] <- "Giant anteater"
mp_data[mp_data$species_common == "Priodontes_maximus","species_common"] <- "Giant armadillo"
mp_data$species_common <- factor(mp_data$species_common, levels = c("Lowland tapir", "Giant armadillo", "Giant anteater"), ordered = TRUE)

#Generate the figure
mid <-
  ggplot(data = mp_data, aes(x = sex,
                             y = mp_ml+1,
                             col = sex,
                             fill = sex,
                             alpha = 0.5)) +
  ggtitle("B") +
  geom_boxplot(size = 0.1, outlier.size = 0.2, outlier.shape = 16, outlier.alpha = 0) +
  geom_jitter(size = 0.5, shape = 16, position=position_jitter(height=0, width=0.2)) +
  scale_fill_manual(values = c("#fca311", "#14213d")) +
  scale_colour_manual(values = c("#fca311", "#14213d")) +
  facet_wrap(vars(species_common),  ncol = 3) +
  ylab("Microplastics (particles/mL)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=4, family = "sans", face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans", face = "bold", color = "black"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold"),
        #strip.text.x = element_text(size=6, family = "sans", face = "bold", color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=5, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_x_discrete(breaks = c("female", "male"), labels = c("Female","Male")) +
  scale_y_log10(breaks = c(1,10+1,100+1,1000+1),
                labels = c(0,10,100,1000),
                limits = c(1,1000)) +
  annotation_logticks(#outside = TRUE,
    sides = "l",
    linewidth = 0.1,
    short = unit(0.05, "cm"),
    mid = unit(0.05, "cm"),
    long = unit(0.08, "cm"))





#---------------------------------------------------------------------
# Figure S3E-G; correlations with body weight


#Generate the figure
bot <-
  ggplot(data = mp_data, aes(x = weight,
                             y = mp_ml+1,
                             col = sex,
                             fill = sex,
                             alpha = 0.5)) +
  ggtitle("E") +
  geom_smooth(method = "gam",
              formula = y ~ x,
              method.args = list(family = tw(link = "log")),
              linewidth = 0.2,
              linetype = "solid",
              show.legend = FALSE) +
  geom_point(size = 0.5, shape = 16) +
  scale_fill_manual(values = c("#fca311", "#14213d")) +
  scale_colour_manual(values = c("#fca311", "#14213d")) +
  facet_wrap(vars(species_common),  ncol = 3, scales = "free_x") +
  ylab("Microplastics (particles/mL)") +
  xlab("Body weight (kg)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=4, family = "sans", face = "bold"),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold"),
        #strip.text.x = element_text(size=6, family = "sans", face = "bold", color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=5, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_y_log10(breaks = c(1,10+1,100+1,1000+1,10000+1),
                labels = c(0,10,100,1000,10000)) +
  annotation_logticks(#outside = TRUE,
    sides = "l",
    linewidth = 0.1,
    short = unit(0.05, "cm"),
    mid = unit(0.05, "cm"),
    long = unit(0.08, "cm"))




# Assemble and save
figure_s3 <-
  grid.arrange(a,mid,bot,
               ncol=1,
               nrow=3,
               heights=c(2,1,1))


#Save the figures
ggsave(figure_s3,
       width = 4.75, height = 5, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/figure_s3.png")



#-------------------------------------------------------------
# Figure S4 Size distributions of the particles in the blank controls
#-------------------------------------------------------------


a <- 
  ggplot() +
  ggtitle("A") +
  geom_histogram(data = control, aes(Length, fill = Sample),
                 alpha = 0.8,
                 bins = 60,
                 col = "black",
                 linewidth = 0.05) +
  scale_fill_manual(values = c("#1c7293", "#1b3b6f")) +
  scale_x_log10(expand = c(0,0.02)) +
  scale_y_continuous(limits = c(0,20), expand = c(0,.1)) +
  ylab("Number of particles") +
  xlab(bquote(bold('Particle length '(µm)))) +
  theme_bw() +
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


b <- 
  ggplot() +
  ggtitle("B") +
  geom_histogram(data = control, aes(Width, fill = Sample),
                 alpha = 0.8,
                 bins = 60,
                 col = "black",
                 linewidth = 0.05) +
  scale_fill_manual(values = c("#1c7293", "#1b3b6f")) +
  scale_x_log10(expand = c(0,0.02)) +
  scale_y_continuous(limits = c(0,20), expand = c(0,.1)) +
  ylab("Number of particles") +
  xlab(bquote(bold('Particle width '(µm)))) +
  theme_bw() +
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




FIG <-
  grid.arrange(a,b,
               ncol=2,
               nrow=1)


#Save the figures
ggsave(FIG,
       width = 4.75, height = 1.8, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/figure_s4.png")





#---------------------------------------------------------------------
# Figure S5 Rarefaction curves on polymer abundances
#---------------------------------------------------------------------


#---------------------------------------------------------------------
# Rarefaction curves for the Cerrado

nReps <- 100


# First subset to just the Cerrado data
Cerrado <- polymers[(polymers$biome == "Cerrado"),4:16]

#empty list for storing results
RES <- list()

for(i in 1:nrow(Cerrado)){
  
  #empty list for storing results
  res <- list()
  
  #loop over desired number of replicates
  for(j in 1:nReps){
    
    #Sample i samples
    data_sub <- Cerrado[sample(1:nrow(Cerrado), i),]
    
    #Get the metrics of interest
    res[[j]] <- data.frame(n = i,
                           n_polymers = vegan::specnumber(colSums(data_sub)),
                           Shannon_diversity = vegan::diversity(colSums(data_sub), index="shannon"),
                           Simpson_diversity = vegan::diversity(colSums(data_sub), index="simpson"))
    
  }# closes the loop over the replicates
  
  
  #Convert the list to a dataframe
  res <- do.call(rbind, res)
  
  #Compile the relevant summary statistics
  RES[[i]] <- data.frame(n = i,
                         n_polymers_mean = mean(res$n_polymers),
                         n_polymers_sd = sd(res$n_polymers),
                         Shannon_mean = mean(res$Shannon_diversity),
                         Shannon_sd = sd(res$Shannon_diversity),
                         Simpson_mean = mean(res$Simpson_diversity),
                         Simpson_sd = sd(res$Simpson_diversity))
  
}

#Convert to dataframe
RES <- do.call(rbind, RES)

RES$n_polymers_min <- RES$n_polymers_mean - 1.96*RES$n_polymers_sd
RES$n_polymers_max <- RES$n_polymers_mean + 1.96*RES$n_polymers_sd


a <- 
  ggplot() +
  geom_hline(yintercept = vegan::specnumber(colSums(Cerrado)), linewidth = 0.6, linetype = "dashed") +
  geom_point(data = RES, aes(x = n, y = n_polymers_mean), col = "#e9c46a", size = 0.3) +
  geom_line(data = RES, aes(x = n, y = n_polymers_mean), col = "#e9c46a", size = 0.3) +
  
  geom_ribbon(data = RES, aes(x = n, ymin = n_polymers_min, ymax = n_polymers_max), alpha = 0.2, fill = "#e9c46a") +
  ggtitle("A") +
  xlab("Number of samples") +
  ylab ("Number of polymers") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.title.y = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans", colour = "black"),
        axis.text.x  = element_text(size=4, family = "sans", colour = "black"),
        strip.text.x = element_text(size=9, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=4, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.position = "right",
        legend.position.inside = c(0.1,0.85),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm"))




RES$Shannon_min <- RES$Shannon_mean - 1.96*RES$Shannon_sd
RES$Shannon_max <- RES$Shannon_mean + 1.96*RES$Shannon_sd

b <- 
  ggplot() +
  geom_hline(yintercept = vegan::diversity(colSums(Cerrado), index="shannon"), linewidth = 0.6, linetype = "dashed") +
  geom_point(data = RES, aes(x = n, y = Shannon_mean), col = "#e9c46a", size = 0.3) +
  geom_line(data = RES, aes(x = n, y = Shannon_mean), col = "#e9c46a", size = 0.3) +
  
  geom_ribbon(data = RES, aes(x = n, ymin = Shannon_min, ymax = Shannon_max), alpha = 0.2, fill = "#e9c46a") +
  ggtitle("B") +
  xlab("Number of samples") +
  ylab ("Polymer Shannon Diversity") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.title.y = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans", colour = "black"),
        axis.text.x  = element_text(size=4, family = "sans", colour = "black"),
        strip.text.x = element_text(size=9, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=4, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.position = "right",
        legend.position.inside = c(0.1,0.85),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm"))



RES$Simpson_min <- RES$Simpson_mean - 1.96*RES$Simpson_sd
RES$Simpson_max <- RES$Simpson_mean + 1.96*RES$Simpson_sd

c <- 
  ggplot() +
  geom_hline(yintercept = vegan::diversity(colSums(Cerrado), index="simpson"), linewidth = 0.6, linetype = "dashed") +
  geom_point(data = RES, aes(x = n, y = Simpson_mean), col = "#e9c46a", size = 0.3) +
  geom_line(data = RES, aes(x = n, y = Simpson_mean), col = "#e9c46a", size = 0.3) +
  
  geom_ribbon(data = RES, aes(x = n, ymin = Simpson_min, ymax = Simpson_max), alpha = 0.2, fill = "#e9c46a") +
  ggtitle("C") +
  xlab("Number of samples") +
  ylab ("Polymer Simpson Diversity") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.title.y = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans", colour = "black"),
        axis.text.x  = element_text(size=4, family = "sans", colour = "black"),
        strip.text.x = element_text(size=9, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=4, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.position = "right",
        legend.position.inside = c(0.1,0.85),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm"))




TOP <-
  grid.arrange(a,b,c,
               ncol=3,
               nrow=1)


#---------------------------------------------------------------------
# Rarefaction curves for the Pantanal

# First subset to just the Pantanal data
Pantanal <- polymers[(polymers$biome == "Pantanal"),4:16]

#empty list for storing results
RES <- list()

for(i in 1:nrow(Pantanal)){
  
  #empty list for storing results
  res <- list()
  
  #loop over desired number of replicates
  for(j in 1:nReps){
    
    #Sample i samples
    data_sub <- Pantanal[sample(1:nrow(Pantanal), i),]
    
    #Get the metrics of interest
    res[[j]] <- data.frame(n = i,
                           n_polymers = vegan::specnumber(colSums(data_sub)),
                           Shannon_diversity = vegan::diversity(colSums(data_sub), index="shannon"),
                           Simpson_diversity = vegan::diversity(colSums(data_sub), index="simpson"))
    
  }# closes the loop over the replicates
  
  
  #Convert the list to a dataframe
  res <- do.call(rbind, res)
  
  #Compile the relevant summary statistics
  RES[[i]] <- data.frame(n = i,
                         n_polymers_mean = mean(res$n_polymers),
                         n_polymers_sd = sd(res$n_polymers),
                         Shannon_mean = mean(res$Shannon_diversity),
                         Shannon_sd = sd(res$Shannon_diversity),
                         Simpson_mean = mean(res$Simpson_diversity),
                         Simpson_sd = sd(res$Simpson_diversity))
  
}

#Convert to dataframe
RES <- do.call(rbind, RES)

RES$n_polymers_min <- RES$n_polymers_mean - 1.96*RES$n_polymers_sd
RES$n_polymers_max <- RES$n_polymers_mean + 1.96*RES$n_polymers_sd


d <- 
  ggplot() +
  geom_hline(yintercept = vegan::specnumber(colSums(Pantanal)), linewidth = 0.6, linetype = "dashed") +
  geom_point(data = RES, aes(x = n, y = n_polymers_mean), col = "#0a9396", size = 0.3) +
  geom_line(data = RES, aes(x = n, y = n_polymers_mean), col = "#0a9396", size = 0.3) +
  
  geom_ribbon(data = RES, aes(x = n, ymin = n_polymers_min, ymax = n_polymers_max), alpha = 0.2, fill = "#0a9396") +
  ggtitle("D") +
  xlab("Number of samples") +
  ylab ("Number of polymers") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.title.y = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans", colour = "black"),
        axis.text.x  = element_text(size=4, family = "sans", colour = "black"),
        strip.text.x = element_text(size=9, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=4, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.position = "right",
        legend.position.inside = c(0.1,0.85),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm"))




RES$Shannon_min <- RES$Shannon_mean - 1.96*RES$Shannon_sd
RES$Shannon_max <- RES$Shannon_mean + 1.96*RES$Shannon_sd

e <- 
  ggplot() +
  geom_hline(yintercept = vegan::diversity(colSums(Pantanal), index="shannon"), linewidth = 0.6, linetype = "dashed") +
  geom_point(data = RES, aes(x = n, y = Shannon_mean), col = "#0a9396", size = 0.3) +
  geom_line(data = RES, aes(x = n, y = Shannon_mean), col = "#0a9396", size = 0.3) +
  
  geom_ribbon(data = RES, aes(x = n, ymin = Shannon_min, ymax = Shannon_max), alpha = 0.2, fill = "#0a9396") +
  ggtitle("E") +
  xlab("Number of samples") +
  ylab ("Polymer Shannon Diversity") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.title.y = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans", colour = "black"),
        axis.text.x  = element_text(size=4, family = "sans", colour = "black"),
        strip.text.x = element_text(size=9, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=4, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.position = "right",
        legend.position.inside = c(0.1,0.85),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm"))



RES$Simpson_min <- RES$Simpson_mean - 1.96*RES$Simpson_sd
RES$Simpson_max <- RES$Simpson_mean + 1.96*RES$Simpson_sd

f <- 
  ggplot() +
  geom_hline(yintercept = vegan::diversity(colSums(Pantanal), index="simpson"), linewidth = 0.6, linetype = "dashed") +
  geom_point(data = RES, aes(x = n, y = Simpson_mean), col = "#0a9396", size = 0.3) +
  geom_line(data = RES, aes(x = n, y = Simpson_mean), col = "#0a9396", size = 0.3) +
  
  geom_ribbon(data = RES, aes(x = n, ymin = Simpson_min, ymax = Simpson_max), alpha = 0.2, fill = "#0a9396") +
  ggtitle("F") +
  xlab("Number of samples") +
  ylab ("Polymer Simpson Diversity") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.title.y = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans", colour = "black"),
        axis.text.x  = element_text(size=4, family = "sans", colour = "black"),
        strip.text.x = element_text(size=9, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=4, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.position = "right",
        legend.position.inside = c(0.1,0.85),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm"))




MID <-
  grid.arrange(d,e,f,
               ncol=3,
               nrow=1)




#---------------------------------------------------------------------
# Rarefaction curves for the Amazon

# First subset to just the Amazon data
Amazon <- polymers[(polymers$biome == "Amazon"),4:16]

#empty list for storing results
RES <- list()

for(i in 1:nrow(Amazon)){
  
  #empty list for storing results
  res <- list()
  
  #loop over desired number of replicates
  for(j in 1:nReps){
    
    #Sample i samples
    data_sub <- Amazon[sample(1:nrow(Amazon), i),]
    
    #Get the metrics of interest
    res[[j]] <- data.frame(n = i,
                           n_polymers = vegan::specnumber(colSums(data_sub)),
                           Shannon_diversity = vegan::diversity(colSums(data_sub), index="shannon"),
                           Simpson_diversity = vegan::diversity(colSums(data_sub), index="simpson"))
    
  }# closes the loop over the replicates
  
  
  #Convert the list to a dataframe
  res <- do.call(rbind, res)
  
  #Compile the relevant summary statistics
  RES[[i]] <- data.frame(n = i,
                         n_polymers_mean = mean(res$n_polymers),
                         n_polymers_sd = sd(res$n_polymers),
                         Shannon_mean = mean(res$Shannon_diversity),
                         Shannon_sd = sd(res$Shannon_diversity),
                         Simpson_mean = mean(res$Simpson_diversity),
                         Simpson_sd = sd(res$Simpson_diversity))
  
}

#Convert to dataframe
RES <- do.call(rbind, RES)

RES$n_polymers_min <- RES$n_polymers_mean - 1.96*RES$n_polymers_sd
RES$n_polymers_max <- RES$n_polymers_mean + 1.96*RES$n_polymers_sd


g <- 
  ggplot() +
  geom_hline(yintercept = vegan::specnumber(colSums(Amazon)), linewidth = 0.6, linetype = "dashed") +
  geom_point(data = RES, aes(x = n, y = n_polymers_mean), col = "#007200", size = 0.3) +
  geom_line(data = RES, aes(x = n, y = n_polymers_mean), col = "#007200", size = 0.3) +
  
  geom_ribbon(data = RES, aes(x = n, ymin = n_polymers_min, ymax = n_polymers_max), alpha = 0.2, fill = "#007200") +
  ggtitle("G") +
  xlab("Number of samples") +
  ylab ("Number of polymers") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.title.y = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans", colour = "black"),
        axis.text.x  = element_text(size=4, family = "sans", colour = "black"),
        strip.text.x = element_text(size=9, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=4, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.position = "right",
        legend.position.inside = c(0.1,0.85),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm"))




RES$Shannon_min <- RES$Shannon_mean - 1.96*RES$Shannon_sd
RES$Shannon_max <- RES$Shannon_mean + 1.96*RES$Shannon_sd

h <- 
  ggplot() +
  geom_hline(yintercept = vegan::diversity(colSums(Amazon), index="shannon"), linewidth = 0.6, linetype = "dashed") +
  geom_point(data = RES, aes(x = n, y = Shannon_mean), col = "#007200", size = 0.3) +
  geom_line(data = RES, aes(x = n, y = Shannon_mean), col = "#007200", size = 0.3) +
  
  geom_ribbon(data = RES, aes(x = n, ymin = Shannon_min, ymax = Shannon_max), alpha = 0.2, fill = "#007200") +
  ggtitle("H") +
  xlab("Number of samples") +
  ylab ("Polymer Shannon Diversity") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.title.y = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans", colour = "black"),
        axis.text.x  = element_text(size=4, family = "sans", colour = "black"),
        strip.text.x = element_text(size=9, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=4, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.position = "right",
        legend.position.inside = c(0.1,0.85),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm"))



RES$Simpson_min <- RES$Simpson_mean - 1.96*RES$Simpson_sd
RES$Simpson_max <- RES$Simpson_mean + 1.96*RES$Simpson_sd

i <- 
  ggplot() +
  geom_hline(yintercept = vegan::diversity(colSums(Amazon), index="simpson"), linewidth = 0.6, linetype = "dashed") +
  geom_point(data = RES, aes(x = n, y = Simpson_mean), col = "#007200", size = 0.3) +
  geom_line(data = RES, aes(x = n, y = Simpson_mean), col = "#007200", size = 0.3) +
  
  geom_ribbon(data = RES, aes(x = n, ymin = Simpson_min, ymax = Simpson_max), alpha = 0.2, fill = "#007200") +
  ggtitle("I") +
  xlab("Number of samples") +
  ylab ("Polymer Simpson Diversity") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.title.y = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans", colour = "black"),
        axis.text.x  = element_text(size=4, family = "sans", colour = "black"),
        strip.text.x = element_text(size=9, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=4, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.position = "right",
        legend.position.inside = c(0.1,0.85),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm"))




BOT <-
  grid.arrange(g,h,i,
               ncol=3,
               nrow=1)



FIG <-
  grid.arrange(TOP,MID,BOT,
               ncol=1,
               nrow=3)

#Save the figures
ggsave(FIG,
       width = 4.75, height = 4.5, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/figure_s5.png")
