library(ggplot2)
library(gridExtra)


#Import the MP datasets
source("scripts/data_import.R")

#-------------------------------------------------------------
# Particle sizes for the different polymers
#-------------------------------------------------------------


aggregate(Length~ polymer, FUN = "length", data = sizes)

#a <- 
  ggplot() +
  ggtitle("A") +
  geom_histogram(data = sizes, aes(Length, fill = polymer),
                 alpha = 0.8,
                 bins = 60,
                 col = "black",
                 linewidth = 0.05) +
  #scale_fill_manual(values = c("#619b8a", "#bb3e03", "#005f73")) +
  scale_x_log10(expand = c(0,0.1)) +
  #scale_y_continuous(limits = c(0,230), expand = c(0,.1)) +
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
  guides(fill = guide_legend(byrow = TRUE)) #+
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





#-------------------------------------------------------------
# Particle sizes in the blank controls
#-------------------------------------------------------------

#Sumary statstics of the particle sizes
mean(control$Length)

sd(control$Length)

range(control$Length)

mean(control$Width)

sd(control$Width)

range(control$Width)

nrow(control)

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
  scale_y_continuous(limits = c(0,200), expand = c(0,.1)) +
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




bot <-
  grid.arrange(a,b,
               ncol=2,
               nrow=1)


#Save the figures
ggsave(bot,
       width = 4.75, height = 1.8, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/particle_size_control.png")


