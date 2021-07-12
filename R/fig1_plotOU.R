
# ####################################
# Description: plotting trajectories of OU process as an ARH(1) process
# Authors: A. López-Pérez and J. Álvarez-Liébana
# Article: «A goodness-of-fit test for functional time series: applications
#          to specification test for stochastic diffusion models» (submitted)
# ####################################

# Packages and working space
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
repos <- "http://cran.us.r-project.org"
if(!require(sde)) install.packages("sde", repos = repos)
if(!require(ggplot2)) install.packages("ggplot2", repos = repos)
if(!require(latex2exp)) install.packages("latex2exp", repos = repos)
if(!require(tidyverse)) install.packages("tidyverse", repos = repos)
if(!require(ggthemes)) install.packages("ggthemes", repos = repos)

# Fixing the seed for the randomness
set.seed(123)

# We simulate paths of continuous-time zero-mean stochastic
# OU process as CKLS process with kappa = 2, sigma = 0.001,
# with 250 grid points, trajectories splitted into subintervals
# [0, h], with h = 50.
t <- 0:250 
x <- sde.sim(X0 = 0, theta = c(0, 2, sqrt(0.001)), N = 250,
              t0 = 0, delta = 1/50, model = "OU")

# We interpolate the trajectories into 10 001 points
interp <- spline(t, x, n = 1e4 + 1)
OU <- as.tibble(data.frame("t" = interp$x, "x" = interp$y,
                           "interv" = as.factor(pmin(floor(interp$x/50), 4))))

# Plot the trajectories
theme_set(theme_bw(base_family = "Poppins"))
fig1 <- 
  ggplot(OU, aes(x = t, y = x, color = interv)) +
  geom_line(size = 1.3) + # width of line
  # Vertical lines
  geom_vline(aes(xintercept = 50), linetype = "dashed", size = 1.2) +
  geom_vline(aes(xintercept = 100), linetype = "dashed", size = 1.2) +
  geom_vline(aes(xintercept = 150), linetype = "dashed", size = 1.2) +
  geom_vline(aes(xintercept = 200), linetype = "dashed", size = 1.2) +
  # Color pattern from tableau appareance according to intervals
  scale_color_tableau() +
  # Annotations: curve + text
  annotate(geom = "curve", x = 15, y = 0.04, xend = 40, yend = 0.015, 
           color = "firebrick", curvature = -0.5, size = 1.1,
           arrow = arrow(length = unit(3.5, "mm"))) +
  annotate(geom = "text", x = 11, y = 0.04,
           label = TeX("$X_{1}(t)$", bold = TRUE),
           hjust = "right", family = "Poppins") +
  annotate(geom = "curve", x = 80, y = 0.035, xend = 60, yend = 0.015, 
           color = "firebrick", curvature = 0.3, size = 1.1,
           arrow = arrow(length = unit(3.5, "mm"))) +
  annotate(geom = "text", x = 81, y = 0.035,
           label = TeX("$X_{2}(t)$", bold = TRUE),
           hjust = "left", family = "Poppins") +
  annotate(geom = "curve", x = 120, y = 0.032, xend = 140, yend = 0.003, 
           color = "firebrick", curvature = -0.7, size = 1.1,
           arrow = arrow(length = unit(3.5, "mm"))) +
  annotate(geom = "text", x = 119, y = 0.032,
           label = TeX("$X_{3}(t)$", bold = TRUE),
           hjust = "right", family = "Poppins") +
  annotate(geom = "curve", x = 180, y = 0.017, xend = 163, yend = 0.005, 
           color = "firebrick", curvature = 0.3, size = 1.1,
           arrow = arrow(length = unit(3.5, "mm"))) +
  annotate(geom = "text", x = 181, y = 0.017,
           label = TeX("$X_{4}(t)$", bold = TRUE),
           hjust = "left", family = "Poppins")  +
  annotate(geom = "curve", x = 220, y = 0.037, xend = 240, yend = 0.004, 
           color = "firebrick", curvature = -0.6, size = 1.1,
           arrow = arrow(length = unit(3.5, "mm"))) +
  annotate(geom = "text", x = 219, y = 0.037,
           label = TeX("$X_{5}(t)$", bold = TRUE),
           hjust = "right", family = "Poppins") +
  # Axis
  labs(x = "t",
       y = TeX("Continuous-time zero-mean process $\\xi_{nh+t}$",
               bold = TRUE)) +
  # Title (splitted in two lines)
  ggtitle("Splitting trajectories\nof an OU process") +
  # Font settings
  theme(axis.title.x = element_text(family = "Poppins", hjust = .5,
                                    size = 11, face = "bold"),
        axis.title.y = element_text(family = "Poppins", hjust = .5,
                                    size = 11, face = "bold"),
        plot.title = element_text(face = "bold", family = "Poppins",
                                  size = 13),
        legend.position = "none") # without legend

# Save as 8x4 png 
ggsave("fig1_plotOU.png", width = 8, height = 4)



