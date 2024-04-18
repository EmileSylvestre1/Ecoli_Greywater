# Clear the environment and close all graphics
graphics.off()
rm(list = ls(all = TRUE))

# Load the necessary library
library(ggplot2)

# Create a data frame with your data
data <- data.frame(
  Category = rep(c("LW", "LR", "SB"), each = 6),
  SubCategory = rep(c("O'Toole et al. (2012) - Sample", "O'Toole et al. (2012) - Lognormal approx.",
                      "O'Toole et al. (2012) - Poisson gamma", "O'Toole et al. (2012) - Poisson lognormal",
                      "Jahne et al. (2017) - Hierarchical LN-PERT model", "Friedler et al. (2004) - Sample"), times = 3),
  Line = factor(paste(rep(c("LW", "LR", "SB"), each = 6), 
                      rep(c("O'Toole et al. (2012) - Sample", "O'Toole et al. (2012) - Lognormal approx.",
                            "O'Toole et al. (2012) - Poisson gamma", "O'Toole et al. (2012) - Poisson lognormal",
                            "Jahne et al. (2017) - Hierarchical LN-PERT model", "Friedler et al. (2004) - Sample"), times = 3), sep = " - ")),
  Mean = c(1.1E+05, 1.34E+03, 1.2E+05, 6.2E+06, 5.5E+04, 4.0E+06, 3.4E+03, 1.67E+02, 3.6E+03, 2.4E+05, 8.3E+03, 4.0E+06, 1.7E+03, 7.62E+02, 1.7E+03, 1.6E+04, 1.3E+04, 4.0E+06),
  Lower_CI = c(NA, NA, 3.5E+04, 3.2E+03, NA, NA, NA, NA, 1.1E+03, 2.9E+02, NA, NA, NA, NA, 7.8E+02, 4.6E+02, NA, NA),
  Upper_CI = c(NA, NA, 3.3E+05, 7.1E+09, NA, NA, NA, NA, 9.3E+03, 1.6E+08, NA, NA, NA, NA, 3.5E+03, 3.9E+05, NA, NA)
)

tiff(file = paste("bar_plot.tiff", sep = ""), width = 7000, height = 4200, 
     units = "px",res = 1000, pointsize = 7)

# Define custom colors
custom_colors <- c("O'Toole et al. (2012) - Sample" = "black", "O'Toole et al. (2012) - Lognormal approx." = "Purple",
                   "O'Toole et al. (2012) - Poisson gamma" = "#007BFF", "O'Toole et al. (2012) - Poisson lognormal" = "#28A745",
                   "Jahne et al. (2017) - Hierarchical LN-PERT model"= "orange","Friedler et al. (2004) - Sample"= "firebrick" )

# Create the plot
p <- ggplot(data, aes(x = Mean, y = Line, color = SubCategory)) +
  geom_point() +  # Use points to represent the mean
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2) +  # Horizontal error bars
  scale_x_log10(limits = c(1e2, NA), breaks = c(1e2, 1e4, 1e6, 1e8, 1e10)) +  # Set x-axis to log scale with specified breaks
  scale_color_manual(values = custom_colors) +  # Use custom colors
  labs(x = expression(paste("Arithmetic mean ", italic("E. coli"), " or faecal coliform", " concentration (#/100 mL)")), y = "") +
  theme_minimal() +
  theme(legend.position = "none") +  # Remove the legend
  theme(plot.title = element_blank()) +  # Remove title
  theme(panel.grid = element_blank()) +  # Remove gridlines
  theme(axis.text.x = element_text(color = "black"),  # Black x-axis text
        axis.text.y = element_text(color = "black", hjust = 0),  # Black x-axis text
        axis.ticks.x = element_line(color = "black"),  # Black x-axis ticks
        axis.title = element_text(size=9,face="bold"),
        axis.line.x = element_line(color = "black", linewidth = 0.5))  # Black line for x-axis
p

dev.off()
