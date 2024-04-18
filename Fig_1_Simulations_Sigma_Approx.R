graphics.off() 
rm(list=ls(all=TRUE)) 

# Load required libraries
library(ggplot2)

# Initialize variables
n_samples <- 1000  # Number of samples per simulation
sigmas <- seq(0.1, 6, by=0.1)  # Different sigma values to try
n_simulations <- 100  # Number of simulations per sigma value
results <- data.frame()

# Run simulation
for (sigma in sigmas) {
  for (i in 1:n_simulations) {
    # Generate lognormally-distributed random variables
    x <- rlnorm(n_samples, meanlog=0, sdlog=sigma)
    
    # Calculate sample median and standard deviation
    sample_median <- median(x)
    sample_sd <- sd(x)
    
    # Estimate sigma using the approximation formula
    estimated_sigma <- sqrt(log(sample_sd / sample_median))
    
    # Calculate the true and estimated arithmetic mean
    true_mean <- exp(0 + sigma^2 / 2)
    estimated_mean <- exp(0 + estimated_sigma^2 / 2)
    
    # Store the results
    results <- rbind(results, data.frame(TrueSigma=sigma, EstimatedSigma=estimated_sigma,
                                         TrueMean=true_mean, EstimatedMean=estimated_mean))
  }
}

tiff(file = paste("Sigma_Approx.tiff", sep = ""), width = 2500, height = 2500, 
     units = "px",res = 1000, pointsize = 7)

sigma_stats <- results %>% group_by(TrueSigma) %>% 
  summarise(Median = median(EstimatedSigma, na.rm = TRUE),
            Lower = quantile(EstimatedSigma, 0.025, na.rm = TRUE), 
            Upper = quantile(EstimatedSigma, 0.975, na.rm = TRUE))

ggplot(results, aes(x = TrueSigma)) +
  geom_point(aes(y = EstimatedSigma), alpha = 0.1) +
  geom_ribbon(data = sigma_stats, aes(ymin = Lower, ymax = Upper), fill = "green", alpha = 0.6) +
  geom_line(data = sigma_stats, aes(y = Median), color = "green", size = 1.0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  xlab(expression(paste("True ", sigma))) +
  ylab(expression(paste("Estimated ", sigma))) +
  scale_y_continuous(limits = c(0, 6)) +
  scale_x_continuous(limits = c(0, 6)) +
  theme_minimal() +
  theme(axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black"), 
        panel.grid.minor = element_blank(), axis.text = element_text(color = "black"))

dev.off()

# Calculate median and 95% interval for each TrueMean
mean_stats <- results %>% group_by(TrueMean) %>% 
  summarise(Median = median(EstimatedMean, na.rm = TRUE),
            Lower = quantile(EstimatedMean, 0.025, na.rm = TRUE), 
            Upper = quantile(EstimatedMean, 0.975, na.rm = TRUE))

# Plot true vs estimated arithmetic mean with median, 95% interval, and dashed 1:1 line
tiff(file = "Mean_Approx.tiff", width = 3000, height = 2500, units = "px",res = 1000, pointsize = 7)
ggplot(results, aes(x = TrueMean)) +
  geom_point(aes(y = EstimatedMean), alpha = 0.2) +
  geom_line(data = mean_stats, aes(y = Median), color = "green", size = 1.0) +  # Increased line width
  geom_ribbon(data = mean_stats, aes(ymin = Lower, ymax = Upper), fill = "green", alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  xlab("True Arithmetic Mean") +
  ylab("Estimated Arithmetic Mean") +
  scale_x_log10(limits = c(1, 10^8)) +
  scale_y_log10(limits = c(1, 10^8)) +
  theme_minimal() +
  theme(axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black"), 
        panel.grid.minor = element_blank(), axis.text = element_text(color = "black"))
dev.off()
