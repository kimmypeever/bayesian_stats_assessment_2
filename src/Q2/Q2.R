# -----------------------------
# Q2 Hierarchical modeling
# -----------------------------

# Load packages and set seed
library(coda)
library(ggplot2)
library(tidyr)
library(rjags)
set.seed(8)

# Load the data into R
data <- read.table(text = "
Group1 Group2 Group3 Group4
626.952 588.610 582.906 640.996
589.944 598.012 604.918 618.210
591.879 605.204 595.259 568.970
590.125 564.930 583.211 580.336
657.360 600.661 604.727 644.002
621.158 652.042 643.959 613.107
570.384 615.596 606.690 609.925
661.516 590.955 633.964 610.026
", header = TRUE)

# Summary statistics
summary_stats <- data.frame(
  Mean = sapply(data, mean),
  Variance = sapply(data, var)
)
summary_stats


# Plotting
# Convert to long format for easy group plotting
data_long <- pivot_longer(data, cols = everything(),
                          names_to = "Group", values_to = "Value")

# Boxplot
ggplot(data_long, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplot of Measurements by Group", y = "Value")

# Density plot
ggplot(data_long, aes(x = Value, color = Group, fill = Group)) +
  geom_density(alpha = 0.3) +
  theme_minimal() +
  labs(title = "Density Plot of Measurements by Group", x = "Value") +
  xlim(530, 700) 


# --------------------------------------------------------------
# Convert group labels to numeric identifiers for indexing in the JAGS model
data_long$GroupID <- as.numeric(factor(data_long$Group))

# Check data
head(data_long)

# Extract and define necessary data
y   <- data_long$Value           # Observed measurements
idx <- data_long$GroupID         # Numeric group identifiers
N <- length(y)                   # Total number of observations
J <- length(unique(idx))         # Number of distinct groups

# Bundle data for JAGS
the_data <- list(y = y, idx = idx, N = N, J = J)

# Specify hierarchical model
model_string <- "
model {
  # Likelihood
  for(i in 1:N) {
    y[i] ~ dnorm(mu[idx[i]], tau)
  }

  # Group-level means drawn from a common prior
  for(j in 1:J) {
    mu[j] ~ dnorm(tau_mu, tau_psi)
  }

  # Hyperpriors for group-level mean and precision
  tau_mu ~ dnorm(600, 1/100^2)     
  tau_psi ~ dgamma(0.5, 250)       
  tau ~ dgamma(0.5, 250)          
}
"

# Initialize JAGS model
jags_model <- jags.model(textConnection(model_string),
                         data = the_data,
                         n.chains = 1,
                         n.adapt = 500)

# Posterior samples
samples <- coda.samples(jags_model,
                        variable.names = c("mu", "tau_mu", "tau_psi", "tau"),
                        n.iter = 1000)

# Diagnostics
plot(samples)
traceplot(samples)            
autocorr.plot(samples)   

summary(samples)

#----------------------------------------------------------------------
# 4 chains with dispersed starting values
set.seed(8)
init_values <- list(
  list(mu = rnorm(J, 590, 5), tau_mu = 590, tau_psi = 1, tau = 1),
  list(mu = rnorm(J, 610, 5), tau_mu = 610, tau_psi = 2, tau = 2),
  list(mu = rnorm(J, 600, 5), tau_mu = 600, tau_psi = 0.5, tau = 0.5),
  list(mu = rnorm(J, 605, 5), tau_mu = 605, tau_psi = 1.5, tau = 1.5)
)

jags_model <- jags.model(textConnection(model_string),
                         data = the_data,
                         inits = init_values,
                         n.chains = 4,
                         n.adapt = 50)  # 50 adaptation iterations

samples_multi <- coda.samples(jags_model,
                              variable.names = c("mu", "tau_mu", "tau_psi", "tau"),
                              n.iter = 50)  # 50 iterations per chain

# Superimposed trace plots
plot(samples_multi)

# Plot tracew
par(mfrow=c(3,3))
traceplot(samples_multi)
par(mfrow=c(1,1))
traceplot(samples_multi)


# Autocorrelation plots
autocorr.plot(samples_multi)

# Gelman-Rubin diagnostic
gelman.diag(samples_multi)

#-------------------------------------------------------------------
# Final 10,000 chain
jags_model_final <- jags.model(
  textConnection(model_string),
  data = the_data,
  n.chains = 1,
  n.adapt = 500
)

samples_final <- coda.samples(
  jags_model_final,
  variable.names = c("mu", "tau_mu", "tau_psi", "tau"),
  n.iter = 10000,
  thin = 5
)

# Discard burn-in
burn_in <- 500
samples_final_burned <- window(samples_final, start = burn_in + 1)

# Summary stats
summary(samples_final_burned)

# Compute 95% equal-tailed credible intervals for mu[1:4]
mu_samples <- as.matrix(samples_final_burned)[, grep("mu\\[", colnames(as.matrix(samples_final_burned)))]

# 2.5% and 97.5% quantiles for each mu
credible_intervals <- apply(mu_samples, 2, quantile, probs = c(0.025, 0.975))
credible_intervals



