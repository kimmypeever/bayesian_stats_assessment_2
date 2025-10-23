# -----------------------------
# Q1 Bayesian analysis of hourly birth rates
# -----------------------------

# Load in packages
library(ggplot2)

# Data
counts <- c(0, 1, 2, 3, 4)      # birth counts per hour category
freqs  <- c(3, 8, 6, 4, 3)      # number of hours with each count
y      <- rep(counts, freqs)    # expand to individual hourly observations

n <- length(y)
sum_y <- sum(y)

# Prior parameters (from moment matching)
alpha <- 4
beta  <- 2

# Simulate prior draws of lambda
lambda_prior <- rgamma(10000, shape = alpha, rate = beta)

# Simulate hourly birth counts from Poisson using prior lambdas
y_prior <- rpois(10000, lambda_prior)

# Summarise simulated counts
summary(y_prior)
table(y_prior)

# Visualise prior predictive distribution
hist(y_prior, breaks = 0:15 - 0.5, col = "lightblue", border = "white",
     main = "Prior Predictive Distribution of Hourly Birth Counts",
     xlab = "Births per Hour", ylab = "Frequency")

# Posterior parameters
alpha_post <- alpha + sum_y
beta_post  <- beta + n

# Posterior summaries
post_mean <- alpha_post / beta_post
post_var  <- alpha_post / beta_post^2
post_sd   <- sqrt(post_var)

# 95% credible interval (central)
ci_lower <- qgamma(0.025, shape = alpha_post, rate = beta_post)
ci_upper <- qgamma(0.975, shape = alpha_post, rate = beta_post)

# Frequentist estimate (MLE)
lambda_mle <- mean(y)

# Print results
posterior_summary <- data.frame(
  Statistic = c("Posterior mean", "Posterior SD", "95% CI lower", "95% CI upper", "MLE (freq.)"),
  Value = round(c(post_mean, post_sd, ci_lower, ci_upper, lambda_mle), 3)
)
print(posterior_summary)

# Plot posterior density
lambda_vals <- seq(0, 5, length.out = 1000)
posterior_density <- dgamma(lambda_vals, shape = alpha_post, rate = beta_post)
posterior_df <- data.frame(lambda_vals, posterior_density)

# Posterior dist with post_mean and MLE
ggplot(data.frame(lambda = lambda_vals, density = posterior_density),
       aes(x = lambda, y = density)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_vline(xintercept = lambda_mle, linetype = "dashed", color = "darkred") +
  geom_vline(xintercept = post_mean, linetype = "dotted", color = "black") +
  labs(title = "Posterior Distribution of Hourly Birth Rate",
       x = expression(lambda), y = "Density") +
  theme_minimal()


# Posterior dist with credible interval
ggplot(data.frame(lambda_vals, posterior_density), aes(x = lambda_vals, y = posterior_density)) +
  geom_area(
    data = subset(posterior_df, lambda_vals >= ci_lower & lambda_vals <= ci_upper),
    aes(x = lambda_vals, y = posterior_density),
    fill = "blue", alpha = 0.2
  ) +
  geom_line(color = "black", size = 1) +
  geom_vline(xintercept = post_mean, color = "red", linetype = "solid", size = 0.8) +
  geom_vline(xintercept = c(ci_lower, ci_upper), color = "red", linetype = "dashed", size = 0.6) +
  annotate("text", x = post_mean, y = 0.75, label = expression(paste("Posterior mean")),
           color = "red", angle = 90, vjust = -0.5, hjust = 0) +
  annotate("text", x = ci_lower+0.05, y = 0.1, label = "95% CI lower", color = "red", angle = 90, vjust = 1.2, size = 4) +
  annotate("text", x = ci_upper-0.1, y = 0.1, label = "95% CI upper", color = "red", angle = 90, vjust = -0.2, size = 4) +
  labs(
    title = "Posterior Distribution of Hourly Birth Rate",
    x = expression(lambda),
    y = "Density"
  ) +
  theme_minimal(base_size = 13)

