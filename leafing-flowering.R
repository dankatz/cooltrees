library(tidyverse)
library(MASS)  # For robust linear modeling using rlm

set.seed(135)

# Number of data points
n <- 150

# Simulate temperature with a relationship to phenological phases
temperature <- rnorm(n, mean = 15, sd = 3)

# Simulate peak leaf out date (day of year) with influence from temperature
leaf_out_date <- rnorm(n, mean = 130, sd = 12) - 2 * temperature

# Create peak flowering date with a linear relationship plus noise
flowering_date <- 0.8 * leaf_out_date + rnorm(n, mean = 22, sd = 7)

# Introduce more extreme leaf outliers
leaf_outliers_indices <- sample(seq_len(n), size = 10, replace = FALSE)
leaf_out_date[leaf_outliers_indices] <- leaf_out_date[leaf_outliers_indices] - rnorm(length(leaf_outliers_indices), mean = 40, sd = 10)

# Introduce more extreme flower outliers
flower_outliers_indices <- sample(seq_len(n), size = 10, replace = FALSE)
flowering_date[flower_outliers_indices] <- flowering_date[flower_outliers_indices] - rnorm(length(flower_outliers_indices), mean = 40, sd = 10)

# Calculate difference between flowering and leafing
difference <- leaf_out_date - flowering_date

# Create a data frame with the simulated data
phenology_data <- tibble(
  temperature = temperature,
  leaf_out_date = leaf_out_date,
  flowering_date = flowering_date,
  difference = difference,
  outlier_type = case_when(
    seq_len(n) %in% leaf_outliers_indices ~ "leaf outlier",
    seq_len(n) %in% flower_outliers_indices ~ "flower outlier",
    TRUE ~ "regular point"
  )
)

# Initial regression to determine weights
initial_fit <- lm(flowering_date ~ leaf_out_date, data = phenology_data)
initial_residuals <- residuals(initial_fit)

# Assign Huber weights based on initial residuals
k <- 1.345 # Common tuning parameter in Huber weighting
weights <- pmin(1, k / abs(initial_residuals))

phenology_data <- phenology_data %>%
  mutate(weights = weights)

# Weighted Regression
weighted_fit <- lm(flowering_date ~ leaf_out_date, data = phenology_data, weights = weights)

# Robust Regression with increased maximum iterations and method set to "MM"
robust_fit <- MASS::rlm(flowering_date ~ leaf_out_date, data = phenology_data, maxit = 100, method = "MM")

# Plot flowering_date vs. leaf_out_date with different fits
ggplot(phenology_data, aes(x = leaf_out_date, y = flowering_date, color = outlier_type)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, color = "blue", se = FALSE, linetype = "solid", linewidth = 1) +
  geom_smooth(method = "lm", formula = y ~ x, aes(weight = weights), color = "green", se = FALSE, linetype = "dashed", linewidth = 1) +
  geom_smooth(method = MASS::rlm, formula = y ~ x, color = "red", se = FALSE, linetype = "dotted", linewidth = 1, n = 150, method.args = list(maxit = 100, method = "MM")) +
  scale_color_manual(values = c("flower outlier" = "cyan", "leaf outlier" = "green", "regular point" = "yellow")) +
  labs(x = "Peak leaf out (day of year)", y = "Peak flowering (day of year)", color = "Data quality") +
  theme_minimal() +
  ggtitle("Peak flowering vs. leaf out with regression fits")

# Plot difference vs. temperature with different fits
ggplot(phenology_data, aes(x = temperature, y = difference, color = outlier_type)) +
  geom_point(size = 3, aes(size = weights)) +
  geom_smooth(method = "lm", formula = y ~ x, color = "blue", se = FALSE, linetype = "solid", linewidth = 1) +
  geom_smooth(method = "lm", formula = y ~ x, aes(weight = weights), color = "green", se = FALSE, linetype = "dashed", linewidth = 1) +
  geom_smooth(method = MASS::rlm, formula = y ~ x, color = "red", linetype = "dotted", se = FALSE, linewidth = 1, n = 150, method.args = list(maxit = 100, method = "MM")) +
  labs(x = "Temperature", y = "Difference between leafing and flowering", color = "Data quality") +
  theme_minimal() +
  scale_size_continuous(name = "Weight") +
  scale_color_manual(values = c("flower outlier" = "cyan", "leaf outlier" = "green", "regular point" = "yellow")) +
  guides(size = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position = "top") +
  ggtitle("Difference between leafing and flowering with regression fits")
