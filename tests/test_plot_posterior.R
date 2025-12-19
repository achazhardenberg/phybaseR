library(because)

# Setup Data
set.seed(123)
N <- 100
Temp_Centered <- rnorm(N, mean = 20, sd = 5) - 20
Growth_g_day <- 0.5 * Temp_Centered + 10 + rnorm(N, sd = 2)
df <- data.frame(Temp_Centered, Growth_g_day)

# 1. Fit Two Models
message("Fitting Model 1 (Default)...")
fit1 <- because(
    equations = list(Growth_g_day ~ Temp_Centered),
    data = df,
    n.iter = 500,
    quiet = TRUE
)

message("Fitting Model 2 (Custom)...")
priors <- list(alphaGrowth_g_day = "dnorm(10, 100)")
fit2 <- because(
    equations = list(Growth_g_day ~ Temp_Centered),
    data = df,
    priors = priors,
    n.iter = 500,
    quiet = TRUE
)

# 2. Test Single Model Plot
message("Testing plot_posterior with single model...")
pdf("test_plot_single.pdf")
plot_posterior(fit1, parameter = "alpha")
dev.off()

# 3. Test Comparison
message("Testing plot_posterior with comparison info...")
pdf("test_plot_compare.pdf")
plot_posterior(list(Default = fit1, Custom = fit2), parameter = "alpha")
dev.off()

# 4. Test Regex (should plot alpha and beta)
message("Testing plot_posterior with regex...")
pdf("test_plot_regex.pdf")
plot_posterior(fit1, parameter = "Growth") # alphaGrowth..., beta_...
dev.off()

message("SUCCESS: All plots generated.")
