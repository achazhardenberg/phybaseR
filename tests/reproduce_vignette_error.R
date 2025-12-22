library(because)

set.seed(123)
N <- 100
# Ambient temperature in Celsius (Mean 20C)
Temp_Raw <- rnorm(N, mean = 20, sd = 5)

# Center it
Temp_Centered <- Temp_Raw - 20

# True relationship
Growth_g_day <- 0.5 * Temp_Centered + 10 + rnorm(N, sd = 2)

df <- data.frame(Temp_Centered, Growth_g_day)

message("Fitting default model...")
fit_default <- because(
    equations = list(Growth_g_day ~ Temp_Centered),
    data = df,
    n.iter = 100, # Short run for speed
    n.adapt = 10,
    quiet = TRUE
)

message("Parameter definitions:")
print(rownames(fit_default$summary$statistics))

p_slope <- "beta_Growth_g_day_Temp_Centered"
p_int <- "alpha_Growth_g_day"

if (p_slope %in% rownames(fit_default$summary$statistics)) {
    message("Slope param found.")
} else {
    message(paste("Slope param NOT found:", p_slope))
}

if (p_int %in% rownames(fit_default$summary$statistics)) {
    message("Intercept param found.")
} else {
    message(paste("Intercept param NOT found:", p_int))
}
