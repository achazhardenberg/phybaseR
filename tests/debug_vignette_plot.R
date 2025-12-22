library(because)

set.seed(123)
N <- 100
Temp_Centered <- rnorm(N, mean = 20, sd = 5) - 20
Growth_g_day <- 0.5 * Temp_Centered + 10 + rnorm(N, sd = 2)
df <- data.frame(Temp_Centered, Growth_g_day)

fit_default <- because(
    equations = list(Growth_g_day ~ Temp_Centered),
    data = df,
    n.iter = 100,
    quiet = TRUE
)

# Mock fit_custom identical for testing structure
fit_custom <- fit_default

get_est <- function(fit, param) {
    if (is.null(fit$summary)) {
        return(c(NA, NA, NA))
    }
    if (!param %in% rownames(fit$summary$statistics)) {
        return(rep(NA, 3))
    }

    stats <- fit$summary$statistics[param, ]
    quant <- fit$summary$quantiles[param, ]
    val <- c(
        mean = stats["Mean"],
        lower = quant["2.5%"],
        upper = quant["97.5%"]
    )
    return(val)
}

p_slope <- "beta_Growth_g_day_Temp_Centered"
p_int <- "alphaGrowth_g_day"

def_slope <- get_est(fit_default, p_slope)
print("Default Slope est:")
print(def_slope)

def_int <- get_est(fit_default, p_int)
print("Default Int est:")
print(def_int)

if (any(is.na(def_slope))) {
    message("Slope has NAs")
}
if (any(is.na(def_int))) {
    message("Int has NAs")
}

# Check range
y_lim_slope <- range(c(def_slope, def_slope), na.rm = TRUE)
message("Ylim Slope:")
print(y_lim_slope)
