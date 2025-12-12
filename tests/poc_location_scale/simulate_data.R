library(ape)
library(MASS)

simulate_plsm_data <- function(N = 100, seed = 123) {
    set.seed(seed)

    # 1. Simulate Phylogeny
    tree <- rcoal(N)
    # Standardize covariance matrix to max 1 (typical in comparative methods)
    C <- vcv(tree)
    C <- C / max(C)

    # 2. Simulate Independent Variables
    # Continuous predictor
    x <- rnorm(N, 0, 1)

    # 3. Simulate Phylogenetic Random Effects
    # We assume independent phylogenetic signals for mean and variance for this POC
    # Parameters for Mean Model
    alpha_m <- 0.5
    beta_m <- 1.2
    lambda_m <- 0.8 # Phylogenetic signal strength for mean (scale of RE)

    # Parameters for Variance Model (Log-Link)
    alpha_v <- -1.0 # Base log-std dev
    beta_v <- 0.5 # Effect of X on log-variance (heteroscedasticity)
    lambda_v <- 0.4 # Phylogenetic signal strength for variance

    # Simulate Random Effects
    # u_mean ~ N(0, lambda_m * C)
    u_mean <- mvrnorm(1, mu = rep(0, N), Sigma = lambda_m * C)

    # u_var ~ N(0, lambda_v * C)
    u_var <- mvrnorm(1, mu = rep(0, N), Sigma = lambda_v * C)

    # 4. Construct Latent Variables
    # True Mean
    mu_true <- alpha_m + beta_m * x + u_mean

    # True Log-Standard Deviation and Sigma
    log_sigma_true <- alpha_v + beta_v * x + u_var
    sigma_true <- exp(log_sigma_true)

    # 5. Simulate Response
    y <- rnorm(N, mean = mu_true, sd = sigma_true)

    # Return list of everything needed
    list(
        data = list(
            y = y,
            x = x,
            N = N,
            C = C
        ),
        params = list(
            alpha_m = alpha_m,
            beta_m = beta_m,
            lambda_m = lambda_m,
            alpha_v = alpha_v,
            beta_v = beta_v,
            lambda_v = lambda_v
        ),
        latent = list(
            u_mean = u_mean,
            u_var = u_var,
            sigma = sigma_true
        ),
        tree = tree
    )
}
