library(ape)
library(MASS)

simulate_plsm_data_correlated <- function(N = 100, seed = 123, rho = 0.7) {
    set.seed(seed)

    # 1. Phylogeny
    tree <- rcoal(N)
    C <- vcv(tree)
    C <- C / max(C)

    # 2. Predictor
    x <- rnorm(N, 0, 1)

    # 3. Correlated Random Effects
    # Trait Covariance Matrix (Between Mean and Log-Variance)
    # variances for u_m and u_v
    sig_m <- 0.8 # signal strength mean
    sig_v <- 0.6 # signal strength variance

    Sigma_traits <- matrix(
        c(
            sig_m^2,
            rho * sig_m * sig_v,
            rho * sig_m * sig_v,
            sig_v^2
        ),
        nrow = 2
    )

    # Total Covariance = Sigma_traits (x) C
    # This matches the "Separable" assumption
    Sigma_tot <- kronecker(Sigma_traits, C)

    # Simulate combined vector [u_m_1, u_v_1, u_m_2, u_v_2, ...] OR [u_m_all, u_v_all]?
    # kronecker(Sigma, C) usually implies:
    # Rows organized by: Species 1 (Trait 1, Trait 2), Species 2... if Sigma (x) C?
    # Wait:
    # kronecker(A, B): A[1,1]B ...
    # If we want [u_m_1... u_m_N, u_v_1 ... u_v_N], then we want Sigma (x) C?
    # Let's check: Cov(u_m_i, u_m_j) = sig_m^2 * C_ij. This is block 1,1.
    # Cov(u_v_i, u_v_j) = sig_v^2 * C_ij. This is block 2,2.
    # Cov(u_m_i, u_v_j) = rho*sig*sig * C_ij.
    # So Yes, Kronecker(Sigma_traits, C) gives a matrix with 4 blocks: MM, MV, VM, VV.
    # The vector structure is: c(u_m_1...N, u_v_1...N).

    u_vec <- mvrnorm(1, mu = rep(0, 2 * N), Sigma = Sigma_tot)

    u_m <- u_vec[1:N]
    u_v <- u_vec[(N + 1):(2 * N)]

    # 4. Parameters
    alpha_m <- 0.5
    beta_m <- 1.2

    alpha_v <- -1.0
    beta_v <- 0.5

    # 5. Response
    mu_true <- alpha_m + beta_m * x + u_m
    log_sigma_true <- alpha_v + beta_v * x + u_v
    sigma_true <- exp(log_sigma_true)

    y <- rnorm(N, mean = mu_true, sd = sigma_true)

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
            alpha_v = alpha_v,
            beta_v = beta_v,
            sig_m = sig_m,
            sig_v = sig_v,
            rho = rho
        ),
        latent = list(
            u_m = u_m,
            u_v = u_v
        )
    )
}
