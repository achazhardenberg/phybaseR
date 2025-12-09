library(testthat)
library(because)

test_that("Full workflow: simple model parameter recovery stress test", {
    # Complete workflow
    set.seed(42)
    N <- 50
    tree <- ape::rtree(N)

    # Simulate correlated traits with strong signal
    # X -> Y
    X <- ape::rTraitCont(tree, model = "BM", sigma = 1)
    # Y = 0.8*X + error
    Y <- 0.8 * X + ape::rTraitCont(tree, model = "BM", sigma = 0.5)

    data_list <- list(X = X, Y = Y, N = N)
    equations <- list(Y ~ X)

    # Run model with MORE iterations
    fit <- because(
        data = data_list,
        tree = tree,
        equations = equations,
        n.iter = 5000,
        n.burnin = 1000,
        n.chains = 2,
        quiet = FALSE,
        optimise = TRUE
    )

    # Check results
    # Check parameter recovery (should be close to 0.8)
    beta_mean <- fit$summary$statistics["beta_Y_X", "Mean"]
    print(paste("Beta mean:", beta_mean))

    # Check that 95% CI includes true value (0.8)
    lower <- fit$summary$quantiles["beta_Y_X", "2.5%"]
    upper <- fit$summary$quantiles["beta_Y_X", "97.5%"]
    print(paste("CI:", lower, "-", upper))

    expect_true(lower < 0.8 && upper > 0.8)
    print("Test PASSED")
})
