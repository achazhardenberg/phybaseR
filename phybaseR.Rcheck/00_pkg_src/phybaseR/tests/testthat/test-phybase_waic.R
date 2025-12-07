test_that("phybase_waic calculates WAIC correctly", {
    skip_on_cran()

    # Setup: Run a minimal model to get a valid phybase object with an rjags model
    set.seed(123)
    tree <- ape::rtree(10)
    X <- rnorm(10)
    Y <- 0.5 * X + rnorm(10)
    data <- list(X = X, Y = Y, N = 10)
    equations <- list(Y ~ X)

    fit <- phybase_run(
        data = data,
        tree = tree,
        equations = equations,
        n.iter = 100,
        n.burnin = 50,
        n.chains = 2,
        quiet = TRUE
    )

    # Test execution
    waic <- phybase_waic(fit, n.iter = 50)

    # Check structure
    expect_type(waic, "double")
    expect_named(waic, c("waic", "p_waic"))
    expect_true(!is.na(waic["waic"]))
    expect_true(!is.na(waic["p_waic"]))
})

test_that("phybase_waic validates input", {
    expect_error(
        phybase_waic("not a model"),
        "Input must be a 'phybase' model object"
    )
})
