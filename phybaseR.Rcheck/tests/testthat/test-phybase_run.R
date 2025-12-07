test_that("phybase_run validates inputs", {
    tree <- ape::rtree(10)

    # Missing data
    expect_error(
        phybase_run(data = NULL, tree = tree, equations = list(Y ~ X)),
        "data"
    )
})

test_that("phybase_run completes a simple run", {
    skip_on_cran() # Skip on CRAN to save time

    # Simulate simple data
    set.seed(123)
    N <- 20
    tree <- ape::rtree(N)
    X <- rnorm(N)
    Y <- 0.5 * X + rnorm(N)

    data_list <- list(X = X, Y = Y, N = N)
    equations <- list(Y ~ X)

    # Run with very few iterations just to check it runs
    fit <- phybase_run(
        data = data_list,
        tree = tree,
        equations = equations,
        n.iter = 100,
        n.burnin = 50,
        n.chains = 2,
        quiet = TRUE
    )

    expect_s3_class(fit, "phybase")
    expect_true(!is.null(fit$model))
    expect_true(!is.null(fit$samples))
})
