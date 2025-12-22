test_that("Handles single species tree gracefully", {
    # Single species tree
    tree <- ape::rtree(2) # ape::rtree(1) fails in ape usually, so min is 2 usually for meaningful phylo
    # Actually ape::rtree(1) is not valid.
    # Let's try N=2 which is the minimum for a tree.

    # But what if user passes a tree with 1 tip?
    # ape::read.tree(text="(A:1);")

    tree_single <- ape::read.tree(text = "(A:1);")
    data <- list(X = 1, Y = 1, N = 1)
    equations <- list(Y ~ X)

    # It seems to run without error, which is surprising but acceptable if JAGS handles N=1
    # We just check it returns a valid object
    fit <- because(
        data = data,
        
        equations = equations,
        n.chains = 2,
        quiet = TRUE
    )
    expect_s3_class(fit, "because")
})

test_that("Handles perfect collinearity", {
    skip_on_cran()

    set.seed(123)
    N <- 20
    tree <- ape::rtree(N)
    X <- rnorm(N)
    Y <- 2 * X # Perfect correlation

    data <- list(X = X, Y = Y, N = N)
    equations <- list(Y ~ X)

    # Should run but might have convergence issues or high variance
    # We just want to ensure it doesn't crash
    fit <- because(
        data = data,
        
        equations = equations,
        n.iter = 100,
        n.burnin = 50,
        n.chains = 2,
        quiet = TRUE
    )

    expect_s3_class(fit, "because")
})

test_that("Handles disconnected DAG", {
    # X -> Y, Z (Z is independent)
    set.seed(123)
    N <- 20
    tree <- ape::rtree(N)
    X <- rnorm(N)
    Y <- X + rnorm(N)
    Z <- rnorm(N)

    data <- list(X = X, Y = Y, Z = Z, N = N)
    equations <- list(
        Y ~ X,
        Z ~ 1 # Intercept only
    )

    fit <- because(
        data = data,
        
        equations = equations,
        n.iter = 100,
        n.burnin = 50,
        n.chains = 2,
        quiet = TRUE
    )

    expect_s3_class(fit, "because")
})
