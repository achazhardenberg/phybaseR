test_that("Gaussian distribution generates correct JAGS code", {
    equations <- list(Y ~ X)

    model_output <- phybase_model(equations)

    # Should use dnorm (normal distribution)
    expect_match(model_output$model, "dnorm")

    # Should NOT use logit or dbern
    expect_false(grepl("logit", model_output$model))
    expect_false(grepl("dbern", model_output$model))

    # Should have standard beta parameters
    expect_match(model_output$model, "betaX")
})

test_that("Binomial distribution generates correct JAGS code", {
    equations <- list(BinaryOutcome ~ X)
    distribution <- c(BinaryOutcome = "binomial")

    model_output <- phybase_model(
        equations,
        distribution = distribution
    )

    # Should use dbern (Bernoulli distribution)
    expect_match(model_output$model, "dbern")

    # Should use logit link function
    expect_match(model_output$model, "logit")

    # Should have beta parameters
    expect_match(model_output$model, "betaX")
})

test_that("Mixed distributions work correctly", {
    # Test that we can have both Gaussian and binomial in same model
    equations <- list(
        Binary ~ X
    )
    distribution <- c(Binary = "binomial")

    model_output <- phybase_model(
        equations,
        distribution = distribution
    )

    # Should have dbern for binomial
    expect_match(model_output$model, "dbern")
    expect_match(model_output$model, "logit")
})

test_that("Binomial distribution runs successfully", {
    skip_on_cran()

    set.seed(123)
    N <- 30
    tree <- ape::rtree(N)
    X <- rnorm(N)

    # Simulate binary outcome
    prob <- plogis(0.5 * X)
    Y <- rbinom(N, 1, prob)

    data <- list(X = X, Y = Y, N = N)
    equations <- list(Y ~ X)
    distribution <- c(Y = "binomial")

    fit <- phybase_run(
        data = data,
        tree = tree,
        equations = equations,
        distribution = distribution,
        n.iter = 200,
        n.burnin = 100,
        n.chains = 2,
        quiet = TRUE
    )

    expect_s3_class(fit, "phybase")
    expect_true("betaX" %in% rownames(fit$summary$statistics))
})

test_that("Default distribution is Gaussian", {
    equations <- list(Y ~ X)

    # Not specifying distribution should default to Gaussian
    model_output <- phybase_model(equations)

    expect_match(model_output$model, "dnorm")
    expect_false(grepl("dbern", model_output$model))
})
