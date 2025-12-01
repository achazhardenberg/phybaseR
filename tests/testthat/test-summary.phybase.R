test_that("summary.phybase handles standard output", {
    # Mock a phybase object
    samples <- coda::mcmc.list(
        coda::mcmc(matrix(
            rnorm(100),
            ncol = 2,
            dimnames = list(NULL, c("alpha", "beta"))
        )),
        coda::mcmc(matrix(
            rnorm(100),
            ncol = 2,
            dimnames = list(NULL, c("alpha", "beta"))
        ))
    )

    fit <- list(
        samples = samples,
        dsep = FALSE,
        DIC = c(deviance = 100, pD = 2, DIC = 102),
        WAIC = c(waic = 100, p_waic = 2)
    )
    class(fit) <- "phybase"

    # Capture output to verify printing
    output <- capture_output({
        summ <- summary(fit)
    })

    expect_s3_class(summ, "summary.mcmc")
    expect_match(output, "DIC")
    expect_match(output, "WAIC")
})

test_that("summary.phybase handles d-sep output", {
    # Mock a d-sep phybase object
    # Use 2 parameters to ensure summary returns a matrix
    samples <- coda::mcmc.list(
        coda::mcmc(matrix(
            rnorm(400, mean = 0),
            ncol = 2,
            dimnames = list(NULL, c("beta_Y_X", "beta_Z_Y"))
        ))
    )

    # Setup d-sep specific fields
    dsep_tests <- list(Y ~ X)
    attr(dsep_tests[[1]], "test_var") <- "X"

    parameter_map <- data.frame(
        response = "Y",
        predictor = "X",
        parameter = "beta_Y_X",
        stringsAsFactors = FALSE
    )

    fit <- list(
        samples = samples,
        dsep = TRUE,
        dsep_tests = dsep_tests,
        parameter_map = parameter_map
    )
    class(fit) <- "phybase"

    output <- capture_output({
        results <- summary(fit)
    })

    expect_match(output, "PhyBaSE d-separation Tests")
    expect_match(output, "Y _\\|\\|_ X")
    expect_s3_class(results, "data.frame")
    expect_true(nrow(results) == 1)
    expect_equal(results$Indep[1], "Yes") # Should be independent as mean is 0
})
