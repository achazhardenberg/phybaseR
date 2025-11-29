# Test script for phybaseR package
# This script tests all core functionality
cat("=== Testing phybaseR Package ===\n\n")

# Load the package functions
source("R/phybase_model.R")
source("R/phybase_run.R")
source("R/phybase_dsep.R")

# Load required libraries
library(ape)
library(rjags)

# Test 1: Load example data
cat("Test 1: Loading example data...\n")
load("data/rhino.tree.rda")
load("data/rhino.dat.rda")

if (exists("rhino.tree") && exists("rhino.dat")) {
    cat("  ✓ Data loaded successfully\n")
    cat("    - Tree has", length(rhino.tree$tip.label), "tips\n")
    cat(
        "    - Data has",
        nrow(rhino.dat),
        "rows and",
        ncol(rhino.dat),
        "columns\n"
    )
} else {
    cat("  ✗ Failed to load data\n")
    stop("Data loading failed")
}

# Test 2: Define equations (from README example)
cat("\nTest 2: Defining structural equations...\n")
equations <- list(LS ~ BM, NL ~ BM + RS, DD ~ NL)
cat("  ✓ Equations defined:\n")
for (i in seq_along(equations)) {
    cat("    ", deparse(equations[[i]]), "\n")
}

# Test 3: Generate JAGS model
cat("\nTest 3: Generating JAGS model with phybase_model()...\n")
tryCatch(
    {
        model_string <- phybase_model(equations, multi.tree = FALSE)
        cat("  ✓ Model generated successfully\n")
        cat("  Model preview (first 500 characters):\n")
        cat("  ", substr(model_string, 1, 500), "...\n", sep = "")

        # Check for key components
        has_model <- grepl("model \\{", model_string)
        has_mvn <- grepl("dmnorm", model_string)
        has_priors <- grepl("alpha.*~ dnorm", model_string)
        has_lambda <- grepl("lambda.*~ dunif", model_string)

        cat("\n  Model checks:\n")
        cat("    - Has model block:", has_model, "\n")
        cat("    - Has multivariate normal:", has_mvn, "\n")
        cat("    - Has priors for alpha:", has_priors, "\n")
        cat("    - Has lambda parameters:", has_lambda, "\n")

        if (!all(has_model, has_mvn, has_priors, has_lambda)) {
            cat("  ⚠ Some expected components missing\n")
        }
    },
    error = function(e) {
        cat("  ✗ Error generating model:", e$message, "\n")
    }
)

# Test 4: Generate model with multi.tree = TRUE
cat("\nTest 4: Generating JAGS model with multi.tree = TRUE...\n")
tryCatch(
    {
        model_string_multi <- phybase_model(equations, multi.tree = TRUE)
        cat("  ✓ Multi-tree model generated\n")

        has_multivev <- grepl("multiVCV", model_string_multi)
        has_K <- grepl("K ~ dcat", model_string_multi)

        cat("  Multi-tree checks:\n")
        cat("    - Has multiVCV:", has_multivev, "\n")
        cat("    - Has tree sampling (K):", has_K, "\n")
    },
    error = function(e) {
        cat("  ✗ Error generating multi-tree model:", e$message, "\n")
    }
)

# Test 5: Test d-separation function
cat("\nTest 5: Testing phybase_dsep()...\n")
tryCatch(
    {
        ind_equations <- phybase_dsep(equations)
        cat("  ✓ D-separation equations generated\n")
        cat("  Number of independence tests:", length(ind_equations), "\n")
        cat("  Independence equations:\n")
        for (i in seq_along(ind_equations)) {
            cat("    ", deparse(ind_equations[[i]]), "\n")
        }
    },
    error = function(e) {
        cat("  ✗ Error in phybase_dsep():", e$message, "\n")
    }
)

# Test 6: Prepare data for JAGS
cat("\nTest 6: Preparing data for JAGS...\n")
tryCatch(
    {
        mod_data <- list(
            BM = rhino.dat$BM,
            LS = rhino.dat$LS,
            NL = rhino.dat$NL,
            DD = rhino.dat$DD,
            RS = rhino.dat$RS
        )

        # Check for missing values
        na_counts <- sapply(mod_data, function(x) sum(is.na(x)))
        cat("  ✓ Data prepared\n")
        cat("  Missing values per variable:\n")
        for (var in names(na_counts)) {
            cat("    ", var, ":", na_counts[var], "\n")
        }
    },
    error = function(e) {
        cat("  ✗ Error preparing data:", e$message, "\n")
    }
)

# Test 7: Check if JAGS is installed
cat("\nTest 7: Checking JAGS installation...\n")
jags_available <- tryCatch(
    {
        rjags::jags.version()
        TRUE
    },
    error = function(e) {
        FALSE
    }
)

if (jags_available) {
    jags_ver <- rjags::jags.version()
    cat(
        "  ✓ JAGS is installed (version:",
        paste(jags_ver, collapse = "."),
        ")\n"
    )

    # Test 8: Run a small test model (if JAGS is available)
    cat("\nTest 8: Running a small MCMC test (100 iterations)...\n")
    tryCatch(
        {
            test_result <- phybase_run(
                data = mod_data,
                tree = rhino.tree,
                equations = equations,
                n.iter = 100,
                n.burnin = 50,
                n.thin = 1,
                n.chains = 2,
                DIC = FALSE,
                quiet = TRUE
            )

            cat("  ✓ Model ran successfully\n")
            cat("  Output class:", class(test_result), "\n")
            cat(
                "  Monitored parameters:",
                paste(test_result$monitor, collapse = ", "),
                "\n"
            )

            # Check samples
            if (!is.null(test_result$samples)) {
                cat("  ✓ Posterior samples obtained\n")
                cat("  Number of chains:", length(test_result$samples), "\n")
            }
        },
        error = function(e) {
            cat("  ✗ Error running model:", e$message, "\n")
        }
    )
} else {
    cat("  ⚠ JAGS is not installed - skipping MCMC tests\n")
    cat("  Install JAGS from: http://mcmc-jags.sourceforge.net\n")
}

cat("\n=== Testing Complete ===\n")
