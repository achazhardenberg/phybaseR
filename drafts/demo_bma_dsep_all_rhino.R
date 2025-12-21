# DEMO: Comprehensive BMA D-separation Test
# -----------------------------------------
# Goal: Test ALL independence claims by looping through missing links.

library(because)
library(rjags)

# 1. Load Data
data(rhino.dat)
data(rhino.tree)

# Base Model (sem8)
sem8_eq <- list(
    LS ~ BM,
    NL ~ BM + RS,
    DD ~ NL
)

# 2. Define the Test Cases (The "Basis Set" of missing links)
# These are the pairs that SHOULD be independent.
# Format: Response ~ Predictor
# IMPORTANT: Must respect topological order (BM, RS -> NL -> DD) to avoid cycles!
test_links <- list(
    c("RS", "BM"), # Test: RS ~ BM (BM->RS)
    c("DD", "RS"), # Test: DD ~ RS (RS->DD) -- Fixed direction
    c("LS", "RS"), # Test: LS ~ RS (RS->LS) -- Fixed direction
    c("DD", "BM"), # Test: DD ~ BM (BM->DD) -- Fixed direction
    c("NL", "LS"), # Test: NL ~ LS (LS->NL)
    c("DD", "LS") # Test: DD ~ LS (LS->DD)
)

results_table <- data.frame(
    Test_Link = character(),
    PIP = numeric(),
    Prob_Indep = numeric(),
    stringsAsFactors = FALSE
)

# Helper for Robust Injection
inject_bma <- function(code, beta_name, z_name) {
    # 1. Add Z prior
    code <- sub(
        "model \\{",
        paste0("model {\n  ", z_name, " ~ dbern(0.5)\n"),
        code
    )

    # 2. Replace beta in structural equations
    # STRICT Regex: beta_name followed by optional space then *
    # We replace "beta *" with "beta_eff *"
    # This avoids matching "beta ~ dnorm"
    beta_eff_name <- paste0("beta_eff_", beta_name)
    pattern <- paste0("(", beta_name, ")\\s*\\*")
    replacement <- paste0(beta_eff_name, " *")
    code <- gsub(pattern, replacement, code)

    # 3. Define beta_eff node at the end
    # Important: Must interpret z*beta as effective coef
    calc_eff <- paste0(
        "  ",
        beta_eff_name,
        " <- ",
        z_name,
        " * ",
        beta_name,
        "\n}\n"
    )

    # Find last brace
    matches <- gregexpr("\\}", code)[[1]]
    last_pos <- matches[length(matches)]
    prefix <- substr(code, 1, last_pos - 1)
    suffix <- substr(code, last_pos + 1, nchar(code))
    code <- paste0(prefix, calc_eff, suffix)

    return(code)
}

# Helper to safe-merge equations (Prevent duplicate definitions)
add_term_to_eqs <- function(eq_list, resp, pred) {
    # Convert to character (lhs ~ rhs)
    # Find if LHS already exists
    for (i in seq_along(eq_list)) {
        f <- eq_list[[i]]
        cols <- as.character(f) # [1]="~", [2]=LHS, [3]=RHS
        if (cols[2] == resp) {
            # Merge: LHS ~ RHS + pred
            new_rhs <- paste0(cols[3], " + ", pred)
            eq_list[[i]] <- as.formula(paste(resp, "~", new_rhs))
            return(eq_list)
        }
    }
    # Not found? Add new.
    eq_list <- c(eq_list, list(as.formula(paste(resp, "~", pred))))
    return(eq_list)
}

# 3. Loop and Test
cat("Running BMA on 6 Independence Claims...\n")

for (link in test_links) {
    resp <- link[1]
    pred <- link[2]
    link_str <- paste(resp, "~", pred)
    cat(sprintf("\nTesting: %s ... ", link_str))

    # Construct Test Model (Base + Forbidden Link)
    # Safely merge
    test_eq_curr <- add_term_to_eqs(sem8_eq, resp, pred)

    # Generate Code
    fit <- because(
        equations = test_eq_curr,
        data = rhino.dat,
        structure = rhino.tree,
        id_col = "SP",
        n.iter = 0,
        quiet = TRUE
    )

    # Identify Parameter Name
    # because() typically names it beta_RESP_PRED
    beta_name <- paste0("beta_", resp, "_", pred)
    z_name <- paste0("z_", resp, "_", pred)

    # Check if parameter exists in code
    if (!grepl(beta_name, fit$model_code)) {
        # Try reverse? No, input order defines name
        cat(" [Error: Param not found] \n")
        next
    }

    # Inject BMA
    mod_code <- inject_bma(fit$model_code, beta_name, z_name)

    # Run JAGS
    mod_jags <- suppressWarnings(jags.model(
        textConnection(mod_code),
        data = fit$data,
        n.chains = 3,
        quiet = TRUE
    ))
    update(mod_jags, n.iter = 500)
    samps <- coda.samples(mod_jags, z_name, n.iter = 2000)

    # Get PIP
    pip <- summary(samps)$statistics["Mean"]
    prob_indep <- 1 - pip

    results_table <- rbind(
        results_table,
        data.frame(
            Test_Link = link_str,
            PIP = pip,
            Prob_Indep = prob_indep
        )
    )

    cat(sprintf("PIP=%.3f (P(Indep)=%.3f)\n", pip, prob_indep))
}

# 4. Final Report
cat("\n--- COMPREHENSIVE BMA D-SEPARATION REPORT ---\n")
print(results_table)
cat(
    "\nInterpretation: High Prob_Indep (>0.5) confirms the model structure is valid.\n"
)
