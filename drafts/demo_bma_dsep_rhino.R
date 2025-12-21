# DEMO: The "Ultimate D-separation Test" with BMA
# ------------------------------------------------
# Goal: Test the conditional independence claim: "LS _||_ NL | BM"
# (Litter Size indep of Neonatal Mass given Body Mass)
#
# Method:
# 1. Take the valid model (sem8).
# 2. Add the "forbidden" link: NL ~ LS.
# 3. Use BMA (z_NL_LS) to calculate the probability this link exists.
# 4. P(Independence) = 1 - Mean(z_NL_LS).

library(because)
library(rjags)

# 1. Load Data
data(rhino.dat)
data(rhino.tree)

# 2. Define Supermodel (sem8 + forbidden link)
# -------------------------------------------
# Valid: NL ~ BM + RS
# Test:  NL ~ BM + RS + LS  <-- Added LS to test dependence
test_eq <- list(
    LS ~ BM,
    NL ~ BM + RS + LS,
    DD ~ NL
)

# 3. Generate Base Code
cat("Generating base model with forbidden link...\n")
fit_base <- because(
    equations = test_eq,
    data = rhino.dat,
    structure = rhino.tree,
    id_col = "SP",
    n.iter = 0,
    quiet = TRUE
)

jags_data <- fit_base$data
base_code <- fit_base$model_code

# 4. Inject BMA Logic for the Forbidden Link ONLY
# -----------------------------------------------
# We only care about z_NL_LS. The other links are assumed true for this test.

bma_prior <- "
  # BMA Prior for the Test Link
  z_NL_LS ~ dbern(0.5)
"

modified_code <- base_code

# Inject z for beta_NL_LS
# Regex: find beta_NL_LS followed by optional space and *
pattern <- "(beta_NL_LS)\\s*\\*"
replacement <- "(z_NL_LS * \\1) *"
modified_code <- gsub(pattern, replacement, modified_code)

# Add Prior
modified_code <- sub("model \\{", paste("model {\n", bma_prior), modified_code)

# Add Monitoring for Effective Coefficient
calc_eff_code <- "
  # Effective Coefficient
  beta_eff_NL_LS <- z_NL_LS * beta_NL_LS
}
"
# Robust injection at end
matches <- gregexpr("\\}", modified_code)[[1]]
if (length(matches) > 0 && matches[1] != -1) {
    last_pos <- matches[length(matches)]
    prefix <- substr(modified_code, 1, last_pos - 1)
    suffix <- substr(modified_code, last_pos + 1, nchar(modified_code))
    modified_code <- paste0(prefix, calc_eff_code, suffix)
}

# 5. Run Test
# -----------
cat("\nCompiling D-sep Test Model...\n")
model_test <- suppressWarnings(jags.model(
    textConnection(modified_code),
    data = jags_data,
    n.chains = 3,
    quiet = TRUE
))

cat("Sampling...\n")
update(model_test, n.iter = 1000)
samples_test <- coda.samples(
    model_test,
    variable.names = c("z_NL_LS", "beta_eff_NL_LS"),
    n.iter = 10000 # Higher iterations for precise probability
)

# 6. Results
# ----------
cat("\n--- D-SEPARATION TEST RESULTS ---\n")
st <- summary(samples_test)$statistics
print(st[, c("Mean", "SD")])

pip <- st["z_NL_LS", "Mean"]
prob_indep <- 1 - pip

cat("\n------------------------------------------------\n")
cat(sprintf("Posterior Inclusion Probability (PIP): %.3f\n", pip))
cat(sprintf("Probability of Independence (1 - PIP): %.3f\n", prob_indep))
cat("------------------------------------------------\n")
if (prob_indep > 0.90) {
    cat("CONCLUSION: Strong evidence for Independence (D-sep holds).\n")
} else if (prob_indep > 0.50) {
    cat("CONCLUSION: Weak evidence for Independence.\n")
} else {
    cat("CONCLUSION: Independence rejected! Link likely exists.\n")
}
