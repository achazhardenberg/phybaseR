# DEMO: Bayesian Model Averaging (BMA) with 'because' + JAGS
# -------------------------------------------------------------
# This script demonstrates how to implement the "Indicator Variable" approach
# for BMA on the 'sem8' Rhino model.
#
# Steps:
# 1. Use because() to generate the base JAGS code and formatted data.
# 2. Modify the JAGS code to inject binary switches (z) for selected paths.
# 3. Run the modified model in JAGS.
# 4. Calculate Posterior Inclusion Probabilities (PIPs).

library(because)
library(rjags)

# 1. Load Data
# ------------
data(rhino.dat)
data(rhino.tree)

# The "Best" Model (sem8)
sem8_eq <- list(
    LS ~ BM,
    NL ~ BM + RS,
    DD ~ NL
)

# 2. Generate Base Code (The "Scaffold")
# --------------------------------------
# We run with n.iter=0 just to get the code and data list
cat("Generating base model...\n")
fit_base <- because(
    equations = sem8_eq,
    data = rhino.dat,
    structure = rhino.tree,
    id_col = "SP",
    n.iter = 0,
    quiet = TRUE
)

# Extract correct data list (phylogenetic matrix already inverted!)
jags_data <- fit_base$data
base_code <- fit_base$model_code

# 3. Inject BMA Logic (The "Surgery")
# -----------------------------------
# We want to test if 'NL -> DD' is real.
# The base code has:  mu_DD[i] <- alpha_DD + beta_DD_NL * NL[i]
# We change it to:    mu_DD[i] <- alpha_DD + (z_DD_NL * beta_DD_NL) * NL[i]

# Add prior for z
bma_prior <- "
  # BMA Indicator Priors
  z_DD_NL ~ dbern(0.5)
  z_NL_RS ~ dbern(0.5)
  z_LS_BM ~ dbern(0.5)
"

# Inject into code (Robust Regex replacement)
modified_code <- base_code

# Helper to inject Z
inject_z <- function(code, beta_name, z_name) {
    # Regex: find beta_name followed by optional space and *
    pattern <- paste0("(", beta_name, ")\\s*\\*")
    replacement <- paste0("(", z_name, " * \\1) *")

    new_code <- gsub(pattern, replacement, code)

    if (new_code == code) {
        warning(paste("Failed to inject", z_name, "for", beta_name))
    }
    return(new_code)
}

# 1. DD ~ NL
modified_code <- inject_z(modified_code, "beta_DD_NL", "z_DD_NL")
# 2. NL ~ RS
modified_code <- inject_z(modified_code, "beta_NL_RS", "z_NL_RS")
# 3. LS ~ BM
modified_code <- inject_z(modified_code, "beta_LS_BM", "z_LS_BM")

# Add the z priors before the likelihood loop
modified_code <- sub("model \\{", paste("model {\n", bma_prior), modified_code)

# 3b. Add Monitoring for Effective Coefficients
# ---------------------------------------------
# When z=0, the raw 'beta' disconnects from data and wanders the wide prior (SD=1000).
# This makes the summary stats look terrible.
# We must monitor the EFFECTIVE beta (z * beta), which is exactly 0 when z=0.

# Add deterministic nodes to track effective betas
calc_eff_code <- "
  # Effective Coefficients (Monitor These!)
  beta_eff_DD_NL <- z_DD_NL * beta_DD_NL
  beta_eff_NL_RS <- z_NL_RS * beta_NL_RS
  beta_eff_LS_BM <- z_LS_BM * beta_LS_BM
}
"
# Replace the closing brace with the new nodes + closing brace
# Replace the LAST closing brace with the new nodes + closing brace
# (Safety: sub() matches the first brace, which might be a for-loop closing brace.
# We must find the absolute last '}' in the string.)

matches <- gregexpr("\\}", modified_code)[[1]]
if (length(matches) > 0 && matches[1] != -1) {
    last_pos <- matches[length(matches)]

    # Inject before the last '}'
    # The replacement code (calc_eff_code) ends with '}'
    # So we replace the character at last_pos with calc_eff_code

    prefix <- substr(modified_code, 1, last_pos - 1)
    suffix <- substr(modified_code, last_pos + 1, nchar(modified_code))

    modified_code <- paste0(prefix, calc_eff_code, suffix)
} else {
    stop("Could not find closing brace '}' in generated JAGS code.")
}

cat("Modified JAGS Code for BMA:\n")
cat(substr(modified_code, 1, 1000), "...\n")

# 4. Run the BMA Model
# --------------------
cat("\nCompiling BMA model...\n")
model_bma <- suppressWarnings(jags.model(
    textConnection(modified_code),
    data = jags_data,
    n.chains = 3,
    quiet = TRUE
))

cat("Sampling (this may take longer due to mixing)....\n")
update(model_bma, n.iter = 1000) # Burn-in
samples_bma <- coda.samples(
    model_bma,
    variable.names = c(
        "z_DD_NL",
        "z_NL_RS",
        "z_LS_BM",
        "beta_eff_DD_NL",
        "beta_eff_NL_RS",
        "beta_eff_LS_BM"
    ),
    n.iter = 5000
)

# 5. Results: Posterior Inclusion Probabilities
# ---------------------------------------------
cat("\n--- BMA RESULTS (Effective Coefficients) ---\n")
summary_bma <- summary(samples_bma)
print(summary_bma$statistics[, c("Mean", "SD")])
summary_bma
# Interpretation:
# Mean(z) = Posterior Inclusion Probability (PIP)
# If Mean(z_DD_NL) = 0.95, we are 95% sure the arrow exists.
# If Mean(z_DD_NL) = 0.40, the evidence is weak.

# Save Plots
pdf("bma_demo_plots.pdf")
plot(samples_bma)
dev.off()
cat("\nPlots saved to 'bma_demo_plots.pdf'\n")
