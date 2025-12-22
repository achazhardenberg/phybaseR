# =============================================================================
# Custom Family Demo: Student-t Distribution for Robust Regression
# =============================================================================
# This demo shows how to create a custom distribution family using because_family()
# The Student-t distribution is useful for robust regression when data has outliers.
# =============================================================================

library(because)

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("Custom Family Demo: Student-t Robust Regression\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# =============================================================================
# STEP 1: Create the Custom Family
# =============================================================================

cat("STEP 1: Creating Student-t family with because_family()\n\n")

# Create a Student-t family - just provide the JAGS likelihood!
student_t <- because_family(
    name = "student_t",
    jags_likelihood = "{response}[{i}] ~ dt({mu}[{i}], {tau}, df_{response})",
    extra_priors = c("df_{response} ~ dunif(2, 100)"), # Degrees of freedom prior
    link = "identity",
    description = "Student-t distribution for robust regression"
)

# =============================================================================
# STEP 2: Simulate Data with Outliers
# =============================================================================

cat("\nSTEP 2: Simulating data with outliers\n\n")

set.seed(42)
n <- 50

# True parameters
beta_true <- 1.5
intercept_true <- 2.0

# Generate data
X <- rnorm(n)
Y_clean <- intercept_true + beta_true * X + rnorm(n, sd = 0.5)

# Add some outliers (10%)
n_outliers <- 5
outlier_idx <- sample(n, n_outliers)
Y <- Y_clean
Y[outlier_idx] <- Y[outlier_idx] + rnorm(n_outliers, mean = 0, sd = 5)

sim_data <- data.frame(Y = Y, X = X)

cat("  N observations:", n, "\n")
cat("  N outliers:", n_outliers, "\n")
cat("  True beta:", beta_true, "\n")
cat("  True intercept:", intercept_true, "\n\n")

# =============================================================================
# STEP 3: Fit Model with Custom Family
# =============================================================================

cat("STEP 3: Fitting model with Student-t family\n\n")

# Create the family object
student_t_family <- student_t()
cat("Family class:", paste(class(student_t_family), collapse = ", "), "\n\n")

cat("Note: Full integration is complete!\n")
cat("Fitting an actual model using the custom family...\n\n")

# Fit the model using because()
# Note: we use our custom family constructor!
fit <- because(
    Y ~ X,
    data = sim_data,
    family = student_t()
)

print(summary(fit))
cat("\nModel fit complete.\n\n")

# Show what the generated JAGS code looks like (for reference)
cat("Generated JAGS likelihood code (just for reference):\n")
likelihood_result <- jags_family_likelihood(student_t_family, "Y", suffix = "")
cat("  ", likelihood_result$likelihood_code, "\n")
cat("\nGenerated priors:\n")
cat("  ", likelihood_result$prior_code, "\n")

# =============================================================================
# CONCLUSION
# =============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("WHAT WE LEARNED\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Creating custom families is as simple as:\n\n")

cat("1. Define your JAGS likelihood:\n\n")
cat("   my_family <- because_family(\n")
cat("       name = 'my_custom',\n")
cat(
    "       jags_likelihood = '{response}[{i}] ~ <distribution>({mu}[{i}], ...)',\n"
)
cat("       extra_priors = c('param ~ prior()')\n")
cat("   )\n\n")

cat("2. Use it in any because() model:\n\n")
cat("   fam <- my_family()\n")
cat("   because(Y ~ X, data = data, family = fam)\n\n")

cat("Available placeholders:\n")
cat("  {response}  - Response variable name\n")
cat("  {mu}        - Linear predictor (mu_<response>)\n")
cat("  {tau}       - Precision parameter (tau_e_<response>)\n")
cat("  {i}         - Loop index\n\n")
