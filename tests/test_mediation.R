# Test Automated Mediation Verification

# Create a clear mediation structure:
# X -> M (beta = 0.8)
# M -> Y (beta = 0.5)
# X -> Y (beta = 0.2) [Direct]
# Expected Indirect = 0.8 * 0.5 = 0.4
# Expected Total = 0.4 + 0.2 = 0.6

set.seed(123)
N <- 200
X <- rnorm(N)
M <- 0.8 * X + rnorm(N, sd = 0.5)
Y <- 0.5 * M + 0.2 * X + rnorm(N, sd = 0.5)

df <- data.frame(X = X, M = M, Y = Y)

# Run because()
library(because)
# Need to source local since not installed
devtools::load_all()

cat("Fitting Mediation Model...\n")
fit <- because(
    list(
        M ~ X,
        Y ~ M + X
    ),
    data = df,
    n.iter = 2000
)

# Run Mediation Analysis
cat("\nRunning Mediation Analysis form X to Y...\n")
med <- because_mediation(fit, exposure = "X", outcome = "Y")

print(med$summary)
print(med$paths)

# Verification Logic
direct_est <- med$summary$Mean[med$summary$Type == "Direct Effect"]
indirect_est <- med$summary$Mean[med$summary$Type == "Total Indirect Effect"]

cat("\n--- Verification ---\n")
cat("Expected Direct: 0.2. Estimated:", round(direct_est, 3), "\n")
cat("Expected Indirect: 0.4. Estimated:", round(indirect_est, 3), "\n")

if (abs(direct_est - 0.2) < 0.1 && abs(indirect_est - 0.4) < 0.1) {
    cat("\nSUCCESS: Estimates are within valid range.\n")
} else {
    cat("\nFAILURE: Estimates deviate too far from truth.\n")
}
