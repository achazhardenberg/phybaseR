# Manual test for ZIP and ZINB distributions in 'because'
tryCatch(
    {
        library(devtools)
        load_all(".")
    },
    error = function(e) {
        cat("Error loading package: ", e$message, "\n")
        # Fallback: Validation requires dependencies, might fail if not installed
        source("R/because.R")
        source("R/because_model.R")
        # Source other helpers
        files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
        sapply(files, source)
    }
)

if (!exists("because")) {
    stop("Function 'because' not found after loading.")
}

set.seed(123)
N <- 200
x <- rnorm(N)
psi_true <- 0.3
beta_0 <- 1.5
beta_1 <- 0.5
mu <- exp(beta_0 + beta_1 * x)

# Simulate ZIP data
z <- rbinom(N, 1, psi_true) # 1 = structural zero
y_zip <- ifelse(z == 1, 0, rpois(N, mu))

cat(
    "Simulated ZIP data:\nMean:",
    mean(y_zip),
    "\nZero Prop:",
    mean(y_zip == 0),
    "\n"
)

# Run ZIP model
cat("Running ZIP model...\n")
mod_zip <- because(
    data = list(y = y_zip, x = x),
    equations = list(y ~ x),
    distribution = c(y = "zip"),
    n.iter = 1000,
    n.burnin = 500,
    n.chains = 1,
    quiet = TRUE
)

print(summary(mod_zip))

cat("\n-------------------------------------------------\n")

# Simulate ZINB data
size_true <- 2
# For ZINB: if z=1 -> 0, else NegBin
y_zinb <- ifelse(z == 1, 0, rnbinom(N, size = size_true, mu = mu))

cat(
    "Simulated ZINB data:\nMean:",
    mean(y_zinb),
    "\nZero Prop:",
    mean(y_zinb == 0),
    "\n"
)

# Run ZINB model
cat("Running ZINB model...\n")
mod_zinb <- because(
    data = list(y = y_zinb, x = x),
    equations = list(y ~ x),
    distribution = c(y = "zinb"),
    n.iter = 1000,
    n.burnin = 500,
    n.chains = 1,
    quiet = TRUE
)

print(summary(mod_zinb))
