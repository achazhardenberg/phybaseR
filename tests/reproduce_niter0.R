library(because)

set.seed(123)
N <- 100
df <- data.frame(
    x = rnorm(N),
    y = rnorm(N)
)

message("Attempting to run because() with n.iter = 0...")

tryCatch(
    {
        fit <- because(
            equations = list(y ~ x),
            data = df,
            n.iter = 0,
            n.adapt = 0,
            quiet = TRUE
        )
        message("SUCCESS: because() completed without error.")
    },
    error = function(e) {
        message("FAILURE: because() failed with error:")
        message(e$message)
    }
)
