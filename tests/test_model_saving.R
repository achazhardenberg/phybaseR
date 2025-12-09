library(becauseR)

data("rhino.dat")
data("rhino.tree")

data_list <- list(
    BM = rhino.dat$BM,
    NL = rhino.dat$NL,
    DD = rhino.dat$DD,
    LS = rhino.dat$LS,
    RS = rhino.dat$RS
)

equations_1 <- list(LS ~ BM, NL ~ BM, DD ~ NL, RS ~ DD)

# Test model saving
cat("\n=== Testing Model Saving for d-sep ===\n")
fit_rhino_dsep <- because(
    data = data_list,
    tree = rhino.tree,
    equations = equations_1,
    dsep = TRUE,
    n.iter = 500,
    n.burnin = 250
)

cat("\n\n=== Checking Models Field ===\n")
cat("Number of models saved:", length(fit_rhino_dsep$models), "\n")
cat("Names of list elements:", names(fit_rhino_dsep), "\n\n")

cat("=== Model 1 (first 30 lines) ===\n")
model_1_lines <- strsplit(fit_rhino_dsep$models[[1]], "\n")[[1]]
cat(paste(model_1_lines[1:min(30, length(model_1_lines))], collapse = "\n"))

cat("\n\n=== Model 3 (first 30 lines) ===\n")
model_3_lines <- strsplit(fit_rhino_dsep$models[[3]], "\n")[[1]]
cat(paste(model_3_lines[1:min(30, length(model_3_lines))], collapse = "\n"))

cat("\n\nTest completed successfully!\n")
