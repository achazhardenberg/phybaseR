pkgname <- "phybaseR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "phybaseR-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('phybaseR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("phybase_compare")
### * phybase_compare

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: phybase_compare
### Title: Compare PhyBaSE Models
### Aliases: phybase_compare

### ** Examples

## Not run: 
##D   # Mode 1: Compare existing fits
##D   phybase_compare(fit1, fit2)
##D 
##D   # Mode 2: Run and compare
##D   specs <- list(m1 = list(equations = list(Y ~ X)), m2 = list(equations = list(Y ~ X + Z)))
##D   res <- phybase_compare(specs, data = df, tree = tr, n.cores = 2)
##D   print(res$comparison)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("phybase_compare", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("phybase_dsep")
### * phybase_dsep

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: phybase_dsep
### Title: Extract d-separation statements from a structural equation model
### Aliases: phybase_dsep

### ** Examples

# Standard DAG
equations <- list(LS ~ BM, NL ~ BM + RS, DD ~ NL)
ind_tests <- phybase_dsep(equations)

# With latent variable
equations_latent <- list(X ~ Quality, Y ~ Quality)
result <- phybase_dsep(equations_latent, latent = "Quality")
# result$tests: m-separation tests
# result$correlations: induced correlation between X and Y



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("phybase_dsep", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("phybase_format_data")
### * phybase_format_data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: phybase_format_data
### Title: Format Data for PhyBaSE Analysis
### Aliases: phybase_format_data

### ** Examples

## Not run: 
##D # Example data in long format
##D data_long <- data.frame(
##D   SP = c("sp1", "sp1", "sp1", "sp2", "sp2", "sp3"),
##D   BM = c(1.2, 1.3, 1.1, 2.1, 2.2, 1.8),
##D   NL = c(0.5, 0.6, NA, 0.7, 0.8, 0.9)
##D )
##D 
##D tree <- ape::read.tree(text = "(sp1:1,sp2:1,sp3:1);")
##D data_list <- phybase_format_data(data_long, species_col = "SP", tree = tree)
##D 
##D # Use with phybase_run
##D fit <- phybase_run(data = data_list, tree = tree, equations = list(NL ~ BM))
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("phybase_format_data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("phybase_loo")
### * phybase_loo

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: phybase_loo
### Title: Calculate LOO-CV for a PhyBaSE Model
### Aliases: phybase_loo

### ** Examples

## Not run: 
##D   fit <- phybase_run(data, tree, equations)
##D   loo_result <- phybase_loo(fit)
##D   print(loo_result)
##D 
##D   # Check for problematic observations
##D   plot(loo_result)
##D 
##D   # Compare models
##D   loo_compare(loo_result1, loo_result2)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("phybase_loo", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("phybase_loo_compare")
### * phybase_loo_compare

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: phybase_loo_compare
### Title: Compare Models Using LOO-CV
### Aliases: phybase_loo_compare

### ** Examples

## Not run: 
##D   loo1 <- phybase_loo(fit1)
##D   loo2 <- phybase_loo(fit2)
##D   phybase_loo_compare(loo1, loo2)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("phybase_loo_compare", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("phybase_model")
### * phybase_model

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: phybase_model
### Title: Generate a JAGS model string for Phylogenetic Bayesian SEM
###   (PhyBaSE)
### Aliases: phybase_model

### ** Examples

eqs <- list(BR ~ BM, S ~ BR, G ~ BR, L ~ BR)
cat(phybase_model(eqs, multi.tree = TRUE)$model)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("phybase_model", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("phybase_waic")
### * phybase_waic

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: phybase_waic
### Title: Calculate WAIC with Standard Errors for a PhyBaSE Model
### Aliases: phybase_waic

### ** Examples

## Not run: 
##D   # Fit model with WAIC monitoring
##D   fit <- phybase_run(data, tree, equations, WAIC = TRUE)
##D 
##D   # View WAIC with standard errors
##D   fit$WAIC
##D   #             Estimate   SE
##D   # elpd_waic   -617.3   12.4
##D   # p_waic        12.3    3.1
##D   # waic        1234.5   24.8
##D 
##D   # Compare two models
##D   fit1$WAIC
##D   fit2$WAIC
##D   # Model with lower WAIC is preferred
##D   # Difference is significant if |WAIC1 - WAIC2| > 2 * sqrt(SE1^2 + SE2^2)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("phybase_waic", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("rhino.dat")
### * rhino.dat

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: rhino.dat
### Title: Rhinograd life-history data
### Aliases: rhino.dat
### Keywords: datasets

### ** Examples

data(rhino.dat)
head(rhino.dat)
summary(rhino.dat)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("rhino.dat", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("rhino.tree")
### * rhino.tree

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: rhino.tree
### Title: Rhinograd phylogenetic tree
### Aliases: rhino.tree
### Keywords: datasets

### ** Examples

data(rhino.tree)
plot(rhino.tree)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("rhino.tree", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
