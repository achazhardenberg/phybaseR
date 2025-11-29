pkgname <- "phybaseR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('phybaseR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("phybase_dsep")
### * phybase_dsep

flush(stderr()); flush(stdout())

### Name: phybase_dsep
### Title: Extract d-separation statements from a structural equation model
### Aliases: phybase_dsep

### ** Examples

# Define a simple path model
equations <- list(LS ~ BM, NL ~ BM + RS, DD ~ NL)

# Get conditional independence tests
ind_tests <- phybase_dsep(equations)

# Each test is a formula
print(ind_tests)




cleanEx()
nameEx("phybase_model")
### * phybase_model

flush(stderr()); flush(stdout())

### Name: phybase_model
### Title: Generate a JAGS model string for Phylogenetic Bayesian SEM
###   (PhyBaSE)
### Aliases: phybase_model

### ** Examples

eqs <- list(BR ~ BM, S ~ BR, G ~ BR, L ~ BR)
cat(phybase_model(eqs, multi.tree = TRUE))




cleanEx()
nameEx("rhino.dat")
### * rhino.dat

flush(stderr()); flush(stdout())

### Name: rhino.dat
### Title: Rhinograd life-history data
### Aliases: rhino.dat
### Keywords: datasets

### ** Examples

data(rhino.dat)
head(rhino.dat)
summary(rhino.dat)




cleanEx()
nameEx("rhino.tree")
### * rhino.tree

flush(stderr()); flush(stdout())

### Name: rhino.tree
### Title: Rhinograd phylogenetic tree
### Aliases: rhino.tree
### Keywords: datasets

### ** Examples

data(rhino.tree)
plot(rhino.tree)




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
