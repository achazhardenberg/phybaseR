#' Rhinograd phylogenetic tree
#'
#' A phylogenetic tree (phylo object) for 100 simulated Rhinograd species.
#' This tree is used to demonstrate phylogenetic Bayesian structural equation
#' modeling with the becauseR package.
#'
#' @format A phylo object (from the ape package) with 100 tips and 99 internal nodes.
#' The tree has been scaled so that the root age is 1.0.
#'
#' @details
#' Rhinograds  are hypothetical mammals used as examples in
#' phylogenetic comparative methods. This simulated tree represents the
#' evolutionary relationships among 100 Rhinograd species.
#'
#' @source Simulated data for Gonzalez-Voyer and von Hardenberg (2014)
#'
#' @references
#' Gonzalez-Voyer, A., & von Hardenberg, A. (2014). An introduction to
#' phylogenetic path analysis. In: Modern phylogenetic comparative methods
#' and their application in evolutionary biology (pp. 201–229). Springer.
#'
#' @examples
#' data(rhino.tree)
#' plot(rhino.tree)
#'
"rhino.tree"
