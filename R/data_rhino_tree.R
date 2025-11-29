#' Rhinograd phylogenetic tree
#'
#' A phylogenetic tree (phylo object) for 100 simulated Rhinograd species.
#' This tree is used to demonstrate phylogenetic Bayesian structural equation
#' modeling with the phybaseR package.
#'
#' @format A phylo object (from the ape package) with 100 tips and 99 internal nodes.
#' The tree has been scaled so that the root age is 1.0.
#'
#' @details
#' Rhinograds (nasobames) are hypothetical mammals used as examples in
#' phylogenetic comparative methods. This simulated tree represents the
#' evolutionary relationships among 100 Rhinograd species.
#'
#' @source Simulated data for von Hardenberg and Gonzalez-Voyer (2025)
#'
#' @references
#' von Hardenberg, A. and Gonzalez-Voyer, A. (2025).
#' PhyBaSE: A Bayesian approach to Phylogenetic Structural Equation Models.
#' Methods in Ecology and Evolution. https://doi.org/10.1111/2041-210X.70044
#'
#' @examples
#' data(rhino.tree)
#' plot(rhino.tree)
#'
"rhino.tree"
