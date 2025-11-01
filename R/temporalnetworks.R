#' temporalnetworks: Miscellaneous convenience functions
#'
#' @description Functions for estimating, analyzing, and visualizing different types of temporal networks.
#'
#' @name temporalnetworks
#' @useDynLib temporalnetworks
#' @md
#' @import magrittr dplyr lmerTest ggplot2 doParallel parallel foreach iterators
#' @import qgraph Rcpp deldir
#' @importFrom grDevices rgb col2rgb dev.off png
#' @importFrom stats as.formula ave coef cor formula logLik model.frame na.omit
#' quantile reformulate sd setNames terms update var lm resid optimize pnorm
#' @importFrom lme4 fixef lmerControl ranef
#' @importFrom igraph graph_from_adjacency_matrix E V delete_edges
#' layout_with_kk layout_with_fr layout_with_drl layout_with_dh
#' make_ring subgraph_isomorphisms
#' @importFrom graphlayouts layout_as_dynamic layout_with_stress layout_with_focus
#' @importFrom glmnet glmnet
#' @importFrom Matrix Matrix
#'
#'
"_PACKAGE"


