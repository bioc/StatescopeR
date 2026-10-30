## usethis namespace: start
## usethis namespace: end
NULL

#' The StatescopeR package
#'
#' framework for discovery of cell states from bulk mRNA profiles
#'
#' @rdname StatescopeR-package
#' @name StatescopeR-package
#' @keywords internal
#' @aliases StatescopeR-package StatescopeR
#'
#' @section Details:
#'  StatescopeR starts from lognormalized cp10k single cell data and a bulk RNA
#'  dataset to be used for deconvolution and cell state analysis. From this
#'  point StatescopeR has some functions to create a signature and select genes
#'  for deconvolution.
#'  \code{\link{create_signature}} creates a signature of lognormalized cp10k
#'      single cell data for deconvolution and cell state analysis.
#'  \code{\link{select_genes}} selects genes which best discriminate cell types
#'      for use in deconvolution.
#'  Alternatively you can also use:
#'   \code{\link{fetch_signature}} to fetch one of the premade signatures +
#'   selected_genes from https://github.com/tgac-vumc/StatescopeData
#'   After these steps and proper cp10k normalization of the bulk mRNA
#'   deconvolution, cell state and ecotype analysis can be done with the
#'    following functions:
#'  \code{\link{BLADE_deconvolution}} estimates cell fractions from bulk mRNA.
#'  \code{\link{Refinement}} refines cell type-specific gene expression profile
#'      estimates to better capture inter-sample variability.
#'  \code{\link{StateDiscovery}} discovery cell states from inferred cell
#'  type-specific gene expression profiles.
#'  \code{\link{StateRetrieval}} retrieve cell states in an external data set.
#'  \code{\link{EcotypeDiscovery}} discover ecotypes from discovered cell state
#'  scores
#'  \code{\link{EcotypeRetrieval}} retrieve ecotypes in an external data set
#'
#'  After these steps the following functions provide some evaluation
#'      and plotting functions:
#'  \code{\link{gather_true_fractions}} gets true cell fractions from
#'      scRNAseq data.
#'  \code{\link{fraction_eval}} plots the correlation of estimated vs true
#'      cell fractionsover all samples.
#'  \code{\link{fraction_heatmap}} plots a heatmap of estimated
#'      cell fractions.
#'  \code{\link{barplot_stateloadings}} plots a barplot showing the top n
#'      most important stateloadings.
#'

"_PACKAGE"
