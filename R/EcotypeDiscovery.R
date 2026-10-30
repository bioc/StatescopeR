#' Run EcotypeDiscovery
#'
#' \code{EcotypeDiscovery.R} Discovers ecotypes from cell state scores
#'
#' @param Statescope SummarizedExperiment object from Statescope StateDiscovery.
#' @param k number of cluster to choose, default is NA for automatic selection
#' @param Ncores number of cores to use for paralellization.
#' @param max_clusters maximum allowed states per cell type.
#' @param n_iter Number of initial cNMF restarts.
#' @param n_final_iter Number of final cNMF restarts.
#' @param min_cophenetic Minimum cophenetic coefficient to determine K.
#'
#' @return SummarizedExperiment object with ecotypescores and loadings
#' @import reticulate basilisk
#' @importFrom Matrix rowSums
#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom SummarizedExperiment assay
#' @importFrom methods is
#' @importFrom StatescopeR concatenate_statescores
#' @export
#'
#' @examples
#' ## Load StateDiscovered Statescope object
#' load(system.file("extdata", "example_Statescope_StateDiscovered.RData",
#'     package = "StatescopeR"
#' ))
#'
#'#' ## Discover ecotypes
#' Statescope <- EcotypeDiscovery(Statescope, k = 2L, Ncores = 2L)
#'
#' ## Look at ecotypescores and ecotypeloadings
#' S4Vectors::metadata(Statescope)$ecotypescores
#' S4Vectors::metadata(Statescope)$ecotypeloadings
#'
EcotypeDiscovery <- function(Statescope, k = NA, max_clusters = 10L, n_iter = 10L,
    n_final_iter = 100L, min_cophenetic = 0.9, Ncores = 1L) {
    if (!is(Statescope, "SummarizedExperiment")) { ## Check Statescope input
        stop("Statescope is not a SummarizedExperiment object")}

    ## start basilisk & run StateDiscovery
    proc <- basiliskStart(statescope)
    Statescope <- basiliskRun(proc, fun <- function(Statescope, max_clusters,
        n_iter, n_final_iter, min_cophenetic, Ncores) {
        ## source cNMF code
        reticulate::source_python(system.file("python/cNMF_functions.py",
            package = "StatescopeR"))
        reticulate::source_python(system.file("python/cNMF_helper_functions.py",
            package = "StatescopeR"))

        ## prepare concatenated statescore matrix
        state_scores_df <- StatescopeR::concatenate_statescores(
            S4Vectors::metadata(Statescope)$statescores)
        state_scores <- t(as.matrix(state_scores_df))

        ## perform ecotype discovery
        ## Run initial NMF runs for k selection
        if (is.na(k)) {nclust <- select_k(state_scores, max_clusters, n_iter,
                Ncores, min_cophenetic)}
        else {nclust <- k} # or use preselected k
        ## Run final model
        final_cNMF_result <- cNMF(state_scores, as.integer(nclust),
            n_final_iter, Ncores)
        final_H <- t(final_cNMF_result[[1]]$H)
        final_H_norm <- S4Vectors::DataFrame(final_H / rowSums(final_H))
        rownames(final_H_norm) <- colnames(Statescope)
        colnames(final_H_norm) = paste0('Ecotype_', seq(ncol(final_H_norm)))

        final_W <- S4Vectors::DataFrame(final_cNMF_result[[1]]$W)
        rownames(final_W) <- rownames(state_scores)
        colnames(final_W) = paste0('Ecotype_', seq(ncol(final_W)))

        ## Add ecotypescores and loadings to Statescope obj
        S4Vectors::metadata(Statescope)$ecotypescores <- final_H_norm
        S4Vectors::metadata(Statescope)$ecotypeloadings <- final_W
        Statescope},Statescope = Statescope, max_clusters = max_clusters,
    n_iter = n_iter, n_final_iter = n_final_iter,
    min_cophenetic = min_cophenetic, Ncores = Ncores)
    ## stop basilisk
    basiliskStop(proc)

    return(Statescope)
}
