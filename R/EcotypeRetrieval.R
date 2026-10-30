#' Run EcotypeRetrieval
#'
#' \code{EcotypeRetrieval.R} Retrieve states in an external dataset
#'
#' @param Statescope_StateRetrieved SummarizedExperiment object from Statescope
#'  StateRetrieval.
#' @param Statescope_EcotypeDiscovered SummarizedExperiment object from
#'  Statescope EcotypeDiscovery
#'
#' @return SummarizedExperiment object with ecotypescores
#' @import reticulate basilisk
#' @importFrom Matrix rowSums
#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom SummarizedExperiment assay
#' @importFrom methods is
#' @export
#'
#' @examples
#' ## Load StateRetrieved Statescope object
#' load(system.file("extdata", "example_Statescope_StateDiscovered.RData",
#'     package = "StatescopeR"
#' ))
#' Statescope_StateRetrieved = Statescope
#'
#' ## Load EcotypeDiscovered Statescope object
#' load(system.file("extdata", "example_Statescope_EcotypeDiscovered.RData",
#'     package = "StatescopeR"
#' ))
#'
#' Statescope_EcotypeDiscovered = Statescope
#'
#' ## Retrieve states in refined Statescope object
#' Statescope_EcotypeRetrieved <- EcotypeRetrieval(Statescope_StateRetrieved,
#'                 Statescope_EcotypeDiscovered)
#'
#' ## Look at ecotypescores in external dataset
#' S4Vectors::metadata(Statescope_EcotypeRetrieved)$ecotypescores
#'
EcotypeRetrieval <- function(Statescope_StateRetrieved,
                             Statescope_EcotypeDiscovered) {
    if (!is(Statescope_StateRetrieved, "SummarizedExperiment") |
        !is(Statescope_EcotypeDiscovered, "SummarizedExperiment")) {
        stop("Statescope objects are not SummarizedExperiment objects")}

    ## Prepare state scores and ecotype loadings
    state_scores_df <- StatescopeR::concatenate_statescores(
        S4Vectors::metadata(Statescope_StateRetrieved)$statescores)
    state_scores <- t(as.matrix(state_scores_df))

    ecotypeloadings =
        S4Vectors::metadata(Statescope_EcotypeDiscovered)$ecotypeloadings
    ## check if retrieved states are the same as in the ecotype loadings
    if (any(colnames(state_scores_df) != rownames(ecotypeloadings))){
        stop("Statescores are not of same states as EcotypeLoadings.
        Check if Statescores are actually recovered from same external dataset")
    }

    ## start basilisk & run StateRetrieval
    proc <- basiliskStart(statescope)
    Statescope_EcotypeRetrieved <- basiliskRun(proc, fun <- function(
        Statescope_StateRetrieved, state_scores, ecotypeloadings) {
        ## source cNMF code
        reticulate::source_python(system.file("python/cNMF_functions.py",
            package = "StatescopeR"))
        reticulate::source_python(system.file("python/cNMF_helper_functions.py",
            package = "StatescopeR"))

        ## Run cNMF based EcotypeRetrieval
        cNMF_model <- cNMF_Retrieval(state_scores, as.data.frame(
            ecotypeloadings))
        final_H <- t(cNMF_model$H)
        final_H_norm <- S4Vectors::DataFrame(final_H / rowSums(final_H))
        rownames(final_H_norm) <- colnames(Statescope_StateRetrieved)
        #colnames(final_H_norm) = paste0('Ecotype_', seq(ncol(final_H_norm)))

        ## Add ecotypescores and loadings to refined Statescope obj
        S4Vectors::metadata(Statescope_StateRetrieved)$ecotypescores <-
            final_H_norm
        S4Vectors::metadata(Statescope_StateRetrieved)$ecotypeloadings <-
            ecotypeloadings
        Statescope_StateRetrieved},
        Statescope_StateRetrieved = Statescope_StateRetrieved,
        state_scores = state_scores, ecotypeloadings= ecotypeloadings)

    ## stop basilisk
    basiliskStop(proc)

    return(Statescope_EcotypeRetrieved)
}
