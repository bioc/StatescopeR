#' Run StateRetrieval
#'
#' \code{StateRetrieval.R} Retrieve states in an external dataset
#'
#' @param Statescope_Refined SummarizedExperiment object from Statescope Refinement.
#' @param Statescope_StateDiscovered SummarizedExperiment object from Statescope StateDiscovery
#'
#' @return SummarizedExperiment object with statescores per celltype added
#' @import reticulate basilisk
#' @importFrom Matrix rowSums
#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom SummarizedExperiment assay
#' @importFrom methods is
#' @export
#'
#' @examples
#' ## Load Refined Statescope object
#' load(system.file("extdata", "example_Statescope_Refined.RData",
#'     package = "StatescopeR"
#' ))
#' Statescope_Refined = Statescope
#'
#' ## Load StateDiscovered Statescope object
#' load(system.file("extdata", "example_Statescope_StateDiscovered.RData",
#'     package = "StatescopeR"
#' ))
#'
#' Statescope_StateDiscovered = Statescope
#'
#' ## Retrieve states in refined Statescope object
#' Statescope_StateRetrieved <-
#'                StateRetrieval(Statescope_Refined, Statescope_StateDiscovered)
#'
#' ## Look at statescores in external dataset
#' S4Vectors::metadata(Statescope_StateRetrieved)$statescores
#'
StateRetrieval <- function(Statescope_Refined, Statescope_StateDiscovered) {
    if (!is(Statescope_Refined, "SummarizedExperiment") |
        !is(Statescope_StateDiscovered, "SummarizedExperiment")) { # Check input
        stop("Statescope object are not SummarizedExperiment objects")}

    ## Find overlapping genes between Statescope objects
    genes = intersect(rownames(assay(metadata(
        Statescope_Refined)$ct_specific_gep[[1]], "weighted_gep")),
        rownames(metadata(Statescope_StateDiscovered)$stateloadings[[1]]))

    ## start basilisk & run StateRetrieval
    proc <- basiliskStart(statescope)
    Statescope_StateRetrieved <- basiliskRun(proc, fun <- function(
        Statescope_Refined, Statescope_StateDiscovered, genes) {
        ## source cNMF code
        reticulate::source_python(system.file("python/cNMF_functions.py",
            package = "StatescopeR"))
        reticulate::source_python(system.file("python/cNMF_helper_functions.py",
            package = "StatescopeR"))

        ## Run cNMF StateRetrieval for each cell type
        statescores <- list()
        for (ct in
             names(S4Vectors::metadata(Statescope_Refined)$ct_specific_gep)) {
            ## get ct_specific_gep for state retrieval
            data_scaled <- as.matrix(SummarizedExperiment::assay(
                S4Vectors::metadata(Statescope_Refined)$ct_specific_gep[[ct]],
                "weighted_gep"))
            StateLoadings = as.data.frame(S4Vectors::metadata( # df for python
                Statescope_StateDiscovered)$stateloadings[[ct]])
            ## Run state retrieval
            cNMF_model <-
                cNMF_Retrieval(data_scaled[genes,], StateLoadings[genes,])
            final_H <- t(cNMF_model$H)
            final_H_norm <- S4Vectors::DataFrame(final_H / rowSums(final_H))
            rownames(final_H_norm) <- colnames(Statescope_Refined)

            ## Add result to ct lists
            statescores[[ct]] <- final_H_norm}

        ## Add statescores and loadings to refined Statescope obj
        S4Vectors::metadata(Statescope_Refined)$statescores <- statescores
        S4Vectors::metadata(Statescope_Refined)$stateloadings <-
            S4Vectors::metadata(Statescope_StateDiscovered)$stateloadings
        Statescope_Refined}, Statescope_Refined = Statescope_Refined,
        Statescope_StateDiscovered = Statescope_StateDiscovered, genes = genes)
    ## stop basilisk
    basiliskStop(proc)

    return(Statescope_StateRetrieved)
}
