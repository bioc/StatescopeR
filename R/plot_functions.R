#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @import ComplexHeatmap
#' @importFrom stats cor setNames
#' @importFrom utils head
#'
NULL

#'
#' Create a barplot of TRUE vs estimated cfs
#'
#' Create a barplot of TRUE vs estimated cfs
#'
#' @param Statescope SummarizedExperiment object from BLADE_deconvolution
#' @param true_fractions s4 Dataframe with true fractions to compare with
#'  estimated fractions
#' @return A barplot rendered to the active graphics device
#' @importFrom S4Vectors metadata
#' @importFrom methods is
#' @export fraction_eval
#' @examples
#' ## Load True fractions
#' load(system.file("extdata", "example_true_fractions.RData",
#'     package = "StatescopeR"
#' ))
#'
#' ## Load Deconvolved Statescope object
#' load(system.file("extdata", "example_Statescope_Deconvolved.RData",
#'     package = "StatescopeR"
#' ))
#'
#' ## ## Plot fraction correlation and RMSE per ct
#' fraction_eval(Statescope, true_fractions)
fraction_eval <- function(Statescope, true_fractions) {
    if (!is(Statescope, "SummarizedExperiment")) { ## Check Statescope input
        stop("Statescope is not a SummarizedExperiment object")}
    ## measure correlation and RMSE per celltype
    eval_results <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),
        c("celltype", "correlation", "RMSE"))
    for (ct in unique(rownames(true_fractions))) {
        cor_ct <- cor(as.matrix(true_fractions)[ct, ], as.matrix(
            metadata(Statescope)$fractions)[ct, names(true_fractions)])

        rmse_ct <- sqrt(mean((as.matrix(
            metadata(Statescope)$fractions)[ct, names(true_fractions)]
            - as.matrix(true_fractions)[ct, ])^2))

        ## add cor to cors
        eval_results[nrow(eval_results) + 1, ] <-
            data.frame(ct, cor_ct, rmse_ct)
    }

    ## plot correlation per celltype
    corplot <- ggplot(eval_results, aes(
        x = celltype, y = correlation,
        fill = celltype)) +
        geom_bar(stat = "identity", width = 0.95) +
        theme_bw() +
        labs(x = NULL) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

    ## plot correlation per celltype
    rmseplot <- ggplot(eval_results, aes(
        x = celltype, y = RMSE,
        fill = celltype)) +
        geom_bar(stat = "identity", width = 0.95) +
        theme_bw() +
        labs(x = NULL) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

    ## Plot fraction evaluation per ct
    plot_grid(corplot, rmseplot)
}

#' Create a heatmap of the estimated fractions
#'
#' Create a heatmap of the estimated fractions
#'
#' @param Statescope SummarizedExperiment object from BLADE_deconvolution
#' @param ... other parameters to pass to  [ComplexHeatmap::Heatmap()]
#' @return A heatmap rendered to the active graphics device
#' @importFrom S4Vectors metadata
#' @importFrom methods is
#' @export fraction_heatmap
#' @examples
#' ## Load Deconvolved Statescope object
#' load(system.file("extdata", "example_Statescope_Deconvolved.RData",
#'     package = "StatescopeR"
#' ))
#'
#' ## Plot fraction heatmap
#' fraction_heatmap(Statescope)
#'
fraction_heatmap <- function(Statescope, ...) {
    if (!is(Statescope, "SummarizedExperiment")) { ## Check Statescope input
        stop("Statescope is not a SummarizedExperiment object")
    }
    Heatmap(as.matrix(metadata(Statescope)$fractions),
        heatmap_legend_param = list(title = "")
    )
}


#' Create a barplot of top stateloadings
#'
#' Create a barplot of top stateloadings
#'
#' @param Statescope SummarizedExperiment object from StateDiscovery
#' @param top_n integer selecting how many genes to show per state
#' @return A barplot rendered to the active graphics device
#' @importFrom S4Vectors metadata
#' @importFrom methods is
#' @export barplot_stateloadings
#' @examples
#' #' ## Load StateDiscovered Statescope object
#' load(system.file("extdata", "example_Statescope_StateDiscovered.RData",
#'     package = "StatescopeR"
#' ))
#'
#' ## Plot fraction heatmap
#' barplot_stateloadings(Statescope, top_n = 1)
#'
barplot_stateloadings <- function(Statescope, top_n = 1) {
    if (!is(Statescope, "SummarizedExperiment")) { ## Check Statescope input
        stop("Statescope is not a SummarizedExperiment object")
    } else if (top_n > length(names(metadata(Statescope)$stateloadings))) {
        stop("n is bigger than number of genes")
    }
    ## init df
    plot_df <- setNames(
        data.frame(matrix(ncol = 4, nrow = 0)),
        c("celltype", "state", "gene", "score")
    )
    ## Gather top genes and scores
    for (ct in names(metadata(Statescope)$stateloadings)) {
        for (state in colnames(metadata(Statescope)$stateloadings[[ct]])) {
            indices <- head(sort(metadata(Statescope)$stateloadings[[ct]]
            [, state], index.return = TRUE, decreasing = TRUE)$ix, top_n)
            genes <- rownames(metadata(
                Statescope
            )$stateloadings[[ct]][indices, ])
            scores <- metadata(Statescope)$stateloadings[[ct]][indices, state]

            ## add to plot df
            for (i in seq(top_n)) {
                plot_df[nrow(plot_df) + 1, ] <-
                    data.frame(ct, state, genes[i], scores[i])
            }
        }
    }

    ## Make barplot
    ggplot(plot_df, aes(x = celltype, y = score, fill = state)) +
        geom_bar(position = "dodge", stat = "identity") +
        geom_text(aes(label = gene),
            position = position_dodge(width = 0.9),
            vjust = 1.3
        )
}


## Make Python and dplyr (.) functions global for check
utils::globalVariables(c(
    "Framework_Iterative", "Purify_AllGenes", "biggest_drop",
    "cNMF", "cNMF_Retrieval", "find_threshold", "select_k", ".",
    "RMSE", "celltype", "correlation", "gene", "score"
))
