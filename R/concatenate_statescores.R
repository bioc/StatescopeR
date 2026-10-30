#' Concatenate statescore lists into one data frame
#'
#' \code{concatenate_statescores.R} Concatenate statescore lists into one
#' data frame
#'
#' @param statescores Named list of cell-type-specific statescore DataFrames,
#' typically obtained from S4Vectors::metadata(Statescope)$statescores.
#' Each DataFrame should contain samples as row names and score
#' components as columns.
#'
#' @return A data frame containing samples as rows and concatenated
#' cell-type-specific scores as columns.
#'
#' @importFrom S4Vectors metadata
#' @export
#'
#' @examples
#' ## Load StateDiscovered Statescope object
#' load(system.file("extdata", "example_Statescope_StateDiscovered.RData",
#'     package = "StatescopeR"
#' ))
#' # Extract state scores from a StateScope object
#' statescores <- S4Vectors::metadata(Statescope)$statescores
#'
#' # Concatenate scores across cell types
#' concatenated_scores <- concatenate_statescores(statescores)
#'
concatenate_statescores <- function(statescores) {
    # Convert each cell type DataFrame to a regular data.frame
    df_list <- lapply(names(statescores), function(celltype) {
        df <- as.data.frame(statescores[[celltype]])

        # Rename columns using celltype_colidx structure
        colnames(df) <- paste0(gsub(" ", "_", celltype), "_",seq_len(ncol(df)))
        df})

    # Combine all dataframes by sample row names
    result <- do.call(cbind, df_list)

    return(result)
}
