#' Compute Trait Dissimilarities
#'
#' For each provided trait combination in a list, returns a list of
#' dissimilarity matrices computed using `funrar`. This function computes
#' euclidean dissimilarity matrices
#' @param trait_comb_list \[`list(character(1+))`\]\cr{}
#'                        a list of character vectors containing the trait
#'                        combination on which dissimilarity should be computed
#' @param trait_df        \[`data.frame()`\]\cr{}
#'                        a trait data.frame with species as rows and traits as
#'                        columns
#'
#' @return a list of dissimilarity matrices with concatenated names from
#' trait combinations
#'
#' @export
compute_trait_dissim = function(trait_comb_list, trait_df) {
    lapply(trait_comb_list, function(single_trait_comb) {

        if (!all(single_trait_comb %in% colnames(trait_df))) {
            stop("A trait combination was absent from trait data.frame")
        }

        funrar::compute_dist_matrix(
            trait_df[, single_trait_comb, drop = FALSE],
            metric = "euclidean", center = TRUE, scale = TRUE)
    }) %>%
        purrr::set_names(nm = lapply(trait_comb_list, function(x) {
            paste(x, sep = "", collapse = "_")
        }))
}
