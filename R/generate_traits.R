#' Generate uncorrelated traits
#'
#' @param n_species \[`integer(1)`\]\cr{}
#'                  the number of species
#' @param n_traits \[`integer(1)`\]\cr{}
#'                  the number of traits
#'
#' @return a matrix with `n_species` rows and `n_traits` column with uniform
#' uncorrelated traits between 0 and 1
#'
#' @examples
#' generate_traits(3, 2)
#'
#' @export
generate_traits = function(n_species, n_traits) {
    traits = matrix(runif(n_species * n_traits), ncol = n_traits)

    rownames(traits) = paste0("sp", seq(nrow(traits)))
    colnames(traits) = paste0("trait", seq(ncol(traits)))

    return(traits)
}
