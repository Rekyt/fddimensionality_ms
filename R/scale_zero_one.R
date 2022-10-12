#' Scale Vector between Zero and One
#'
#' With a numeric vector scale between zero and one.
#'
#' @param vec a numeric vector
#'
#' @examples
#'
#' vec = stats::rnorm(1000)
#'
#' zero_one = scale_zero_one(vec)
#'
#' range(vec)
#' range(zero_one)
#'
#' @importFrom stats rnorm
#' @export
scale_zero_one = function(vec) {
    (vec - min(vec))/diff(range(vec))
}
