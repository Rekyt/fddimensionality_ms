compute_global_di = function(global_dist) {
    global_com = matrix(1, ncol = nrow(global_dist),
                        dimnames = list(site = "glob",
                                        species = rownames(global_dist)))
    global_di = funrar::distinctiveness(global_com, global_dist)

    global_di_df = funrar::matrix_to_stack(global_di, "global_di")

    return(global_di_df)
}


compute_all_di = function(sim_df, global_dist = global_dist,
                          single_trait_dist   = single_trait_dist,
                          other_trait_dist    = other_trait_dist,
                          global_di_df        = global_di_df) {

    other_cols = setdiff(colnames(sim_df), c("abund", "rel_abund", "species"))
    join_cols = c("species", other_cols)
    # Local Di no abundance
    local_di_no_abund = compute_specific_di(sim_df, global_dist,
                                            "local_di_no_abund", NULL)
    local_di_abund    = compute_specific_di(sim_df, global_dist,
                                            "local_di_abund",
                                            "rel_abund")
    single_trait_di   = compute_local_global_di(sim_df, single_trait_dist,
                                                "single_trait_di", NULL)
    other_trait_di    = compute_local_global_di(sim_df, other_trait_dist,
                                                "other_trait_di", NULL)

    # Combine all Functional distinctiveness information
    local_di_abund %>%
        rename(Di = local_di_abund) %>%
        inner_join(global_di_df %>%
                       select(-site),
                   by = "species") %>%
        inner_join(local_di_no_abund, by = join_cols) %>%
        inner_join(single_trait_di, by = join_cols) %>%
        inner_join(other_trait_di, by = join_cols) %>%
        inner_join(sim_df %>%
                       rename(abundance = abund),
                   by = join_cols)
}


compute_specific_di = function(sim_df, dist_matrix, di_type = "global_di",
                               rel_abund_col = NULL) {

    di_name = rlang::sym(di_type)
    rel_abund_name = rel_abund_col

    sim_df %>%
        tidyr::nest(-seed, -sigma, -lim_sim) %>%
        mutate(data = purrr::map(data,
                                 ~funrar::distinctiveness_tidy(
                                     .x, "species", "env", abund = rel_abund_name,
                                     dist_matrix = dist_matrix) %>%
                                     rename(!!di_name := Di) %>%
                                     select(-abund, -rel_abund))) %>%
        tidyr::unnest(data)
}


compute_local_global_di = function(sim_df, dist_matrix, di_name = NULL,
                                   rel_abund_col = NULL) {

    if (is.null(di_name)) {
        stop("please provide a distinctiveness name")
    }

    global_name = paste0("global_", di_name)
    local_name  = paste0("local_", di_name)

    global_com = matrix(1, ncol = nrow(dist_matrix),
                        dimnames = list(site = "glob",
                                        species = rownames(dist_matrix)))

    global_di = funrar::distinctiveness(global_com, dist_matrix) %>%
        funrar::matrix_to_stack(global_name) %>%
        select(-site)

    local_di = compute_specific_di(sim_df, dist_matrix, local_name,
                                   rel_abund_col) %>%
        inner_join(global_di, by = "species")
}

#' Compute Functional Diversity Using All Dissimilarities
#'
#' @param trait_dissim_list \[`list(list(matrix(1+)))`\]\cr{}
#'                          a list containing list of dissimilarity matrices
#'                          with their names
#' @param site_sp_df        \[`data.frame()`\]\cr{}
#'                          a tidy species data.frame with a row defining a
#'                          species in a site with its abundance
#' @param trait_comb_df     \[`data.frame()`\]\cr{}
#'                          a data.frame with the list of trait combinations
#'                          considered as well as the number of traits included
#'
#' @export
compute_dissim_di = function(trait_dissim_list, site_sp_df, trait_comb_df) {
    trait_dissim_list %>%
        unlist(recursive = FALSE) %>%
        purrr::imap(function(dissim, dissim_name) {
            site_sp_df %>%
                compute_local_global_di(dissim %>%
                                            as.matrix(),
                                        di_name = dissim_name,
                                        rel_abund_col = "rel_abund")
        }) %>%
        {Reduce(function(x, y) inner_join(x, y), .)} %>%
        tidyr::gather("di_type", "di_value", starts_with("local"),
                      starts_with("global")) %>%
        tidyr::extract(di_type, c("di_scale", "trait_comb"),
                       "([a-z]+)_(.*)") %>%
        tidyr::spread(di_scale, di_value) %>%
        inner_join(trait_comb_df, by = "trait_comb") %>%
        mutate(
            contains_trait1 = ifelse(
                grepl("trait1", trait_comb, fixed = TRUE), TRUE, FALSE))
}
