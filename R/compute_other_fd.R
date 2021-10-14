#' Compute FD indices without abundances
#'
#' Wrapper around [`FD::dbFD()`]
#' @inheritParams compute_dissim_di
#' @param trait_comb \[`list(character(1+))`\]\cr{}
#'                   a list contaning combinations of trait combinations
#'
#' @export
compute_other_fd = function(site_sp_df, trait_comb, trait_df,
                            var_type = "env") {

    var_sym = rlang::sym(var_type)

    site_sp_df = site_sp_df %>%
        select(seed, !!var_sym, species, abund) %>%
        mutate(new_var = paste0(seed, "__", !!var_sym)) %>%
        select(-seed, -!!var_sym) %>%
        tidyr::spread(species, abund) %>%
        mutate_all(~ifelse(is.na(.), 0, .)) %>%
        as.data.frame() %>%
        tibble::column_to_rownames("new_var")

    names(trait_comb) = lapply(trait_comb, function(x) {
        paste(x, sep = "", collapse = "_")
    })

    all_fd = purrr::map_dfr(trait_comb, function(given_comb) {
        trait_df = trait_df[colnames(site_sp_df), given_comb, drop = FALSE]

        simul_fd = FD::dbFD(trait_df, site_sp_df, w.abun = TRUE,
                            calc.FRic = TRUE, stand.FRic = TRUE,
                            scale.RaoQ = TRUE, calc.FGR = FALSE,
                            calc.CWM = FALSE, messages = FALSE) %>%
            as.data.frame() %>%
            tibble::rownames_to_column("new_var") %>%
            tidyr::separate(new_var, c("seed", var_type), sep = "__",
                            convert = TRUE)
    }, .id = "trait_comb")
}

compute_flood_other_fd = function(site_sp_df, trait_comb, trait_df) {

    names(trait_comb) = lapply(trait_comb, function(x) {
        paste(x, sep = "", collapse = "_")
    })

    all_fd = purrr::map_dfr(trait_comb, function(given_comb) {

        # Keep only species in common between traits and site-species matrix
        common_species = intersect(rownames(trait_df), colnames(site_sp_df))

        trait_df_sub = trait_df[common_species, given_comb, drop = FALSE]

        site_sp_df_sub = site_sp_df[, common_species]

        zero_abund_site = rowSums(site_sp_df_sub) == 0

        site_sp_df_sub = site_sp_df_sub[!zero_abund_site,]

        # Compute functional diversity
        simul_fd = FD::dbFD(trait_df_sub, site_sp_df_sub, w.abun = FALSE,
                            calc.FRic = TRUE, stand.FRic = TRUE,
                            scale.RaoQ = TRUE, calc.FGR = FALSE,
                            calc.CWM = FALSE, messages = FALSE) %>%
            as.data.frame() %>%
            tibble::rownames_to_column("quadrat")
    }, .id = "trait_comb")
}
