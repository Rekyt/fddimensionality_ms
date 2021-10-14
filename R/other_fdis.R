other_fdis = function(simul_df, given_traits, trait_comb_df_rows) {

    purrr::map_dfr(seq(nrow(trait_comb_df_rows)), function(row_number) {
        # Name of the trait combination
        trait_comb_name = trait_comb_df_rows[row_number,][["trait_comb"]]

        # Actual traits that compose this combination
        given_trait_comb = trait_comb_name %>%
            strsplit("_", fixed = TRUE) %>%
            .[[1]]

        simul_df = simul_df %>%
            mutate(single_site = paste(seed, pool, m, sigma, lim_sim, env,
                                       sep = "-__-")) %>%
            select(single_site, species, abund) %>%
            tidyr::spread(species, abund) %>%
            as.data.frame() %>%
            tibble::column_to_rownames("single_site")

        trait_subset = given_traits[colnames(simul_df), given_trait_comb,
                                    drop = FALSE]
        FD::fdisp(dist(trait_subset), as.matrix(simul_df)) %>%
            .[["FDis"]] %>%
            tibble::enframe("single_site", "FDis") %>%
            tidyr::separate(single_site,
                            c("seed", "pool", "m", "sigma", "lim_sim", "env"),
                            sep = "-__-") %>%
            mutate(trait_comb = trait_comb_name) %>%
            select(everything(), trait_comb, FDis)
    })
}
