other_fric = function(simul_df, given_traits, trait_comb_df_rows) {

    purrr::map_dfr(seq(nrow(trait_comb_df_rows)), function(row_number) {
        # Name of the trait combination
        trait_comb_name = trait_comb_df_rows[row_number,][["trait_comb"]]

        # Actual traits that compose this combination
        given_trait_comb = trait_comb_name %>%
            strsplit("_", fixed = TRUE) %>%
            .[[1]]

        # Does combination contains a single trait
        single_dim = length(given_trait_comb) == 1

        if (single_dim) {
            fric_pool = diff(range(given_traits[, given_trait_comb,
                                                drop = FALSE]))
        } else {
            fric_pool = geometry::convhulln(
                given_traits[, given_trait_comb, drop = FALSE], "FA")[["vol"]]
        }

        simul_df %>%
            tidyr::nest(sp_df = c(species, abund, rel_abund)) %>%
            mutate(
                trait_comb = trait_comb_name,
                FRic_other = purrr::map_dbl(
                    sp_df, function(df) {
                        if (single_dim) {
                            diff(range(given_traits[df[["species"]],
                                                    given_trait_comb,
                                                    drop = FALSE]))
                        } else {
                            geometry::convhulln(given_traits[df[["species"]],
                                                             given_trait_comb,
                                                             drop = FALSE],
                                                "FA")[["vol"]]
                        }
                    }),
                FRic_stand = FRic_other/fric_pool) %>%
            select(-sp_df)
    })
}
