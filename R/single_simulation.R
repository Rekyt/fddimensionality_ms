single_ecolottery_simulation = function(param_df, pool, n_indiv) {

    param_df %>%
        rowwise() %>%
        mutate(
            simulation = list(
                ecolottery::coalesc(
                    J = n_indiv, m = 1, pool = pool,
                    filt = function(x) {
                        given_seed = seed
                        filt_gaussian(x, env, sigma = sigma)
                    }
                )
            )
        ) %>%
        ungroup()

}

filt_gaussian = function(x, topt, sigma) {
    exp(-(x[1] - topt)^2/(2*sigma^2))
}

format_single_simulation = function(
        var_sigma_simul, lim_sim_coef = FALSE, pool_type = ""
) {

    # Extract simulation parameters
    var_sigma_simul %>%
        rowwise() %>%
        mutate(
            update_simul = list(
                simulation$com %>%
                    rename(species = sp) %>%
                    mutate(
                        seed = seed, env = env, sigma = sigma,
                        pool = pool_type, m = simulation$call$m,
                        lim_sim = ifelse(
                            !lim_sim_coef, 0, simulation$call$coeff.lim.sim
                        )
                    )
            )
        ) %>%
        select(update_simul) %>%
        tidyr::unnest(update_simul) %>%
        group_by(seed, pool, m, sigma, lim_sim, env, species) %>%
        summarise(abund = n()) %>%
        mutate(rel_abund = abund/sum(abund)) %>%
        ungroup()

}


generate_null_traits = function(traits, n_null) {

    null_traits = traits

    lapply(seq(n_null), function(x) {
        rownames(null_traits) = sample(
            rownames(null_traits), nrow(null_traits), replace = FALSE
        )

        null_traits
    }) %>%
        setNames(paste0("null_trait_", seq(n_null)))

}


compute_fd = function(site_sp_df, traits, trait_comb_df) {

    site_sp_df = site_sp_df %>%
        select(seed, env, species, abund) %>%
        mutate(new_var = paste0(seed, "__", env)) %>%
        select(-seed, -env) %>%
        tidyr::spread(species, abund) %>%
        mutate(across(where(is.numeric), ~ifelse(is.na(.x), 0, .x))) %>%
        as.data.frame() %>%
        tibble::column_to_rownames("new_var") %>%
        as.matrix()

    lapply(seq(nrow(trait_comb_df)), function(n_comb) {

        given_comb = trait_comb_df$trait_comb[[n_comb]] %>%
            strsplit("_", fixed = TRUE) %>%
            .[[1]]

        fric = fundiversity::fd_fric(
            traits[, given_comb, drop = FALSE], site_sp_df
        )

        fdis = fundiversity::fd_fdis(
            traits[, given_comb, drop = FALSE], site_sp_df
        )

        merge(fric, fdis, by = "site")

    }) %>%
        setNames(nm = trait_comb_df$trait_comb) %>%
        bind_rows(.id = "trait_comb") %>%
        tibble::rownames_to_column("new_var") %>%
        tidyr::separate(site, c("seed", "env"), sep = "__", convert = TRUE)

}
