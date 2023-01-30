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
