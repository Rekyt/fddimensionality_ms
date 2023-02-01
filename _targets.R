# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed. # nolint

# Set target options:
tar_option_set(
  packages = c("tibble", "dplyr"), # packages that your targets need to run
  format = "rds" # default storage format
)

# tar_make_clustermq() configuration
options(clustermq.scheduler = "multiprocess")

# Load all functions in the R/ folder
tar_source()

list(
    # Simulation Parameters ----------------------------------------------------
    # Constants
    tar_target(n_seed, 10L),
    tar_target(n_traits, 11L),
    tar_target(n_species, 200L),
    tar_target(n_com, 10L),
    tar_target(n_indiv, 500L),
    tar_target(null_n, seq(80)),

    # Environmental Gradients
    tar_target(env_gradient, seq(0.1,   0.9, length.out = n_com)),
    tar_target(var_sigmas,   seq(0.2,  0.03, length.out = n_com)),
    tar_target(str_sigmas,   seq(0.1, 0.005, length.out = n_com)),

    # Migration rate
    tar_target(mig_rates, seq(1, 0.1, length.out = 4)),

    # Simulation Seeds
    tar_target(
        all_seeds,
        {
            set.seed(20190724)
            sample(1e6, size = n_seed)
        }
    ),
    tar_target(
        phylo_seeds,
        {
            set.seed(20190719)
            sample(1e6, size = 1)
        }
    ),

    # Parameter sets
    tar_target(
        single_trait_param,
        data.frame(env   = env_gradient,
                   sigma = var_sigmas) %>%
            tidyr::crossing(seed = all_seeds)
    ),
    tar_target(
        two_traits_param,
        data.frame(first_sigma = var_sigmas,
                   second_sigma = rev(var_sigmas),
                   env = env_gradient) %>%
            tidyr::crossing(seed = all_seeds)
    ),
    tar_target(
        mig_param,
        data.frame(env   = env_gradient,
                   sigma = var_sigmas) %>%
            tidyr::crossing(seed = all_seeds, mig = mig_rates)
    ),

    # Traits and Pools ---------------------------------------------------------
    # Traits
    tar_target(traits, generate_traits(n_species, n_traits)),
    tar_target(traits_sc, scale(traits)),
    tar_target(
        traits_df,
        tibble::rownames_to_column(data.frame(traits), "species")
    ),
    tar_target(traits_df_sc, mutate_if(traits_df, is.numeric, ~scale(.))),

    # Species Pools
    tar_target(smallest_pool, create_pool(traits_df, n_indiv/50)),
    tar_target(smaller_pool,  create_pool(traits_df, n_indiv/5)),
    tar_target(big_pool,      create_pool(traits_df, n_indiv)),

    # Compute trait dissimilarities
    tar_target(
        trait_dissim_all,
        compute_all_trait_combinations_dissimilarity(traits_sc)
    ),
    tar_target(
        trait_comb_df,
        get_trait_comb_length(trait_dissim_all)
    ),
    tar_target(
        trait_comb_df_sliced_grouped,
        tar_group(slice_smaller(trait_comb_df)),
        iteration = "group"
    ),

    # Null Traits
    tar_target(
        null_traits,
        generate_null_traits(traits, null_n)
    ),


    # Single Trait simulations -------------------------------------------------
    tar_target(
        var_sigma_simul,
        single_ecolottery_simulation(single_trait_param, big_pool, n_indiv)
    ),
    tar_target(
        var_sigma_df,
        format_single_simulation(var_sigma_simul, pool = "big_pool")
    ),
    tar_target(
        fd_var_sigma_obs,
        compute_fd(var_sigma_df, traits, trait_comb_df_sliced_grouped),
        pattern = map(trait_comb_df_sliced_grouped),
        iteration = "list"
    )
)


#     # Presence-Absence Null Models
    # var_sigma_null = var_sigma_df %>%
    #     select(seed, env, species, abund) %>%
    #     mutate(seed_env = paste0(seed, "__", env)) %>%
    #     select(-seed, -env) %>%
    #     tidyr::spread(species, abund) %>%
    #     mutate_all(~ifelse(is.na(.), 0, .)) %>%
    #     as.data.frame() %>%
    #     tibble::column_to_rownames("seed_env") %>%
    #     vegan::nullmodel("curveball") %>%
    #     simulate(nsim = 100) %>%
    #     narray::split(3),
#     # Trait Null Models
#     null_traits = target({
#
#         null_traits = traits
#         rownames(null_traits) = sample(rownames(null_traits),
#                                        nrow(null_traits),
#                                        replace = FALSE)
#         null_traits
#     },
#     transform = map(null_n = !!null_n)
#     ),
#     # Functional Diversity
#     var_sigma_other_fd_row = target(
#         compute_other_fd(var_sigma_df, trait_comb, traits),
#         transform = map(trait_comb)
#     ),
#     var_sigma_other_fd = target(
#         dplyr::bind_rows(var_sigma_other_fd_row),
#         transform = combine(var_sigma_other_fd_row)
#     ),
#     var_sigma_fric_row = target(
#         other_fric(var_sigma_df, traits, trait_comb_df),
#         transform = split(trait_comb_df, slices = 25)
#     ),
#     var_sigma_fric = target(
#         dplyr::bind_rows(var_sigma_fric_row),
#         transform = combine(var_sigma_fric_row)
#     ),
#     var_sigma_fdis_row = target(
#         other_fdis(var_sigma_df, traits, trait_comb_df),
#         transform = split(trait_comb_df, slices = 25)
#     ),
#     var_sigma_fdis = target(
#         dplyr::bind_rows(var_sigma_fdis_row),
#         transform = combine(var_sigma_fdis_row)
#     ),
#     # Functional Diversity on null models
#     var_sigma_null_trait_other_fd_row = target(
#         compute_other_fd(var_sigma_df, trait_comb, null_traits) %>%
#             mutate(null_trait = null_n),
#         transform = cross(trait_comb, null_traits)
#     ),
#     var_sigma_null_trait_fric_row = target(
#         other_fric(var_sigma_df, null_traits, trait_comb_df_row) %>%
#             mutate(null_trait = null_n),
#         transform = cross(trait_comb_df_row, null_traits)
#     ),
#     var_sigma_null_trait_fdis_row = target(
#         other_fdis(var_sigma_df, null_traits, trait_comb_df_row) %>%
#             mutate(null_trait = null_n),
#         transform = cross(trait_comb_df_row, null_traits)
#     ),
#     var_sigma_null_trait_other_fd = target(
#         dplyr::bind_rows(var_sigma_null_trait_other_fd_row),
#         transform = combine(var_sigma_null_trait_other_fd_row)
#     ),
#     var_sigma_null_trait_fric = target(
#         dplyr::bind_rows(var_sigma_null_trait_fric_row),
#         transform = combine(var_sigma_null_trait_fric_row)
#     ),
#     var_sigma_null_trait_fdis = target(
#         dplyr::bind_rows(var_sigma_null_trait_fdis_row),
#         transform = combine(var_sigma_null_trait_fdis_row)
#     ),
#     var_sigma_null_trait_fd_ses = var_sigma_null_trait_other_fd %>%
#         tidyr::gather("fd_index", "obs_value", nbsp:FDiv, -null_trait) %>%
#         filter(fd_index %in% c("FRic", "FDis", "RaoQ")) %>%
#         rename(null_value = obs_value) %>%
#         group_by(trait_comb, seed, env, fd_index) %>%
#         summarise(null_mean = mean(null_value),
#                   null_sd   = sd(null_value),
#                   null_ecdf = lst(ecdf(null_value))) %>%
#         inner_join(var_sigma_other_fd %>%
#                        tidyr::gather("fd_index", "obs_value", nbsp:FDiv),
#                    by = c("trait_comb", "seed", "env", "fd_index")) %>%
#         mutate(fd_ses = (obs_value - null_mean)/null_sd,
#                empirical_p_fd = purrr::map2_dbl(null_ecdf, obs_value,
#                                                 ~.x(.y))) %>%
#         ungroup(),
#     # Single-Trait Variable Migration Rate ---------------------------------
#     var_mig_simul = target(
#         ecolottery::coalesc(J = 100, m = mig, pool = big_pool,
#                             filt = function(x) {
#                                 given_seed = seed
#                                 filt_gaussian(x, env, sigma = sigma)
#                             }),
#         transform = map(.data = !!mig_param)
#     ),
#     var_mig_df_single = target(
#         format_simulation_big_df(var_mig_simul),
#         transform = map(var_mig_simul)
#     ),
#     var_mig_df = target(
#         bind_rows(var_mig_df_single),
#         transform = combine(var_mig_df_single)
#     ),
#     var_mig_other_fd_row = target(
#         var_mig_df %>%
#             mutate(m_env = paste(m, env, sep = "-_-")) %>%
#             compute_other_fd(trait_comb, traits, "m_env"),
#         transform = map(trait_comb)
#     ),
#     var_mig_other_fd = target(
#         dplyr::bind_rows(var_mig_other_fd_row) %>%
#             tidyr::separate(m_env, c("m", "env"), sep = "-_-",
#                             convert = TRUE),
#         transform = combine(var_mig_other_fd_row)
#     ),
#     var_mig_null_trait_other_fd_row = target(
#         var_mig_df %>%
#             mutate(m_env = paste(m, env, sep = "-_-")) %>%
#             compute_other_fd(trait_comb, null_traits, "m_env") %>%
#             mutate(null_trait = null_n),
#         transform = cross(trait_comb, null_traits)
#     ),
#     var_mig_null_trait_other_fd = target(
#         dplyr::bind_rows(var_mig_null_trait_other_fd_row) %>%
#             tidyr::separate(m_env, c("m", "env"), sep = "-_-",
#                             convert = TRUE),
#         transform = combine(var_mig_null_trait_other_fd_row)
#     ),
#     var_mig_all_fd = var_mig_other_fd %>%
#         tidyr::gather("fd_index", "obs_value", nbsp:FDiv) %>%
#         inner_join(var_mig_null_trait_other_fd %>%
#                        tidyr::gather("fd_index", "null_value",
#                                      nbsp:FDiv) %>%
#                        group_by(trait_comb, seed, m, env, fd_index) %>%
#                        summarise(null_mean = mean(null_value, na.rm = TRUE),
#                                  null_sd   = sd(null_value,
#                                                 na.rm = TRUE)) %>%
#                        ungroup()) %>%
#         mutate(fd_ses = (obs_value - null_mean)/(null_sd)),
#     var_mig_sub_fd = var_mig_all_fd %>%
#         filter(fd_index %in% c("FRic", "FDis")) %>%
#         inner_join(trait_comb_df, by = "trait_comb") %>%
#         select(-null_mean, -null_sd) %>%
#         tidyr::gather("fd_type", "fd_value", obs_value, fd_ses) %>%
#         mutate(fd_type = factor(fd_type, c("obs_value", "fd_ses")),
#                fd_index = factor(fd_index, c("FRic", "FDis")),
#                contains_trait1 = cw_trait1(trait_comb)),
#     var_mig_sub_fd_reg = var_mig_sub_fd %>%
#         filter(!is.infinite(fd_value), !is.null(fd_value)) %>%
#         group_by(fd_index, fd_type, m, trait_comb, contains_trait1) %>%
#         tidyr::nest() %>%
#         ungroup() %>%
#         mutate(lm_mod = purrr::map(data, ~lm(fd_value ~ env, data = .x)),
#                lm_tidy = purrr::map(lm_mod, broom::tidy),
#                lm_glan = purrr::map(lm_mod, broom::glance)) %>%
#         tidyr::unnest(lm_glan) %>%
#         filter(p.value < 0.05) %>%
#         tidyr::unnest(lm_tidy, names_repair = "universal") %>%
#         select_at(vars(-(std.error:p.value...12))) %>%
#         tidyr::spread(term, estimate),
#     # Single-Trait Variable Pool Size --------------------------------------
#     small_pool_simul = target(
#         ecolottery::coalesc(J = 100, m = 1, pool = smallest_pool,
#                             filt = function(x) {
#                                 given_seed = seed
#                                 filt_gaussian(x, env, sigma = sigma)
#                             }),
#         transform = map(.data = !!single_trait_param)
#     ),
#     # Assemble simulation in tidy data.frame
#     small_pool_df_single = target(
#         format_simulation_big_df(small_pool_simul),
#         transform = map(small_pool_simul)
#     ),
#     small_pool_df = target(
#         bind_rows(small_pool_df_single),
#         transform = combine(small_pool_df_single)
#     ),
#     small_pool_other_fd_row = target(
#         compute_other_fd(small_pool_df, trait_comb, traits),
#         transform = map(trait_comb)
#     ),
#     small_pool_other_fd = target(
#         dplyr::bind_rows(small_pool_other_fd_row),
#         transform = combine(small_pool_other_fd_row)
#     ),
#     small_pool_null_trait_other_fd_row = target(
#         compute_other_fd(small_pool_df, trait_comb, null_traits) %>%
#             mutate(null_trait = null_n),
#         transform = cross(trait_comb, null_traits)
#     ),
#     small_pool_null_trait_other_fd = target(
#         dplyr::bind_rows(small_pool_null_trait_other_fd_row),
#         transform = combine(small_pool_null_trait_other_fd_row)
#     ),
#     small_pool_null_trait_fd_ses = small_pool_null_trait_other_fd %>%
#         tidyr::gather("fd_index", "obs_value", nbsp:FDiv, -null_trait) %>%
#         filter(fd_index %in% c("FRic", "FDis", "RaoQ")) %>%
#         rename(null_value = obs_value) %>%
#         group_by(trait_comb, seed, env, fd_index) %>%
#         summarise(null_mean = mean(null_value),
#                   null_sd   = sd(null_value),
#                   null_ecdf = lst(ecdf(null_value))) %>%
#         inner_join(small_pool_other_fd %>%
#                        tidyr::gather("fd_index", "obs_value", nbsp:FDiv),
#                    by = c("trait_comb", "seed", "env", "fd_index")) %>%
#         mutate(fd_ses = (obs_value - null_mean)/null_sd,
#                empirical_p_fd = purrr::map2_dbl(null_ecdf, obs_value,
#                                                 ~.x(.y))) %>%
#         ungroup(),
#     # Power analysis on single trait ---------------------------------------
#     var_sigma_ses = var_sigma_null_trait_fd_ses %>%
#         select(-null_mean, -null_sd, -null_ecdf) %>%
#         inner_join(trait_comb_df) %>%
#         mutate(contains_trait1 = cw_trait1(trait_comb)) %>%
#         select(-empirical_p_fd) %>%
#         # Compare SES to a bunche of threshold
#         mutate(ses_comp = purrr::map(fd_ses, function(x) {
#             tibble(
#                 ses_threshold = seq(0, 2.75, length.out = 250),
#                 ses_comp = abs(x) >= ses_threshold)
#         })) %>%
#         tidyr::unnest(ses_comp) %>%
#         group_by(fd_index, trait_comb, trait_num, contains_trait1,
#                  ses_threshold) %>%
#         summarise(n_detected = sum(ses_comp)) %>%
#         ungroup(),
#     var_sigma_ses_prop = var_sigma_ses %>%
#         filter(contains_trait1) %>%
#         mutate(prop_detected = n_detected/30),
#     var_sigma_ses_power = var_sigma_ses_prop %>%
#         filter(contains_trait1) %>%
#         group_by(fd_index, ses_threshold, trait_num) %>%
#         summarise(ses_power = mean(prop_detected, na.rm = TRUE)),
#     var_sigma_ses_alpha = var_sigma_ses %>%
#         select(-trait_comb) %>%
#         filter(!contains_trait1) %>%
#         mutate(prop_detected = n_detected/30) %>%
#         group_by(fd_index, ses_threshold, trait_num) %>%
#         summarise(ses_alpha = mean(prop_detected, na.rm = TRUE)),
#     # Two Opposing Traits --------------------------------------------------
#     two_var_simul = target(
#         ecolottery::coalesc(J = 100, m = 1, pool = big_pool,
#                             filt = function(x) {
#                                 given_seed = seed
#                                 filt_gaussian(x, env, first_sigma) *
#                                     filt_gaussian(x[2], env, second_sigma)
#                             }),
#         transform = map(.data = !!two_traits_param)
#     ),
#     two_var_df = target(
#         format_simulation_big_df_two_traits(two_var_simul),
#         transform = combine(two_var_simul)
#     ),
#     # Perform null models both in presence-absence and abundance
#     two_var_null = two_var_df %>%
#         select(seed, env, species, abund) %>%
#         mutate(seed_env = paste0(seed, "__", env)) %>%
#         select(-seed, -env) %>%
#         tidyr::spread(species, abund) %>%
#         mutate_all(~ifelse(is.na(.), 0, .)) %>%
#         as.data.frame() %>%
#         tibble::column_to_rownames("seed_env") %>%
#         vegan::nullmodel("curveball") %>%
#         simulate(nsim = 100) %>%
#         narray::split(3),
#     # Compute Functional Diversity
#     two_var_other_fd_row = target(
#         compute_other_fd(two_var_df, trait_comb, traits),
#         transform = map(trait_comb)
#     ),
#     two_var_other_fd = target(
#         dplyr::bind_rows(two_var_other_fd_row),
#         transform = combine(two_var_other_fd_row)
#     ),
#     # Functional Diversity on Null Models
#     two_var_null_other_fd_row = target(
#         compute_null_var_sigma_other_fd(two_var_null, trait_comb,
#                                         traits, "new_var") %>%
#             tidyr::separate(new_var, c("seed", "env"), sep = "__",
#                             convert = TRUE),
#         transform = map(trait_comb)
#     ),
#     two_var_null_other_fd = target(
#         bind_rows(two_var_null_other_fd_row),
#         transform = combine(two_var_null_other_fd_row)
#     ),
#     # Functional Diversity on Trait Null Models
#     two_var_null_trait_other_fd_row = target(
#         compute_other_fd(two_var_df, trait_comb, null_traits) %>%
#             mutate(null_trait = null_n),
#         transform = cross(trait_comb, null_traits)
#     ),
#     two_var_null_trait_other_fd = target(
#         dplyr::bind_rows(two_var_null_trait_other_fd_row),
#         transform = combine(two_var_null_trait_other_fd_row)
#     ),
#     two_var_null_trait_fd_ses = two_var_null_trait_other_fd %>%
#         tidyr::gather("fd_index", "obs_value", nbsp:FDiv, -null_trait) %>%
#         filter(fd_index %in% c("FRic", "FDis", "RaoQ")) %>%
#         rename(null_value = obs_value) %>%
#         group_by(trait_comb, seed, env, fd_index) %>%
#         summarise(null_mean = mean(null_value),
#                   null_sd   = sd(null_value),
#                   null_ecdf = lst(ecdf(null_value))) %>%
#         inner_join(two_var_other_fd %>%
#                        tidyr::gather("fd_index", "obs_value", nbsp:FDiv),
#                    by = c("trait_comb", "seed", "env", "fd_index")) %>%
#         mutate(fd_ses = (obs_value - null_mean)/null_sd,
#                empirical_p_fd = purrr::map2_dbl(null_ecdf, obs_value,
#                                                 ~.x(.y)),
#                contains_trait = cw_trait1_trait2(trait_comb)) %>%
#         ungroup(),
#     # Combine Observed & Null Functional Diversity
#     two_var_other_fd_full = two_var_other_fd %>%
#         select(-qual.FRic, -FDiv) %>%
#         tidyr::gather("fd_index", "fd_obs", nbsp:RaoQ) %>%
#         inner_join(two_var_null_other_fd %>%
#                        select(-qual.FRic, -FDiv) %>%
#                        tidyr::gather("fd_index", "fd_null", nbsp:RaoQ) %>%
#                        group_by(trait_comb, seed, env, fd_index) %>%
#                        summarise(null_mean = mean(fd_null),
#                                  null_sd   = sd(fd_null))) %>%
#         mutate(fd_ses = (fd_obs - null_mean)/null_sd) %>%
#         inner_join(trait_comb_df) %>%
#         mutate(contains_trait = cw_trait1_trait2(trait_comb)),
#     # Figures from Simulations ---------------------------------------------
#     fig_theoretical = data.frame(
#         avg       = seq(0.1,  0.9, length.out = 5),
#         filter_sd = seq(0.2, 0.03, length.out = 5)) %>%
#         as_tibble() %>%
#         mutate(all_values = purrr::map2(
#             avg, filter_sd,
#             ~tibble(env_value = seq(0, 1, length.out = 1000),
#                     filt_value = dnorm(env_value, mean = .x, sd = .y),
#                     sc_abund = scale_zero_one(filt_value)))) %>%
#         tidyr::unnest(all_values) %>%
#         ggplot(aes(env_value, sc_abund, color = as.factor(avg), group = avg)) +
#         geom_line(size = 2) +
#         scale_color_viridis_d() +
#         labs(x     = "Environment",
#              y     = "Abundance (scaled)",
#              color = "Optimal\nTrait") +
#         theme_bw(12) +
#         theme(aspect.ratio = 1,
#               panel.grid = element_blank(),
#               axis.text = element_blank(),
#               axis.ticks = element_blank()),
#     fig_theoretical_two = cowplot::plot_grid(
#         fig_theoretical,
#         fig_theoretical +
#             scale_x_reverse() +
#             scale_color_viridis_d(option = "C",
#                                   direction = -1),
#         nrow = 2, ncol = 1
#     ),
#     fig_var_sigma_fd_null = var_sigma_null_trait_fd_ses %>%
#         select(-null_mean, -null_sd, -null_ecdf) %>%
#         inner_join(trait_comb_df) %>%
#         mutate(contains_trait1 = cw_trait1(trait_comb)) %>%
#         tidyr::gather("fd_type", "fd_value", obs_value, fd_ses,
#                       empirical_p_fd) %>%
#         filter(fd_index %in% c("FDis", "FRic"),
#                fd_type != "empirical_p_fd") %>%
#         plot_env_fd_obs_ses_single(),
#     fig_two_var_fd_indices = two_var_null_trait_fd_ses %>%
#         inner_join(trait_comb_df) %>%
#         mutate(contains_trait = cw_trait1_trait2(trait_comb)) %>%
#         tidyr::gather("fd_type", "fd_value", obs_value, fd_ses,
#                       empirical_p_fd) %>%
#         filter(fd_index %in% c("FDis", "FRic"),
#                fd_type != "empirical_p_fd") %>%
#         plot_env_fd_obs_ses_two(),
#     fig_var_mig_fd_null = var_mig_sub_fd %>%
#         filter(m %in% c(0.1, 1)) %>%
#         group_by(fd_index) %>%
#         group_map(function(grouped_data, group_name) {
#             index_name = as.character(pull(group_name))
#             grouped_data %>%
#                 ggplot(aes(env, fd_value, color = contains_trait1,
#                            group = trait_comb)) +
#                 geom_hline(data = data.frame(
#                     fd_type = factor("fd_ses", levels = c("obs_value",
#                                                           "fd_ses")),
#                     trait_num = 1:7,
#                     yint = 0),
#                     aes(yintercept = yint), linetype = 2, size = 1/2,
#                     alpha = 2/3) +
#                 # Axes over and under 0 for SES
#                 geom_hline(
#                     data = data.frame(
#                         fd_type = factor("fd_ses", levels = c("obs_value",
#                                                               "fd_ses")),
#                         trait_num = rep(1:7, 2),
#                         yint = rep(c(-1.96, 1.96), each = 7)),
#                     aes(yintercept = yint), linetype = 3, size = 1/2,
#                     alpha = 2/3) +
#                 geom_point(size = 1/2, alpha = 1/6) +
#                 # Regression Line
#                 geom_abline(
#                     data = var_mig_sub_fd_reg %>%
#                         inner_join(trait_comb_df) %>%
#                         semi_join(
#                             grouped_data,
#                             by = c("fd_index","fd_type", "m",
#                                    "trait_comb")),
#                     size = 0.8, alpha = 1/3,
#                     aes(slope = env, intercept = `(Intercept)`,
#                         color = contains_trait1)) +
#                 # Facet Specification
#                 facet_grid(vars(fd_type, m), vars(trait_num),
#                            labeller = labeller(
#                                fd_type = as_labeller(
#                                    c(fd_ses = "SES",
#                                      obs_value = "Observed")),
#                                trait_num = trait_num_lab,
#                                m = function(x) paste0("m = ", x)),
#                            scales = "free_y") +
#                 labs(x = "Environment", y = "Functional Diversity",
#                      color = "Contains Filtered Trait",
#                      title = index_name) +
#                 # Cosmetic changes
#                 guides(color = guide_legend(override.aes = list(alpha = 1,
#                                                                 size = 1))) +
#                 theme_bw() +
#                 theme(text = element_text(size = 12),
#                       axis.text.x = element_text(angle = 35, vjust = 0.7),
#                       panel.grid = element_blank(),
#                       strip.background = element_blank(),
#                       panel.spacing = unit(0.1, "cm"),
#                       aspect.ratio = 1,
#                       legend.position = "top")
#         }, keep = TRUE),
#     fig_smallest_pool_ses = small_pool_null_trait_fd_ses %>%
#         mutate(pool = "smallest_pool") %>%
#         bind_rows(var_sigma_null_trait_fd_ses %>%
#                       mutate(pool = "big_pool")) %>%
#         inner_join(trait_comb_df) %>%
#         filter(fd_index %in% c("FRic", "FDis")) %>%
#         mutate(cw1 = cw_trait1(trait_comb)) %>%
#         select(-starts_with("null"), -empirical_p_fd) %>%
#         tidyr::gather("fd_type", "fd_value", obs_value:fd_ses) %>%
#         ggplot(aes(env, fd_value, color = cw1, group = trait_comb)) +
#         geom_point(alpha = 1/6, size = 1/2) +
#         geom_hline(aes(yintercept = y), linetype = 2,
#                    data = data.frame(fd_type = rep("fd_ses", 2),
#                                      y = c(-1.96, 1.96))) +
#         stat_smooth(method = "lm", geom = "line", size = 1, alpha = 1/5) +
#         facet_grid(vars(fd_type, fd_index), vars(trait_num, pool),
#                    scales = "free_y") +
#         labs(x = "Environment",
#              y = "Functional Diversity",
#              color = "Contains Filtered Trait") +
#         theme(aspect.ratio = 1,
#               panel.grid = element_blank(),
#               strip.background = element_blank(),
#               legend.position = "top"),
#     fig_power_ses_threshold = var_sigma_ses_prop %>%
#         filter(fd_index != "RaoQ") %>%
#         ggplot(aes(ses_threshold, prop_detected,
#                    color = as.factor(trait_num))) +
#         geom_vline(xintercept = 1.96, linetype = 2, size = 0.9) +
#         geom_smooth(se = FALSE) +
#         facet_wrap(~fd_index) +
#         labs(x = "SES Threshold",
#              y = expression("Power"~(1 - beta)),
#              color = "Number of traits") +
#         scale_color_viridis_d() +
#         theme(aspect.ratio = 1,
#               legend.position = "top",
#               panel.grid = element_blank()),
#     fig_ses_power_alpha = var_sigma_ses_power %>%
#         full_join(var_sigma_ses_alpha) %>%
#         ggplot(aes(ses_alpha, ses_power, color = as.factor(trait_num))) +
#         geom_abline(slope = 1, intercept = 0, linetype = 2) +
#         geom_smooth(se = FALSE, size = 1) +
#         labs(x = expression(alpha),
#              y = expression("Power"~(1-beta)),
#              color = "Number of traits") +
#         scale_color_viridis_d() +
#         facet_wrap(vars(fd_index)) +
#         theme(aspect.ratio = 1,
#               legend.position = "top",
#               panel.grid = element_blank()),
#     fig_ses_distribution = var_sigma_null_trait_fd_ses %>%
#         mutate(contains_trait1 = cw_trait1(trait_comb)) %>%
#         inner_join(trait_comb_df) %>%
#         filter(fd_index != "RaoQ") %>%
#         group_by(fd_index) %>%
#         group_map(~.x %>%
#                       ggplot(aes(factor(trait_num, levels = 7:1), fd_ses,
#                                  color = contains_trait1)) +
#                       geom_violin(position = position_dodge()) +
#                       geom_point(alpha = 1/6,
#                                  position = position_jitterdodge()) +
#                       geom_hline(yintercept = 0, linetype = 2) +
#                       geom_hline(yintercept = c(-1.96, 1.96), linetype = 2,
#                                  size = 1) +
#                       coord_flip() +
#                       labs(y = "Standard Effect Size",
#                            x = "Number of Traits",
#                            color = "Contains Filtered Trait",
#                            title = .y) +
#                       theme(aspect.ratio = 1,
#                             legend.position = "top")) %>%
#         {filtered_trait_legend = cowplot::get_legend(.[[1]])
#         hu = purrr::map(., ~.x + theme(legend.position = "none"))
#         cowplot::plot_grid(
#             filtered_trait_legend,
#             cowplot::plot_grid(plotlist = hu, labels = NULL, ncol = 2),
#             nrow = 2, ncol = 1, rel_heights = c(0.05, 0.95))
#         },
#     fig_ses_distribution_pool_size = readd(small_pool_null_trait_fd_ses) %>%
#         mutate(pool = "smallest_pool") %>%
#         bind_rows(readd(var_sigma_null_trait_fd_ses) %>%
#                       mutate(pool = "big_pool")) %>%
#         inner_join(readd(trait_comb_df)) %>%
#         filter(fd_index %in% c("FRic", "FDis")) %>%
#         mutate(cw1 = cw_trait1(trait_comb)) %>%
#         select(-starts_with("null"), -empirical_p_fd) %>%
#         group_by(fd_index) %>%
#         group_map(
#             ~.x %>%
#                 ggplot(aes(factor(rev(trait_num)), fd_ses, color = cw1)) +
#                 geom_violin() +
#                 geom_point(size = 1/3, alpha = 1/6,
#                            position = position_jitterdodge()) +
#                 geom_hline(yintercept = 0, linetype = 2) +
#                 geom_hline(yintercept = c(-1.96, 1.96), linetype = 2,
#                            size = 1) +
#                 facet_grid(
#                     rows = vars(pool),
#                     labeller = labeller(
#                         pool = c(big_pool = "Big Species Pool",
#                                  smallest_pool = "Small Species Pool"))) +
#                 labs(x = "Trait Number",
#                      y = "Functional Diversity SES",
#                      color = "Contains Filtered Trait",
#                      title = .y) +
#                 coord_flip() +
#                 theme(aspect.ratio = 1,
#                       legend.position = "top")) %>%
#         cowplot::plot_grid(plotlist = ., ncol = 2)
# )
