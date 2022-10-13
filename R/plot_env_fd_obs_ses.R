#' Figure of FD vs. Environment Comparing Observed and SES values
#'
#' @param full_fd_df \[`data.frame()`\]\cr{}
#'                   data.frame containing both observed and SES values of
#'                   functional diversity
#'
#' @export
plot_env_fd_obs_ses_two = function(full_fd_df) {
    full_fd_split = full_fd_df %>%
        filter(fd_index %in% c("FRic", "FDis", "RaoQ")) %>%
        mutate(fd_type = factor(fd_type, c("obs_value", "fd_ses",
                                           "empirical_p_fd")),
               fd_index = factor(fd_index, c("FRic", "FDis", "RaoQ")),
               contains_trait = contains_trait %>%
                   as.factor() %>%
                   forcats::fct_relevel("trait1 & trait2", after = 3L)) %>%
        group_split(fd_index)

    max_number_trait = max(full_fd_split[[1]]$trait_num)

    fig_split = lapply(full_fd_split, function(x) {

        fd_ind = unique(x$fd_index)

        ggplot(x, aes(env, fd_value, color = contains_trait,
                      group = trait_comb)) +
            # x-axis for SES
            geom_hline(
                data = data.frame(
                    fd_type = factor("fd_ses", levels = c("fd_obs", "fd_ses")),
                    trait_num = 1:max_number_trait,
                    yint = 0),
                aes(yintercept = yint), linetype = 2, size = 1/2,
                alpha = 2/3) +
            # Axes over and under 0 for SES
            geom_hline(
                data = data.frame(
                    fd_type = factor("fd_ses", levels = c("fd_obs", "fd_ses")),
                    trait_num = rep(1:max_number_trait, 2),
                    yint = rep(c(-1.96, 1.96), max_number_trait, each = TRUE)),
                aes(yintercept = yint), linetype = 3, size = 1/2,
                alpha = 2/3) +
            # Other geoms
            geom_point(alpha = 1/5, size = 1/2) +
            stat_smooth(geom = "line", size = 2/3, alpha = 1/3) +
            # Facetting + labels
            facet_grid(vars(fd_type), vars(trait_num), scales = "free_y",
                       labeller = labeller(
                           fd_type = c(fd_ses = "SES",
                                       obs_value = "Observed",
                                       empirical_p_fd = "p-value"),
                           trait_num = trait_num_lab
                       )) +
            labs(x = "Environment",
                 y = "Functional Diversity",
                 color = "One or Two Filtered traits?",
                 title = fd_ind) +
            guides(color = guide_legend(override.aes = list(alpha = 1,
                                                            size = 1))) +
            theme_bw() +
            theme(text = element_text(size = 12),
                  axis.text.x = element_text(angle = 35, vjust = 0.7),
                  panel.grid = element_blank(),
                  strip.background = element_blank(),
                  panel.spacing = unit(0.1, "cm"),
                  aspect.ratio = 1,
                  legend.position = "top")
    })
}
plot_env_fd_obs_ses_single = function(full_fd_df) {

    max_number_trait = max(full_fd_df$trait_num)

    group_plots = full_fd_df %>%
        mutate(fd_type = factor(fd_type, c("obs_value", "fd_ses",
                                           "empirical_p_fd")),
               fd_index = factor(fd_index, c("FRic", "FDis", "RaoQ"))) %>%
        group_by(fd_index) %>%
        group_map(function(grouped_data, group_name) {

            # Name of the Functional Diversity Index
            index_name = as.character(pull(group_name))

            # Subselect significant regression lines
            signif_reg = grouped_data %>%
                tidyr::nest(all_fd = c(seed, env, fd_value)) %>%
                mutate(lm_reg = purrr::map(all_fd, ~lm(fd_value ~ env,
                                                       data = .x)),
                       lm_sum = purrr::map(lm_reg, broom::glance)) %>%
                tidyr::unnest(lm_sum) %>%
                # Select signifiant regression
                # P-value (0.05) corrected for 127 combinations times 2 indices
                # with Bonferroni correction
                filter(p.value < 2e-4)

            # Actual Plot
            grouped_data %>%
                ggplot(aes(env, fd_value, color = contains_trait1,
                           group = trait_comb)) +
                # Central trend for SES
                geom_hline(data = data.frame(
                    fd_type = factor("fd_ses", levels = c("obs_value",
                                                          "fd_ses")),
                    trait_num = 1:max_number_trait,
                    yint = 0), aes(yintercept = yint), linetype = 3, size = 1/2,
                    alpha = 2/3) +
                # Axes over and under 0 for SES
                geom_hline(
                    data = data.frame(
                        fd_type = factor("fd_ses",
                                         levels = c("obs_value", "fd_ses")),
                        trait_num = rep(1:max_number_trait, 2),
                        yint = rep(c(1.96, -1.96), max_number_trait,
                                   each = TRUE)),
                    aes(yintercept = yint), linetype = 2, size = 1/2,
                    alpha = 2/3) +
                # Actual Data
                geom_point(alpha = 1/5, size = 1/2) +
                # Regression lines
                stat_smooth(data = grouped_data %>%
                                semi_join(signif_reg),
                            method = "lm", se = FALSE, geom = "line",
                            size = 2/3, alpha = 1/3) +
                # Facet Spec
                facet_grid(vars(fd_type), vars(trait_num),
                           labeller = labeller(
                               fd_type = as_labeller(
                                   c(fd_ses = "SES",
                                     obs_value = "Observed",
                                     empirical_p_fd = "p-value")),
                               trait_num = trait_num_lab),
                           scales = "free_y") +
                labs(x = "Environment", y = "Functional Diversity",
                     color = "Contains Filtered Trait",
                     title = index_name) +
                # Cosmetic changes
                guides(color = guide_legend(override.aes = list(alpha = 1,
                                                                size = 1))) +
                theme_bw() +
                theme(text = element_text(size = 12),
                      axis.text.x = element_text(angle = 35, vjust = 0.5),
                      panel.grid = element_blank(),
                      strip.background = element_blank(),
                      panel.spacing = unit(0.1, "cm"),
                      aspect.ratio = 1,
                      legend.position = "top")
        }, keep = TRUE)
}

