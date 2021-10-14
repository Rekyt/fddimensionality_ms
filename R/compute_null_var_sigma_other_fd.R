compute_null_var_sigma_other_fd = function(null_com_list, trait_comb, trait_df,
                                           var_type) {
    purrr::map_dfr(null_com_list, function(null_com) {

        names(trait_comb) = lapply(trait_comb, function(x) {
            paste(x, sep = "", collapse = "_")
        })

        all_fd = purrr::map_dfr(trait_comb, function(given_comb) {
            trait_df = trait_df[colnames(null_com), given_comb, drop = FALSE]

            simul_fd = FD::dbFD(trait_df, null_com, w.abun = FALSE,
                                calc.FRic = TRUE, stand.FRic = TRUE,
                                scale.RaoQ = TRUE, calc.FGR = FALSE,
                                calc.CWM = FALSE, messages = FALSE) %>%
                as.data.frame() %>%
                tibble::rownames_to_column(var_type)
        }, .id = "trait_comb")
    }, .id = "sim")
}
