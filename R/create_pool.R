create_pool = function(traits_df, number_of_individuals) {

    traits_df %>%
        slice(rep(row_number(), number_of_individuals)) %>%
        mutate(indiv = row_number()) %>%
        select(indiv, species, starts_with("trait"))

}

get_trait_comb_length = function(trait_dissim_all) {

    trait_comb = names(trait_dissim_all)
    trait_num = unlist(lapply(strsplit(trait_comb, "_", fixed = TRUE), length))

    data.frame(
        trait_comb = trait_comb,
        trait_num  = trait_num
    )

}


slice_smaller = function(trait_comb_df, cutoff = 100) {

    set.seed(20230201)

    split_groups = trait_comb_df %>%
        group_by(trait_num) %>%
        group_split()

    lapply(split_groups, function(x) {

        if (nrow(x) <= cutoff) {

            x

        } else {

            slice_sample(x, n = cutoff, replace = FALSE)
        }

    }) %>%
        bind_rows() %>%
        group_by(trait_num)
}
