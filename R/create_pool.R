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
