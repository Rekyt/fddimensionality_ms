Error message obtained lass moment

```
fail fig_two_var_fd
Error: target fig_two_var_fd failed.
diagnose(fig_two_var_fd)$error$message:
  Problem while computing `..1 = fd_index == "Q"`.
diagnose(fig_two_var_fd)$error$calls:
  two_var_df %>% filter(fd_index == "Q") %>% mutate(contains_trait = contains_trait %>% 
    as.factor() %>% forcats::fct_relevel("trait1 & trait2", after = 3L)) %>% 
    ggplot(aes(as.numeric(sigma), fd_value, color = contains_trait, 
        group = trait_comb))
  ggplot(., aes(as.numeric(sigma), fd_value, color = contains_trait, 
    group = trait_comb))
  mutate(., contains_trait = contains_trait %>% as.factor() %>% 
    forcats::fct_relevel("trait1 & trait2", after = 3L))
  filter(., fd_index == "Q")
  filter.data.frame(., fd_index == "Q")
  filter_rows(.data, ..., caller_env = caller_env())
  filter_eval(dots, mask = mask, error_call = error_call)
  withCallingHandlers({
    mask$eval_all_filter(dots, env_filter)
}, error = function(e) {
    local_error_context(dots = dots, .index = env_filter$current_expression, 
```

