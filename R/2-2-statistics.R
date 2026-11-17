freq_table <- function(df, hap_col = "hap", by_col = "population") {
    df %>% count(.data[[by_col]], .data[[hap_col]], name = "n") %>% group_by(.data[[by_col]]) %>% mutate(freq = n/sum(n)) %>% 
        ungroup()
}
make_freq_matrix <- function(freq_df, by_col = "population", hap_col = "hap", val_col = "freq") {
    freq_df %>% select(all_of(c(by_col, hap_col, val_col))) %>% tidyr::pivot_wider(names_from = all_of(hap_col), 
        values_from = all_of(val_col), values_fill = 0) %>% as.data.frame()
}
shannon <- function(p) {
    p <- p[p > 0]
    -sum(p * log(p))
}
simpson <- function(p) {
    1 - sum(p^2)
}
heff_shannon <- function(p) {
    exp(shannon(p))
}
heff_simpson <- function(p) {
    1/sum(p^2)
}
hap_enrichment <- function(df, hap_col = "hap", pop_col = "population") {
    all_pops <- unique(df[[pop_col]])
    all_haps <- unique(df[[hap_col]])
    res <- list()
    for (h in all_haps) {
        df$flag <- df[[hap_col]] == h
        for (g in all_pops) {
            a <- sum(df$flag & df[[pop_col]] == g)
            b <- sum(!df$flag & df[[pop_col]] == g)
            c <- sum(df$flag & df[[pop_col]] != g)
            d <- sum(!df$flag & df[[pop_col]] != g)
            mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
            if (all(mat >= 0)) {
                ft <- suppressWarnings(fisher.test(mat))
                res[[length(res) + 1]] <- data.frame(hap = h, population = g, a = a, b = b, c = d[1], d = d[2], 
                  p = ft$p.value, odds = unname(ft$estimate))
            }
        }
    }
    out <- do.call(rbind, res)
    out$padj <- p.adjust(out$p, method = "BH")
    out[order(out$padj, out$p), ]
}
compute_hap_spread_stats <- function(count_df, freq_df) {
    base_ct <- count_df %>% group_by(hap) %>% summarise(total_n = sum(n), n_pop = n(), .groups = "drop")
    spread <- count_df %>% group_by(hap) %>% mutate(p = n/sum(n)) %>% summarise(shannon_pop = shannon(p), evenness_pop = ifelse(n() > 
        1, shannon_pop/log(n()), 1), .groups = "drop")
    pop_detail <- count_df %>% group_by(hap) %>% arrange(desc(n), population, .by_group = TRUE) %>% summarise(pop_detail = paste0(population, 
        ":", n, collapse = ";"), top_population = population[1], top_count = n[1], .groups = "drop")
    base_ct %>% left_join(spread, by = "hap") %>% left_join(pop_detail, by = "hap") %>% arrange(desc(total_n))
}
