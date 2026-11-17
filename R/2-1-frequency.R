source("R/1-2-runtime.R")
df_in <- read_table(input)
HH_range <- HH_range[HH_range <= max(nchar(df_in$hap_clean))]
if (!length(HH_range)) stop("No applicable levels")
for (HH in HH_range) {
    df <- df_in %>% mutate(hap = hap_substr(hap_clean, HH))
    freq_df <- freq_table(df, hap_col = "hap", by_col = "population")
    freq_mat_df <- make_freq_matrix(freq_df, by_col = "population", hap_col = "hap", val_col = "freq")
    count_df <- df %>% count(population, hap, name = "n")
    count_mat_df <- make_freq_matrix(count_df, by_col = "population", hap_col = "hap", val_col = "n")
    write_result(freq_mat_df, file.path(output, sprintf("freq_matrix_HH%02d.tsv", HH)), row.names = FALSE, fileEncoding = "UTF-8")
    write_result(count_mat_df, file.path(output, sprintf("count_matrix_HH%02d.tsv", HH)), row.names = FALSE, fileEncoding = "UTF-8")
    hap_spread_tab <- compute_hap_spread_stats(count_df, freq_df)
    write_result(hap_spread_tab, file.path(output, sprintf("hap_subclade_diversity_HH%02d.tsv", HH)), row.names = FALSE, 
        fileEncoding = "UTF-8")
    bip_edges <- count_df %>% rename(count = n) %>% left_join(freq_df %>% select(population, hap, freq), by = c("population", 
        "hap")) %>% arrange(desc(freq))
    write_result(bip_edges, file.path(output, sprintf("bipartite_edges_HH%02d.tsv", HH)), row.names = FALSE, fileEncoding = "UTF-8")
    div_list <- lapply(split(freq_df, freq_df$population), function(x) {
        p <- x$freq
        data.frame(population = x$population[1], k = nrow(x), shannon = shannon(p), simpson = simpson(p), heff_expH = heff_shannon(p), 
            heff_1sumP2 = heff_simpson(p))
    })
    div_tab <- do.call(rbind, div_list)
    write_result(div_tab, file.path(output, sprintf("diversity_indices_HH%02d.tsv", HH)), row.names = FALSE, fileEncoding = "UTF-8")
    classify_prevalence <- function(freq_df, count_df) {
        x <- left_join(freq_df, count_df, by = c("population", "hap"))
        x$class <- with(x, ifelse(n.y == 1, "Singleton", ifelse(n.y == 2, "Doubleton", ifelse(freq <= thr_rare, 
            "Very rare", ifelse(freq <= thr_common, "Rare", "Common")))))
        detail <- x %>% select(population, hap, count = n.y, freq, class)
        summary <- detail %>% count(population, class, name = "hap_n") %>% pivot_wider(names_from = class, values_from = hap_n, 
            values_fill = 0)
        list(detail = detail, summary = summary)
    }
    pv <- classify_prevalence(freq_df, count_df)
    write_result(pv$detail, file.path(output, sprintf("hap_prevalence_classes_HH%02d.tsv", HH)), row.names = FALSE, 
        fileEncoding = "UTF-8")
    write_result(pv$summary, file.path(output, sprintf("hap_prevalence_summary_by_population_HH%02d.tsv", HH)), 
        row.names = FALSE, fileEncoding = "UTF-8")
    enr <- hap_enrichment(df, hap_col = "hap", pop_col = "population")
    write_result(enr, file.path(output, sprintf("enrichment_Fisher_HH%02d.tsv", HH)), row.names = FALSE, fileEncoding = "UTF-8")
}
