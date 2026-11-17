get_pop_order_by_ORpeak <- function(df_boot) {
    df_boot <- df_boot %>% mutate(Haplogroup_up = toupper(Haplogroup), is_OR = grepl("^(O|R)", Haplogroup_up))
    pop_tot <- df_boot %>% count(population, name = "pop_total")
    or_cnt <- df_boot %>% filter(is_OR) %>% count(population, Haplogroup, name = "or_count") %>% left_join(pop_tot, 
        by = "population") %>% mutate(or_freq = or_count/pop_total)
    best_or <- or_cnt %>% group_by(population) %>% summarise(max_or_freq = if (length(or_freq)) 
        max(or_freq, na.rm = TRUE)
    else -1, best_or_hap = Haplogroup[which.max(or_freq)], .groups = "drop")
    pop_info <- pop_tot %>% left_join(best_or, by = "population") %>% mutate(max_or_freq = ifelse(is.na(max_or_freq), 
        -1, max_or_freq)) %>% arrange(desc(max_or_freq), desc(pop_total), population)
    pop_info$population
}
bootstrap_once_with_div <- function(rep_id, df, n_total, n_per_sample) {
    idx <- sample.int(n_total, size = n_per_sample, replace = TRUE)
    sample_df <- df[idx, c("population", "Haplogroup"), drop = FALSE]
    counts <- sample_df %>% count(population, Haplogroup, name = "count")
    pop_counts <- sample_df %>% count(population, name = "pop_total")
    freq_tbl <- counts %>% left_join(pop_counts, by = "population") %>% mutate(freq = count/pop_total, rep = rep_id) %>% 
        select(rep, population, Haplogroup, freq)
    div_tbl <- freq_tbl %>% group_by(population, rep) %>% summarise(shannon = {
        p <- freq
        p <- p[p > 0]
        -sum(p * log(p))
    }, simpson = {
        p <- freq
        1 - sum(p^2)
    }, .groups = "drop")
    list(freq_tbl = freq_tbl, div_tbl = div_tbl)
}
get_sig_table_for_hap <- function(subdat, hap_name) {
    pop_counts <- subdat %>% group_by(population) %>% summarise(n_points = n(), .groups = "drop")
    pops_ok <- pop_counts %>% filter(n_points >= 2) %>% pull(population) %>% as.character()
    if (length(unique(pops_ok)) < 2) {
        return(tibble(Haplogroup = character(0), group1 = character(0), group2 = character(0), p = numeric(0), 
            p.adj = numeric(0), p.adj.signif = character(0)))
    }
    sub_for_test <- subdat %>% filter(population %in% pops_ok)
    stat_test <- sub_for_test %>% rstatix::pairwise_wilcox_test(freq ~ population, p.adjust.method = "BH")
    if (nrow(stat_test) == 0) {
        return(tibble(Haplogroup = character(0), group1 = character(0), group2 = character(0), p = numeric(0), 
            p.adj = numeric(0), p.adj.signif = character(0)))
    }
    stat_test %>% mutate(p.adj.signif = ifelse(is.na(p.adj.signif), "ns", p.adj.signif), Haplogroup = hap_name) %>% 
        select(Haplogroup, group1, group2, p, p.adj, p.adj.signif)
}
plot_hap_single <- function(subdat, pop_levels, pop_colors, n_per_sample, n_reps) {
    subdat <- subdat %>% group_by(population) %>% filter(any(freq > 0)) %>% ungroup()
    if (nrow(subdat) == 0) 
        return(ggplot() + ggtitle("no data for this hap in any population"))
    used_pops <- intersect(pop_levels, unique(subdat$population))
    subdat$population <- factor(subdat$population, levels = used_pops, ordered = TRUE)
    y_lower <- min(subdat$freq, na.rm = TRUE)
    y_upper <- max(subdat$freq, na.rm = TRUE)
    base_plot <- ggplot(subdat, aes(x = population, y = freq)) + ggbeeswarm::geom_quasirandom(aes(color = population), 
        method = qr_method, width = qr_width, alpha = qr_alpha, size = qr_size, show.legend = FALSE) + geom_boxplot(aes(color = population), 
        width = 0.5, outlier.shape = NA, fill = NA, size = 0.6, show.legend = FALSE) + stat_summary(fun = median, 
        geom = "crossbar", width = 0.5, fatten = 0, color = "red", size = 0.6) + scale_color_manual(values = pop_colors[used_pops]) + 
        labs(x = "Population", y = "Frequency", title = paste0(unique(subdat$Haplogroup), " (n=", n_per_sample, 
            " /rep, ", n_reps, " reps)")) + theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1, 
        vjust = 1, size = 8), axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9), 
        plot.title = element_text(hjust = 0.5, size = 9)) + coord_cartesian(ylim = c(y_lower, y_upper))
    if (!isTRUE(freqplot_show_sig)) 
        return(base_plot)
    stat_tbl <- tryCatch(subdat %>% rstatix::pairwise_wilcox_test(freq ~ population, p.adjust.method = "BH"), error = function(e) NULL)
    if (is.null(stat_tbl) || nrow(stat_tbl) == 0) 
        return(base_plot)
    if (!is.null(freqplot_sig_ref_pop) && freqplot_sig_ref_pop %in% levels(subdat$population)) {
        stat_tbl <- stat_tbl %>% filter(group1 == freqplot_sig_ref_pop | group2 == freqplot_sig_ref_pop)
    }
    stat_tbl <- stat_tbl %>% filter(!is.na(p.adj) & p.adj <= freqplot_sig_alpha) %>% arrange(p.adj)
    if (nrow(stat_tbl) == 0) 
        return(base_plot)
    if (!is.null(freqplot_sig_max_pairs) && nrow(stat_tbl) > freqplot_sig_max_pairs) {
        stat_tbl <- stat_tbl %>% slice_head(n = freqplot_sig_max_pairs)
    }
    comps <- lapply(seq_len(nrow(stat_tbl)), function(i) c(as.character(stat_tbl$group1[i]), as.character(stat_tbl$group2[i])))
    ymax <- max(subdat$freq, na.rm = TRUE)
    step <- max(0.06 * ymax, 0.02)
    y_positions <- seq(ymax + step, by = step, length.out = length(comps))
    labels <- ifelse(is.na(stat_tbl$p.adj.signif), "ns", stat_tbl$p.adj.signif)
    base_plot + ggsignif::geom_signif(comparisons = comps, annotations = labels, y_position = y_positions, tip_length = 0.01, 
        vjust = 0.5, textsize = 3)
}
plot_and_save_facets <- function(df_in, suffix_label, make_title, outdir_facets, pop_levels, nature_palette, ncol_facets, 
    per_col_w, per_row_h) {
    if (nrow(df_in) == 0) {
        message("facet skipped: ", suffix_label, " has no data")
        return(invisible(NULL))
    }
    df_in <- df_in %>% mutate(population = factor(population, levels = pop_levels, ordered = TRUE), HaplogroupF = factor(Haplogroup, 
        ordered = FALSE))
    median_by_pop <- df_in %>% group_by(HaplogroupF, population) %>% summarise(freq_median = median(freq, na.rm = TRUE), 
        .groups = "drop") %>% mutate(pop_index = as.numeric(factor(population, levels = pop_levels, ordered = TRUE)))
    corr_table <- median_by_pop %>% group_by(HaplogroupF) %>% summarise(cor_trend = ifelse(n() >= 2 && sd(pop_index) > 
        0 && sd(freq_median) > 0, cor(pop_index, freq_median, method = "pearson"), NA_real_), .groups = "drop") %>% 
        arrange(is.na(cor_trend), cor_trend)
    hap_levels_sorted <- as.character(corr_table$HaplogroupF)
    df_in <- df_in %>% mutate(HaplogroupF = factor(HaplogroupF, levels = hap_levels_sorted, ordered = TRUE))
    pop_colors <- setNames(rep(nature_palette, length.out = length(pop_levels)), pop_levels)
    p_facets <- ggplot(df_in, aes(x = population, y = freq)) + ggbeeswarm::geom_quasirandom(aes(color = population), 
        method = qr_method, width = qr_width, alpha = qr_alpha, size = qr_size, show.legend = FALSE) + geom_boxplot(aes(color = population), 
        width = 0.5, outlier.shape = NA, fill = NA, size = 0.4, show.legend = FALSE) + stat_summary(fun = median, 
        geom = "crossbar", width = 0.5, fatten = 0, color = "red", size = 0.4) + scale_color_manual(values = pop_colors) + 
        labs(x = "Population", y = "Frequency", title = make_title) + theme_bw() + theme(strip.text = element_text(size = 6, 
        face = "bold"), axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 5), axis.text.y = element_text(size = 5), 
        axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7), plot.title = element_text(hjust = 0.5, 
            size = 10), panel.spacing.x = unit(0.4, "lines"), panel.spacing.y = unit(0.4, "lines")) + facet_wrap(~HaplogroupF, 
        ncol = ncol_facets, scales = "free_y")
    n_haps <- length(unique(df_in$HaplogroupF))
    n_rows <- ceiling(n_haps/ncol_facets)
    fig_w <- per_col_w * ncol_facets
    fig_h <- per_row_h * n_rows
    pdf_path <- file.path(outdir_facets, paste0("facet_", suffix_label, ".pdf"))
    save_plot(pdf_path, plot = p_facets, width = fig_w, height = fig_h, limitsize = FALSE)
    corr_outfile <- file.path(outdir_facets, paste0("haplogroup_correlation_order_", suffix_label, ".tsv"))
    write_result(corr_table %>% rename(Haplogroup = HaplogroupF), corr_outfile)
}
run_bootstrap_for_hh <- function(df_hh, HH, out_dir, n_per_sample = bootstrap_n_per_sample, n_reps = bootstrap_n_reps, 
    pop_levels_opt = bootstrap_pop_levels, ncol_facets = bootstrap_ncol_facets, per_col_w = bootstrap_per_col_w, 
    per_row_h = bootstrap_per_row_h, keep_only = bootstrap_keep_only_levels) {
    message2("Start Bootstrap (HH=%d)", HH)
    outdir_main <- file.path(out_dir, "bootstrap")
    outdir_facets <- file.path(outdir_main, "facet_pages")
    dir.create(outdir_main, showWarnings = FALSE, recursive = TRUE)
    dir.create(outdir_facets, showWarnings = FALSE, recursive = TRUE)
    df_boot <- df_hh %>% transmute(id = sample_id %||% NA, Haplogroup_full = haplogroup_full %||% NA, population = population, 
        Haplogroup = hap) %>% filter(!is.na(Haplogroup), !is.na(population))
    if (is.null(pop_levels_opt) || length(pop_levels_opt) == 0) 
        pop_levels <- get_pop_order_by_ORpeak(df_boot)
    else pop_levels <- pop_levels_opt
    if (keep_only) 
        df_boot <- df_boot %>% filter(population %in% pop_levels)
    n_total <- nrow(df_boot)
    if (n_total == 0) {
        message2("HH=%d has no rows for bootstrap, skip.", HH)
        return(invisible(NULL))
    }
    message2("HH=%d bootstrap total rows: %d", HH, n_total)
    nature_palette <- c("#4C72B0", "#55A868", "#C44E52", "#8172B2", "#CCB974", "#64B5CD", "#8C8C8C", "#E17C05", 
        "#937860", "#5A9E6F")
    pop_colors <- setNames(rep(nature_palette, length.out = length(pop_levels)), pop_levels)
    boot_list <- vector("list", n_reps)
    for (r in seq_len(n_reps)) boot_list[[r]] <- bootstrap_once_with_div(r, df_boot, n_total, n_per_sample)
    freq_long <- bind_rows(lapply(boot_list, `[[`, "freq_tbl"))
    diversity_long <- bind_rows(lapply(boot_list, `[[`, "div_tbl"))
    freq_table_tsv <- file.path(outdir_main, "haplogroup_bootstrap_freq_long.tsv")
    sig_table_tsv <- file.path(outdir_main, "haplogroup_significance_results.tsv")
    diversity_table_tsv <- file.path(outdir_main, "population_diversity_bootstrap_long.tsv")
    write_result(freq_long, freq_table_tsv, row.names = FALSE)
    write_result(diversity_long, diversity_table_tsv, row.names = FALSE)
    hap_levels_thisHH <- sort(unique(freq_long$Haplogroup))
    all_sig_rows <- list()
    for (hap in hap_levels_thisHH) {
        subdat <- freq_long %>% filter(Haplogroup == hap)
        if (nrow(subdat) == 0) 
            next
        sig_tab <- get_sig_table_for_hap(subdat, hap)
        all_sig_rows[[hap]] <- sig_tab
        p_hap <- plot_hap_single(subdat = subdat, pop_levels = pop_levels, pop_colors = pop_colors, n_per_sample = n_per_sample, 
            n_reps = n_reps)
        w <- max(3.5, 0.4 * length(pop_levels))
        h <- 3.5
        save_plot(file.path(outdir_main, paste0("freqplot_", hap, ".pdf")), plot = p_hap, width = w, height = h, 
            limitsize = FALSE)
    }
    sig_results_all <- bind_rows(all_sig_rows)
    write_result(sig_results_all, sig_table_tsv, row.names = FALSE)
    p_shannon <- plot_diversity_metric(diversity_long, metric = "shannon", pop_levels = pop_levels, pop_colors = pop_colors)
    p_simpson <- plot_diversity_metric(diversity_long, metric = "simpson", pop_levels = pop_levels, pop_colors = pop_colors)
    save_plot(file.path(outdir_main, "diversity_shannon.pdf"), plot = p_shannon, width = 6, height = 4, limitsize = FALSE)
    save_plot(file.path(outdir_main, "diversity_simpson.pdf"), plot = p_simpson, width = 6, height = 4, limitsize = FALSE)
    freq_long_forFacet <- freq_long %>% mutate(Haplogroup = as.character(Haplogroup), hap_len = nchar(Haplogroup), 
        population = factor(population, levels = pop_levels, ordered = TRUE))
    freq_lenK <- freq_long_forFacet %>% filter(hap_len == HH)
    plot_and_save_facets(df_in = freq_lenK, suffix_label = paste0("length", sprintf("%02d", HH), "_sorted"), make_title = paste0("Haplogroup bootstrap frequency (length=", 
        HH, ", sorted)"), outdir_facets = outdir_facets, pop_levels = pop_levels, nature_palette = nature_palette, 
        ncol_facets = bootstrap_ncol_facets, per_col_w = bootstrap_per_col_w, per_row_h = bootstrap_per_row_h)
    n_pops_total <- length(pop_levels)
    haps_with_all <- freq_lenK %>% group_by(Haplogroup) %>% summarise(n_pop = n_distinct(population), .groups = "drop") %>% 
        filter(n_pop == n_pops_total) %>% pull(Haplogroup)
    freq_lenK_all <- freq_lenK %>% filter(Haplogroup %in% haps_with_all)
    plot_and_save_facets(df_in = freq_lenK_all, suffix_label = paste0("length", sprintf("%02d", HH), "_sorted_allPops"), 
        make_title = paste0("Haplogroup bootstrap frequency (length=", HH, ", present in all ", n_pops_total, " groups)"), 
        outdir_facets = outdir_facets, pop_levels = pop_levels, nature_palette = nature_palette, ncol_facets = bootstrap_ncol_facets, 
        per_col_w = bootstrap_per_col_w, per_row_h = bootstrap_per_row_h)
    message2("Bootstrap done (HH=%d), output: %s", HH, basename(outdir_main))
}
