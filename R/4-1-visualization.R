source("R/1-2-runtime.R")
meta_in <- if (nzchar(metadata)) read_table(metadata) else NULL
df_in <- read_table(input)
make_multilayer_sankey <- TRUE
for (HH in HH_range) {
    f <- file.path(frequency, sprintf("freq_matrix_HH%02d.tsv", HH))
    if (!file.exists(f)) 
        next
    freq_mat_df <- read_table(f)
    count_mat_df <- read_table(file.path(frequency, sprintf("count_matrix_HH%02d.tsv", HH)))
    df <- df_in %>% mutate(hap = hap_substr(hap_clean, HH))
    for (metric in c("TVD", "JSDist", "NeiD", "DA")) {
        tab <- read_table(file.path(distance, sprintf("pairwise_%s_HH%02d.tsv", metric, HH)))
        D <- as.matrix(tab[, -1, drop = FALSE])
        rownames(D) <- tab$population
        assign(c(TVD = "D_TVD", JSDist = "D_JSD", NeiD = "D_NEI", DA = "D_DA")[[metric]], D)
    }
    plot_distance_heatmap(D_TVD, file.path(figure, sprintf("heatmap_TVD_HH%02d.pdf", HH)), title = sprintf("TVD heatmap (HH=%d)", 
        HH))
    plot_distance_heatmap(D_JSD, file.path(figure, sprintf("heatmap_JSD_HH%02d.pdf", HH)), title = sprintf("JSDist heatmap (HH=%d)", 
        HH))
    plot_distance_heatmap(D_NEI, file.path(figure, sprintf("heatmap_NeiD_HH%02d.pdf", HH)), title = sprintf("Nei's D heatmap (HH=%d)", 
        HH))
    plot_distance_heatmap(D_DA, file.path(figure, sprintf("heatmap_DA_HH%02d.pdf", HH)), title = sprintf("DA heatmap (HH=%d)", 
        HH))
    plot_dendrogram(D_TVD, file.path(figure, sprintf("dendrogram_TVD_HH%02d.pdf", HH)))
    plot_dendrogram(D_NEI, file.path(figure, sprintf("dendrogram_NeiD_HH%02d.pdf", HH)))
    plot_dendrogram(D_DA, file.path(figure, sprintf("dendrogram_DA_HH%02d.pdf", HH)))
    plot_nj_tree(D_NEI, file.path(figure, sprintf("NJ_NeiD_HH%02d.pdf", HH)), title = sprintf("NJ tree (Nei's D, HH=%d)", 
        HH))
    plot_nj_tree(D_DA, file.path(figure, sprintf("NJ_DA_HH%02d.pdf", HH)), title = sprintf("NJ tree (DA, HH=%d)", 
        HH))
    run_mds_plot(D_TVD, file.path(figure, sprintf("MDS_TVD_HH%02d.pdf", HH)))
    run_mds_plot(D_JSD, file.path(figure, sprintf("MDS_JSD_HH%02d.pdf", HH)))
    run_mds_plot(D_NEI, file.path(figure, sprintf("MDS_NeiD_HH%02d.pdf", HH)))
    run_mds_plot(D_DA, file.path(figure, sprintf("MDS_DA_HH%02d.pdf", HH)))
    run_pca_plot(freq_mat_df, file.path(figure, sprintf("PCA_HH%02d.pdf", HH)), file.path(output, sprintf("PCA_scores_HH%02d.tsv", 
        HH)))
    run_ca_plot(count_mat_df, file.path(figure, sprintf("CA_HH%02d.pdf", HH)), file.path(output, sprintf("CA_coords_HH%02d.tsv", 
        HH)))
    plot_stacked_bars(freq_mat_df, file.path(figure, sprintf("StackedBars_HH%02d.pdf", HH)))
    pop_order_opt <- order_pop_by_dominant_hap(freq_mat_df)
    sankey_pure_file <- file.path(figure, sprintf("Sankey_pure_HH%02d.pdf", HH))
    plot_sankey_pop_series_pure(freq_mat_df, out_pdf = sankey_pure_file, pop_order = pop_order_opt, label_threshold = 0.01, 
        bar_halfwidth = 0.32, plot_title = sprintf("Haplogroup flows across populations (pure, HH=%d)", HH))
    sankey_bar_file <- file.path(figure, sprintf("Sankey_withBars_HH%02d.pdf", HH))
    plot_sankey_pop_series_withbars(freq_mat_df, out_pdf = sankey_bar_file, pop_order = pop_order_opt, label_threshold = 0.01, 
        bar_halfwidth = 0.33, plot_title = sprintf("Haplogroup flows across populations (with bars, HH=%d)", HH))
    plot_alluvial_simple(df, hap_col = "hap", pop_col = "population", out_pdf = file.path(figure, sprintf("Alluvial_simple_HH%02d.pdf", 
        HH)))
    plot_alluvial_simple_freq(freq_mat_df, out_pdf = file.path(figure, sprintf("Alluvial_simpleFreq_HH%02d.pdf", 
        HH)))
    if (make_multilayer_sankey && HH >= 2) {
        HC <- sort(unique(pmin(seq_len(HH), multilayer_max_k)))
        mdat <- make_multilayer_long(df, HH_chain = HC)
        plot_alluvial_multilayer_counts_from_long(mdat, out_pdf = file.path(figure, sprintf("Alluvial_multilayer_chain_toHH%02d.pdf", 
            HH)))
        plot_alluvial_multilayer_freq_from_long(mdat, out_pdf = file.path(figure, sprintf("Alluvial_multilayerFreq_chain_toHH%02d.pdf", 
            HH)))
        root_freq_tab <- mdat$long %>% group_by(.data[[mdat$axes[1]]]) %>% summarise(global_freq = sum(freq), .groups = "drop") %>% 
            arrange(desc(global_freq))
        root_to_plot <- root_freq_tab %>% filter(global_freq >= perhap_min_freq) %>% slice_head(n = perhap_max_plot) %>% 
            pull(1)
        if (length(root_to_plot)) {
            for (rh in root_to_plot) {
                out_raw <- file.path(figure, sprintf("Alluvial_multilayerFreq_chain_toHH%02d_root_%s_raw.pdf", 
                  HH, rh))
                plot_alluvial_multilayer_freq_singlehap_raw(mdat, rh, out_raw)
                out_scaled <- file.path(figure, sprintf("Alluvial_multilayerFreq_chain_toHH%02d_root_%s_scaled.pdf", 
                  HH, rh))
                plot_alluvial_multilayer_freq_singlehap_scaled(mdat, rh, out_scaled)
            }
        }
    }
}
if (file.exists(file.path(distance, "pairwise_Rst.tsv"))) {
    tab <- read_table(file.path(distance, "pairwise_Rst.tsv"))
    D <- as.matrix(tab[, -1, drop = FALSE])
    rownames(D) <- tab$population
    plot_distance_heatmap(D, file.path(figure, "heatmap_Rst.pdf"), title = "Y-STR Rst")
    plot_dendrogram(D, file.path(figure, "dendrogram_Rst.pdf"))
    plot_nj_tree(D, file.path(figure, "NJ_Rst.pdf"), title = "Y-STR Rst")
    run_mds_plot(D, file.path(figure, "MDS_Rst.pdf"))
}
