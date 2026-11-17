source("R/1-2-runtime.R")
if (!nzchar(metadata)) {
    message("Metadata input omitted; skipping geography")
    quit(status = 0)
}
meta <- read_table(metadata)
if (!all(c("population", "latitude", "longitude") %in% names(meta))) stop("Missing geographic fields")
if (anyDuplicated(meta$population) || anyNA(meta[c("latitude", "longitude")]) || any(abs(meta$latitude) > 90) || 
    any(abs(meta$longitude) > 180)) stop("Invalid geographic metadata")
if (nzchar(population_order)) {
    pop_order <- unique(trimws(readLines(population_order, warn = FALSE)))
    pop_order <- pop_order[nzchar(pop_order)]
    if (!all(pop_order %in% meta$population)) 
        stop("Unknown populations in population order")
    meta <- meta[match(c(pop_order, setdiff(meta$population, pop_order)), meta$population), , drop = FALSE]
}
pop_coords <- meta
pop_coords$n_samples <- as.numeric(table(factor(read_table(input)$population, levels = meta$population)))
pop_coords <- adjust_bar_positions(pop_coords, bar_halfwidth_deg = default_bar_halfwidth_deg, bar_height_deg = default_bar_height_deg)
pop_coords <- separate_close_longitudes(pop_coords, min_diff = 1, frac_shift = 0.05, shift_deg = NA_real_, passes = 2)
all_pops <- pop_coords$population
geoD <- geo_dist_matrix_km(pop_coords$longitude, pop_coords$latitude)
rownames(geoD) <- colnames(geoD) <- all_pops
for (HH in HH_range) {
    f <- file.path(frequency, sprintf("freq_matrix_HH%02d.tsv", HH))
    if (!file.exists(f)) 
        next
    freq_mat_df2 <- read_table(f)
    if (!all(freq_mat_df2$population %in% all_pops)) 
        stop("Missing population coordinates")
    M <- as.matrix(freq_mat_df2[, setdiff(names(freq_mat_df2), "population"), drop = FALSE])
    rownames(M) <- freq_mat_df2$population
    M <- row_normalize(M)
    pops_now <- rownames(M)
    pop_df_now <- pop_coords %>% dplyr::filter(population %in% pops_now) %>% dplyr::arrange(match(population, pops_now))
    within_df <- within_pop_consistency(M, eps = eps0)
    message2("HH=%d computing similarity matrices (n=%d, K=%d)...", HH, nrow(M), ncol(M))
    sim_cos <- cosine_sim_matrix(M)
    sim_js <- js_similarity_matrix(M, eps = eps0)
    rownames(sim_cos) <- pops_now
    colnames(sim_cos) <- pops_now
    rownames(sim_js) <- pops_now
    colnames(sim_js) <- pops_now
    geoD_now <- geoD[pops_now, pops_now, drop = FALSE]
    knn_cos <- knn_mean_similarity(sim_cos, geoD_now, k = k_geo)
    knn_js <- knn_mean_similarity(sim_js, geoD_now, k = k_geo)
    metrics_df <- pop_df_now %>% dplyr::bind_cols(within_df) %>% dplyr::mutate(kNN_mean_cosine = knn_cos, kNN_mean_JS = knn_js, 
        HH = HH, n_components = ncol(M))
    write_result(metrics_df, file.path(output, sprintf("pop_metrics_HH%02d.tsv", HH)), row.names = FALSE, fileEncoding = "UTF-8")
    write_result(sim_cos, file.path(output, sprintf("similarity_cosine_HH%02d.tsv", HH)), fileEncoding = "UTF-8")
    write_result(sim_js, file.path(output, sprintf("similarity_JS_HH%02d.tsv", HH)), fileEncoding = "UTF-8")
    plot_heatmap_matrix_clustered(sim_cos, out_prefix = file.path(output, sprintf("Heatmap_cosine_clustered_HH%02d", 
        HH)), title = sprintf("Cosine similarity (clustered, HH=%d)", HH), is_similarity = TRUE, cluster_method = cluster_method, 
        show_dendrogram = TRUE, expected_names = pops_now)
    plot_heatmap_matrix_clustered(sim_js, out_prefix = file.path(output, sprintf("Heatmap_JS_clustered_HH%02d", 
        HH)), title = sprintf("Jensen-Shannon similarity (clustered, HH=%d)", HH), is_similarity = TRUE, cluster_method = cluster_method, 
        show_dendrogram = TRUE, expected_names = pops_now)
    plot_metric_map(metrics_df, "Consistency_1minusEntropy", out_prefix = file.path(output, sprintf("Map_1minusEntropy_HH%02d", 
        HH)), title = sprintf("Within-pop consistency (1 - normalized entropy) (HH=%d)", HH))
    plot_metric_map(metrics_df, "kNN_mean_cosine", out_prefix = file.path(output, sprintf("Map_kNNmeanCosine_HH%02d", 
        HH)), title = sprintf("kNN mean cosine similarity (k=%d, HH=%d)", k_geo, HH))
    plot_metric_map(metrics_df, "kNN_mean_JS", out_prefix = file.path(output, sprintf("Map_kNNmeanJS_HH%02d", HH)), 
        title = sprintf("kNN mean Jensen-Shannon similarity (k=%d, HH=%d)", k_geo, HH))
    plot_edges_fierce(pop_df = metrics_df, sim_mat = sim_cos, out_prefix = file.path(output, sprintf("Map_edges_FIERCE_cosine_HH%02d", 
        HH)), title = sprintf("Edges (cosine): color=sim, width=pair consistency (HH=%d)", HH), edge_mode = edge_mode, 
        top_n = top_edge_per_pop, min_sim = edge_min_sim, max_edges = max_edges, cons_col = cons_col_for_edges)
    plot_edges_fierce(pop_df = metrics_df, sim_mat = sim_js, out_prefix = file.path(output, sprintf("Map_edges_FIERCE_JS_HH%02d", 
        HH)), title = sprintf("Edges (JS): color=sim, width=pair consistency (HH=%d)", HH), edge_mode = edge_mode, 
        top_n = top_edge_per_pop, min_sim = edge_min_sim, max_edges = max_edges, cons_col = cons_col_for_edges)
    allpairs_file <- file.path(output, sprintf("GeoSankey_allpairs_curved_HH%02d.pdf", HH))
    plot_geo_sankey_allpairs(freq_mat_df = freq_mat_df2, pop_coords = pop_coords, out_pdf = allpairs_file, bar_halfwidth_deg = default_bar_halfwidth_deg, 
        bar_height_deg = default_bar_height_deg, flow_min_freq = 0.03, curve_factor = 0.35, bezier_points = 16, 
        plot_title = sprintf("Geo-Sankey haplogroup flows (all-pairs, HH = %d, freq >= %.2f)", HH, sankey_component_min_freq), 
        component_min_freq = sankey_component_min_freq, label_min_freq = sankey_label_min_freq)
    for (tp in all_pops) {
        star_file <- file.path(output, sprintf("GeoSankey_star_%s_HH%02d.pdf", tp, HH))
        plot_geo_sankey_star(freq_mat_df = freq_mat_df2, pop_coords = pop_coords, target_population = tp, out_pdf = star_file, 
            bar_halfwidth_deg = default_bar_halfwidth_deg, bar_height_deg = default_bar_height_deg, flow_min_freq = 0.03, 
            curve_factor = 0.35, bezier_points = 16, plot_title = sprintf("Geo-Sankey star flows (target = %s, HH = %d, freq >= %.2f)", 
                tp, HH, sankey_component_min_freq), component_min_freq = sankey_component_min_freq, label_min_freq = sankey_label_min_freq)
    }
    message2("HH = %d done.", HH)
}
