`%||%` <- function(a, b) if (!is.null(a)) a else b
message2 <- function(fmt, ...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(fmt, ...), "\n")
sanitize_label_ascii <- function(x) {
    x <- as.character(x)
    y <- iconv(x, from = "", to = "ASCII//TRANSLIT", sub = "_")
    y[is.na(y)] <- "_"
    y
}
gsub_bytes <- function(pattern, replacement, x) gsub(pattern, replacement, x, useBytes = TRUE)
sanitize_hap <- function(x) {
    x <- stringr::str_replace_all(x, "[^A-Za-z0-9]", "")
    x <- tolower(x)
    if (nchar(x[1]) > 0) 
        substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
}
hap_substr <- function(h, HH) substr(h, 1, HH)
safe_group_from_meta <- function(meta_in) {
    if (is.null(meta_in) || !is.data.frame(meta_in) || ncol(meta_in) == 0) 
        return(NULL)
    meta_lc <- meta_in
    colnames(meta_lc) <- tolower(colnames(meta_lc))
    if (!("population" %in% colnames(meta_lc))) {
        cn <- tolower(names(meta_lc))
        pid <- which(cn %in% c("population", "pop", "group_name"))
        if (length(pid) == 1) 
            names(meta_lc)[pid] <- "population"
    }
    if ("group" %in% colnames(meta_lc)) {
        out <- meta_lc[, c("population", "group"), drop = FALSE]
    }
    else if ("region" %in% colnames(meta_lc)) {
        out <- meta_lc[, c("population", "region"), drop = FALSE]
        names(out)[2] <- "group"
    }
    else return(NULL)
    out$population <- sanitize_label_ascii(out$population)
    out$group <- sanitize_label_ascii(out$group)
    out
}
.okabe_ito <- c("#0072B2", "#56B4E9", "#009E73", "#E69F00", "#CC79A7", "#D55E00", "#999999", "#F0E442")
.nature_seeds <- c("#4477AA", "#66CCEE", "#228833", "#CCBB44", "#EE6677", "#AA3377", "#332288", "#88CCEE", "#44AA99", 
    "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499", "#DDDDDD", "#1F77B4", "#FF7F0E", "#2CA02C", 
    "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF")
.seeds <- unique(c(.okabe_ito, .nature_seeds))
.hex_to_lab <- function(hex) {
    hex <- as.character(hex)
    rgb <- t(grDevices::col2rgb(hex)/255)
    lab <- grDevices::convertColor(rgb, from = "sRGB", to = "Lab", scale.in = 1)
    lab <- as.matrix(lab)
    colnames(lab) <- c("L", "a", "b")
    lab
}
.gen_candidates <- function() {
    levL <- c(80, 72, 64, 56, 48)
    levC <- c(60, 66, 62, 58, 54)
    out <- character()
    for (i in seq_along(levL)) {
        cols <- grDevices::hcl(h = seq(0, 359, by = 1.5), c = levC[i], l = levL[i])
        out <- c(out, cols)
    }
    unique(out)
}
glasbey_like <- function(n, seed_colors = .seeds, max_n = 2048) {
    n <- min(n, max_n)
    seed_colors <- unique(seed_colors)
    if (n <= length(seed_colors)) 
        return(seed_colors[seq_len(n)])
    cand <- setdiff(.gen_candidates(), seed_colors)
    sel <- seed_colors
    sel_lab <- .hex_to_lab(sel)
    cand_lab <- .hex_to_lab(cand)
    saved_rng <- .Random.seed
    on.exit(assign(".Random.seed", saved_rng, envir = .GlobalEnv))
    set.seed(seed)
    ord <- sample.int(nrow(cand_lab))
    cand_lab <- cand_lab[ord, , drop = FALSE]
    cand <- cand[ord]
    min_dist_to_sel <- function(row, sel_lab) {
        min(sqrt(rowSums((t(sel_lab) - row)^2)))
    }
    dist_min <- apply(cand_lab, 1, min_dist_to_sel, sel_lab = sel_lab)
    while (length(sel) < n && nrow(cand_lab) > 0) {
        pick <- which.max(dist_min)
        sel <- c(sel, cand[pick])
        new_row <- cand_lab[pick, , drop = FALSE]
        sel_lab <- rbind(sel_lab, new_row)
        dnew <- sqrt(rowSums((cand_lab - matrix(new_row, nrow = nrow(cand_lab), ncol = 3, byrow = TRUE))^2))
        dist_min <- pmin(dist_min, dnew)
        cand_lab <- cand_lab[-pick, , drop = FALSE]
        cand <- cand[-pick]
        dist_min <- dist_min[-pick]
    }
    sel[seq_len(n)]
}
get_palette_for_levels <- function(fac_levels, min_n = 0) {
    lv <- as.character(unique(fac_levels))
    n <- max(length(lv), min_n)
    pal <- glasbey_like(n)
    setNames(pal[seq_len(length(lv))], lv)
}
nature_hm_seq <- function(n = 100) {
    (grDevices::colorRampPalette(c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", 
        "#08519C", "#08306B")))(n)
}
save_pdf <- function(p, file_path_noext, width, height) {
    suppressWarnings(save_plot(paste0(file_path_noext, ".pdf"), p, width = width, height = height, limitsize = FALSE))
}
order_pop_by_dominant_hap <- function(freq_mat_df) {
    long <- freq_mat_df %>% rename(population = 1) %>% pivot_longer(-population, names_to = "hap", values_to = "freq")
    long$hap <- droplevels(factor(long$hap))
    pop_top <- long %>% group_by(population) %>% slice_max(freq, n = 1, with_ties = FALSE) %>% ungroup() %>% rename(top_hap = hap, 
        top_freq = freq)
    hap_rank <- pop_top %>% group_by(top_hap) %>% summarise(max_freq_any_pop = max(top_freq, na.rm = TRUE), .groups = "drop") %>% 
        arrange(desc(max_freq_any_pop)) %>% mutate(hap_order = row_number())
    pop_order_df <- pop_top %>% left_join(hap_rank, by = "top_hap") %>% arrange(hap_order, desc(top_freq), population)
    pop_order_df$population
}
plot_distance_heatmap <- function(D, out_file, title = "Distance heatmap") {
    D <- sanitize_dist_matrix(D, name = title)
    cols <- nature_hm_seq(100)
    base <- tools::file_path_sans_ext(out_file)
    hc <- hclust(as.dist(D), method = "ward.D2")
    ord <- hc$order
    D2 <- D[ord, ord, drop = FALSE]
    open_pdf(paste0(base, ".pdf"), width = plot_base_width, height = plot_base_height)
    suppressWarnings(pheatmap::pheatmap(D2, color = cols, cluster_rows = FALSE, cluster_cols = FALSE, main = title, 
        border_color = NA, legend = TRUE, show_rownames = TRUE, show_colnames = TRUE))
    dev.off()
}
plot_dendrogram <- function(D, out_file) {
    D <- sanitize_dist_matrix(D, name = out_file)
    hc <- hclust(as.dist(D), method = "ward.D2")
    base <- tools::file_path_sans_ext(out_file)
    open_pdf(paste0(base, ".pdf"), width = plot_base_width, height = plot_base_height)
    plot(hc, xlab = "", sub = "", main = "Hierarchical clustering (ward.D2)")
    dev.off()
    phy <- ape::as.phylo(hc)
    phy$tip.label <- sanitize_label_ascii(gsub_bytes("[^A-Za-z0-9_.-]", "_", phy$tip.label))
    ape::write.tree(phy, file = paste0(base, ".nwk"))
    message2("Exported Newick tree: %s", basename(paste0(base, ".nwk")))
}
plot_nj_tree <- function(D, out_file, title = "Neighbor-Joining tree") {
    D <- sanitize_dist_matrix(D, name = out_file)
    tr <- ape::nj(as.dist(D))
    base <- tools::file_path_sans_ext(out_file)
    open_pdf(paste0(base, ".pdf"), width = plot_base_width, height = plot_base_height)
    plot(tr, main = title, cex = 0.8)
    dev.off()
    tr$tip.label <- sanitize_label_ascii(gsub_bytes("[^A-Za-z0-9_.-]", "_", tr$tip.label))
    ape::write.tree(tr, file = paste0(base, ".nwk"))
    message2("Exported NJ Newick: %s", basename(paste0(base, ".nwk")))
}
plot_stacked_bars <- function(freq_mat_df, out_file) {
    long <- freq_mat_df %>% rename(population = 1) %>% pivot_longer(-population, names_to = "hap", values_to = "freq")
    long$hap <- droplevels(factor(long$hap))
    pal <- get_palette_for_levels(levels(long$hap), min_n = max(256, nlevels(long$hap)))
    p <- ggplot(long, aes(x = population, y = freq, fill = hap)) + geom_bar(stat = "identity", width = 0.8) + scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
        scale_fill_manual(values = pal, limits = names(pal), drop = FALSE, na.value = "#BDBDBD") + theme_bw(base_size = 12) + 
        theme(axis.text.x = element_text(angle = 60, hjust = 1), panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) + 
        labs(x = "Population", y = "Frequency", fill = "Haplogroup")
    base <- tools::file_path_sans_ext(out_file)
    save_pdf(p, base, width = plot_base_width * 1.6, height = plot_base_height * 1.2)
}
run_pca_plot <- function(freq_mat_df, out_pdf, out_scores_tsv) {
    mat <- as.matrix(freq_mat_df[, -1])
    rownames(mat) <- sanitize_label_ascii(freq_mat_df[[1]])
    pca <- prcomp(mat, center = TRUE, scale. = FALSE)
    var_exp <- (pca$sdev^2)/sum(pca$sdev^2)
    scores <- data.frame(population = rownames(pca$x), pca$x[, seq_len(min(5, ncol(pca$x))), drop = FALSE], check.names = FALSE, 
        row.names = NULL)
    grp <- safe_group_from_meta(if (exists("meta_in", inherits = TRUE)) 
        meta_in
    else NULL)
    if (!is.null(grp)) {
        scores <- dplyr::left_join(scores, grp, by = "population")
    }
    else {
        hc <- hclust(dist(mat), method = "ward.D2")
        cl <- cutree(hc, k = min(4, nrow(mat)))
        scores$group <- factor(cl)
    }
    write_result(scores, out_scores_tsv, row.names = FALSE, fileEncoding = "UTF-8")
    pc1 <- sprintf("PC1 (%.1f%%)", 100 * var_exp[1])
    pc2 <- sprintf("PC2 (%.1f%%)", 100 * var_exp[2])
    cols <- get_palette_for_levels(levels(factor(scores$group)))
    p <- ggplot(scores, aes(PC1, PC2, label = population, color = group)) + geom_point(size = 3) + geom_text_repel(hjust = 0, 
        max.overlaps = 100) + scale_color_manual(values = cols, limits = names(cols), drop = FALSE, na.value = "#BDBDBD") + 
        theme_bw(base_size = 12) + theme(panel.grid = element_blank()) + labs(title = "PCA of haplogroup frequencies", 
        x = pc1, y = pc2, color = "Group")
    base <- tools::file_path_sans_ext(out_pdf)
    save_pdf(p, base, width = plot_base_width, height = plot_base_height)
}
run_ca_plot <- function(count_mat_df, out_pdf, out_coords_tsv) {
    mat <- as.data.frame(count_mat_df)
    rownames(mat) <- sanitize_label_ascii(mat[[1]])
    mat <- mat[, -1, drop = FALSE]
    ca <- FactoMineR::CA(mat, graph = FALSE)
    row_coords <- data.frame(population = rownames(ca$row$coord), ca$row$coord, check.names = FALSE, row.names = NULL)
    col_coords <- data.frame(hap = rownames(ca$col$coord), ca$col$coord, check.names = FALSE, row.names = NULL)
    eig <- as.data.frame(ca$eig)
    row_contrib <- data.frame(population = rownames(ca$row$contrib), ca$row$contrib, check.names = FALSE, row.names = NULL)
    col_contrib <- data.frame(hap = rownames(ca$col$contrib), ca$col$contrib, check.names = FALSE, row.names = NULL)
    write_result(row_coords, sub("coords", "row_coords", out_coords_tsv), row.names = FALSE, fileEncoding = "UTF-8")
    write_result(col_coords, sub("coords", "col_coords", out_coords_tsv), row.names = FALSE, fileEncoding = "UTF-8")
    write_result(eig, sub("coords", "eigenvalues", out_coords_tsv), row.names = FALSE, fileEncoding = "UTF-8")
    write_result(row_contrib, sub("coords", "row_contrib", out_coords_tsv), row.names = FALSE, fileEncoding = "UTF-8")
    write_result(col_contrib, sub("coords", "col_contrib", out_coords_tsv), row.names = FALSE, fileEncoding = "UTF-8")
    scores <- row_coords
    var_exp <- eig[, 2]
    grp <- safe_group_from_meta(if (exists("meta_in", inherits = TRUE)) 
        meta_in
    else NULL)
    if (!is.null(grp)) {
        scores <- dplyr::left_join(scores, grp, by = "population")
    }
    else {
        hc <- hclust(dist(as.matrix(mat)), method = "ward.D2")
        cl <- cutree(hc, k = min(4, nrow(mat)))
        scores$group <- factor(cl)
    }
    rows12 <- scores[, c("population", "Dim 1", "Dim 2", "group")]
    colnames(rows12) <- c("population", "Dim1", "Dim2", "group")
    cols12 <- col_coords[, c("hap", "Dim 1", "Dim 2")]
    colnames(cols12) <- c("hap", "Dim1", "Dim2")
    pc1 <- sprintf("Dim1 (%.1f%%)", var_exp[1])
    pc2 <- sprintf("Dim2 (%.1f%%)", var_exp[2])
    cols_grp <- get_palette_for_levels(levels(factor(rows12$group)))
    p <- ggplot() + geom_point(data = rows12, aes(Dim1, Dim2, color = group), size = 3) + geom_text_repel(data = rows12, 
        aes(Dim1, Dim2, label = population, color = group), hjust = 0, max.overlaps = 100) + geom_segment(data = cols12, 
        aes(x = 0, y = 0, xend = Dim1, yend = Dim2), arrow = arrow(length = unit(0.02, "npc")), alpha = 0.6) + 
        geom_text_repel(data = cols12, aes(Dim1, Dim2, label = hap), size = 3) + scale_color_manual(values = cols_grp, 
        limits = names(cols_grp), drop = FALSE, na.value = "#BDBDBD") + theme_bw(base_size = 12) + theme(panel.grid = element_blank()) + 
        labs(title = "CA biplot (rows=populations, cols=haplogroups)", x = pc1, y = pc2, color = "Group")
    base <- tools::file_path_sans_ext(out_pdf)
    save_pdf(p, base, width = plot_base_width, height = plot_base_height)
}
run_mds_plot <- function(D, out_pdf) {
    D <- sanitize_dist_matrix(D, name = out_pdf)
    conf <- cmdscale(as.dist(D), k = 2, eig = TRUE)
    pts <- as.data.frame(conf$points)
    pts$population <- sanitize_label_ascii(rownames(pts))
    colnames(pts)[1:2] <- c("Dim1", "Dim2")
    eig <- conf$eig
    pos <- eig[eig > 0]
    var_exp <- if (length(pos)) 
        pos/sum(pos)
    else c(NA_real_, NA_real_)
    var_lab <- c(sprintf("Dim1 (%.1f%%)", 100 * ifelse(length(var_exp) >= 1, var_exp[1], NA)), sprintf("Dim2 (%.1f%%)", 
        100 * ifelse(length(var_exp) >= 2, var_exp[2], NA)))
    grp <- safe_group_from_meta(if (exists("meta_in", inherits = TRUE)) 
        meta_in
    else NULL)
    if (!is.null(grp)) {
        pts <- dplyr::left_join(pts, grp, by = "population")
    }
    else {
        hc <- hclust(as.dist(D), method = "ward.D2")
        cl <- cutree(hc, k = min(4, nrow(D)))
        pts$group <- factor(cl)
    }
    base_tsv <- tools::file_path_sans_ext(out_pdf)
    write_result(pts, paste0(base_tsv, "_MDS_points.tsv"), row.names = FALSE, fileEncoding = "UTF-8")
    write_result(data.frame(eigenvalue = eig, component = seq_along(eig), row.names = NULL), paste0(base_tsv, "_MDS_eigenvalues.tsv"), 
        row.names = FALSE, fileEncoding = "UTF-8")
    cols <- get_palette_for_levels(levels(factor(pts$group)))
    p <- ggplot(pts, aes(Dim1, Dim2, label = population, color = group)) + geom_point(size = 3) + geom_text_repel(hjust = 0, 
        max.overlaps = 100) + scale_color_manual(values = cols, limits = names(cols), drop = FALSE, na.value = "#BDBDBD") + 
        theme_bw(base_size = 12) + theme(panel.grid = element_blank()) + labs(title = "MDS (cmdscale) on distance matrix", 
        x = var_lab[1], y = var_lab[2], color = "Group")
    base <- tools::file_path_sans_ext(out_pdf)
    save_pdf(p, base, width = plot_base_width, height = plot_base_height)
}
plot_alluvial_simple <- function(df, hap_col, pop_col, out_pdf) {
    tab <- df %>% count(.data[[hap_col]], .data[[pop_col]], name = "n") %>% rename(hap = 1, population = 2)
    tab$hap <- droplevels(factor(tab$hap))
    tab$population <- factor(sanitize_label_ascii(tab$population))
    pal <- get_palette_for_levels(levels(tab$hap), min_n = max(512, nlevels(tab$hap)))
    p <- ggplot(tab, aes(y = n, axis1 = hap, axis2 = population)) + geom_alluvium(aes(fill = hap), width = 1/12, 
        knot.pos = 0.4, alpha = 0.9) + geom_stratum(aes(fill = hap), width = 1/12, color = "grey30") + geom_text(stat = "stratum", 
        aes(label = after_stat(stratum)), hjust = 0.5, vjust = 0.5, size = 3) + scale_x_discrete(limits = c("Haplogroup", 
        "Population"), expand = c(0.1, 0.05)) + scale_fill_manual(values = pal, limits = names(pal), drop = TRUE, 
        na.value = "#BDBDBD") + theme_bw(base_size = 12) + theme(panel.grid = element_blank()) + labs(y = "Count", 
        x = "", fill = "Haplogroup", title = "Alluvial: Haplogroup → Population (Counts)")
    base <- tools::file_path_sans_ext(out_pdf)
    save_pdf(p, base, width = plot_base_width * 1.8, height = plot_base_height * 1.3)
}
plot_alluvial_simple_freq <- function(freq_mat_df, out_pdf) {
    long <- freq_mat_df %>% rename(population = 1) %>% pivot_longer(-population, names_to = "hap", values_to = "freq")
    long$hap <- droplevels(factor(long$hap))
    long$population <- factor(sanitize_label_ascii(long$population))
    pal <- get_palette_for_levels(levels(long$hap), min_n = max(512, nlevels(long$hap)))
    p <- ggplot(long, aes(y = freq, axis1 = hap, axis2 = population)) + geom_alluvium(aes(fill = hap), width = 1/12, 
        knot.pos = 0.4, alpha = 0.9) + geom_stratum(aes(fill = hap), width = 1/12, color = "grey30") + geom_text(stat = "stratum", 
        aes(label = after_stat(stratum)), hjust = 0.5, vjust = 0.5, size = 3) + scale_x_discrete(limits = c("Haplogroup", 
        "Population"), expand = c(0.1, 0.05)) + scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
        scale_fill_manual(values = pal, limits = names(pal), drop = TRUE, na.value = "#BDBDBD") + theme_bw(base_size = 12) + 
        theme(panel.grid = element_blank()) + labs(y = "Frequency", x = "", fill = "Haplogroup", title = "Alluvial: Haplogroup → Population (Frequency)")
    base <- tools::file_path_sans_ext(out_pdf)
    save_pdf(p, base, width = plot_base_width * 1.8, height = plot_base_height * 1.3)
}
make_multilayer_long <- function(df, HH_chain = c(1, 2, 3, 4)) {
    tmp <- df
    for (k in HH_chain) tmp[[paste0("hap", k)]] <- hap_substr(tmp$hap_clean, k)
    axes <- c(paste0("hap", HH_chain), "population")
    long <- tmp %>% dplyr::count(dplyr::across(dplyr::all_of(axes)), name = "n") %>% dplyr::group_by(population) %>% 
        dplyr::mutate(freq = n/sum(n)) %>% dplyr::ungroup()
    for (ax in axes) long[[ax]] <- as.character(long[[ax]])
    stopifnot(all(vapply(long[axes], is.character, TRUE)))
    long$population <- as.character(long$population)
    deepest_col <- axes[length(axes) - 1]
    deepest_levels <- sort(unique(long[[deepest_col]]))
    pal_deep <- get_palette_for_levels(deepest_levels, min_n = max(1024, length(deepest_levels)))
    lighten_col <- function(col, f = 0.3) {
        v <- grDevices::col2rgb(col)/255
        v <- v + (1 - v) * f
        rgb(v[1], v[2], v[3])
    }
    hap_axes <- axes[-length(axes)]
    n_hap_axes <- length(hap_axes)
    pal_all <- pal_deep
    if (n_hap_axes > 1) {
        for (i in seq_len(n_hap_axes - 1)) {
            ax <- hap_axes[i]
            levs_ax <- sort(unique(long[[ax]]))
            for (lv in levs_ax) {
                cand <- deepest_levels[startsWith(deepest_levels, lv)]
                base_col <- if (!length(cand)) 
                  "#BDBDBD"
                else if (length(cand) == 1) {
                  pal_deep[cand]
                }
                else {
                  cols <- pal_deep[cand]
                  rgbs <- sapply(cols, function(z) grDevices::col2rgb(z)/255)
                  m <- rowMeans(rgbs)
                  rgb(m[1], m[2], m[3])
                }
                frac <- (n_hap_axes - i)/(n_hap_axes - 1 + 1e-09)
                pal_all[[lv]] <- lighten_col(base_col, f = 0.35 * frac)
            }
        }
    }
    pop_levels <- sort(unique(long$population))
    pal_all <- c(pal_all, setNames(rep("#E0E0E0", length(pop_levels)), pop_levels))
    list(long = long, axes = axes, deepest = deepest_col, deepest_levels = deepest_levels, pal_all = pal_all, pop_levels = pop_levels)
}
.compute_label_keep <- function(long, hap_axes, label_min_freq = 0.01) {
    keep <- character(0)
    for (ax in hap_axes) {
        tmp <- long %>% dplyr::group_by(population, lvl = .data[[ax]]) %>% dplyr::summarise(f = sum(freq, na.rm = TRUE), 
            .groups = "drop") %>% dplyr::group_by(lvl) %>% dplyr::summarise(avg_freq = mean(f, na.rm = TRUE), .groups = "drop")
        keep <- union(keep, as.character(tmp$lvl[tmp$avg_freq >= label_min_freq]))
    }
    keep[!is.na(keep) & nzchar(keep)]
}
plot_alluvial_multilayer_counts_from_long <- function(mdat, out_pdf, label_min_freq = 0.01) {
    long <- mdat$long
    axes <- mdat$axes
    deepest <- mdat$deepest
    pal_all <- mdat$pal_all
    hap_axes <- axes[-length(axes)]
    pop_levels <- mdat$pop_levels
    deep_levels <- mdat$deepest_levels
    label_keep <- .compute_label_keep(long, hap_axes, label_min_freq)
    mapping <- do.call(ggplot2::aes, c(list(y = quote(n)), setNames(lapply(axes, as.name), paste0("axis", seq_along(axes)))))
    p <- ggplot2::ggplot(long, mapping) + ggalluvial::geom_alluvium(ggplot2::aes(fill = .data[[deepest]]), width = 1/14, 
        knot.pos = 0.38, alpha = 0.85) + ggalluvial::geom_stratum(ggplot2::aes(fill = after_stat(stratum)), width = 1/14, 
        color = "grey30", alpha = 0.98, show.legend = FALSE) + ggplot2::geom_text(stat = "stratum", ggplot2::aes(label = dplyr::if_else(after_stat(stratum) %in% 
        !!pop_levels, after_stat(stratum), dplyr::if_else(after_stat(stratum) %in% !!label_keep, after_stat(stratum), 
        NA_character_))), size = 3, vjust = 0.5, hjust = 0.5) + ggplot2::scale_x_discrete(limits = c(paste0("HH=", 
        sub("hap", "", hap_axes)), "Population"), expand = c(0.1, 0.05)) + ggplot2::scale_fill_manual(values = pal_all, 
        breaks = deep_levels, drop = TRUE, na.value = "#BDBDBD") + ggplot2::guides(fill = ggplot2::guide_legend(title = "Haplogroup")) + 
        ggplot2::theme_bw(base_size = 12) + ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 10), 
        plot.margin = ggplot2::margin(10, 40, 10, 10)) + ggplot2::coord_cartesian(clip = "off") + ggplot2::labs(y = "Count", 
        x = "", title = "Multilayer Alluvial (Counts)")
    base <- tools::file_path_sans_ext(out_pdf)
    suppressWarnings(save_plot(paste0(base, ".pdf"), p, width = plot_base_width * 2, height = plot_base_height * 
        1.25, limitsize = FALSE))
}
plot_alluvial_multilayer_freq_from_long <- function(mdat, out_pdf, label_min_freq = 0.01) {
    long <- mdat$long
    axes <- mdat$axes
    deepest <- mdat$deepest
    pal_all <- mdat$pal_all
    hap_axes <- axes[-length(axes)]
    pop_levels <- mdat$pop_levels
    deep_levels <- mdat$deepest_levels
    label_keep <- .compute_label_keep(long, hap_axes, label_min_freq)
    mapping <- do.call(ggplot2::aes, c(list(y = quote(freq)), setNames(lapply(axes, as.name), paste0("axis", seq_along(axes)))))
    p <- ggplot2::ggplot(long, mapping) + ggalluvial::geom_alluvium(ggplot2::aes(fill = .data[[deepest]]), width = 1/14, 
        knot.pos = 0.38, alpha = 0.85) + ggalluvial::geom_stratum(ggplot2::aes(fill = after_stat(stratum)), width = 1/14, 
        color = "grey30", alpha = 0.98, show.legend = FALSE) + ggplot2::geom_text(stat = "stratum", ggplot2::aes(label = dplyr::if_else(after_stat(stratum) %in% 
        !!pop_levels, after_stat(stratum), dplyr::if_else(after_stat(stratum) %in% !!label_keep, after_stat(stratum), 
        NA_character_))), size = 3, vjust = 0.5, hjust = 0.5) + ggplot2::scale_x_discrete(limits = c(paste0("HH=", 
        sub("hap", "", hap_axes)), "Population"), expand = c(0.1, 0.05)) + ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
        ggplot2::scale_fill_manual(values = pal_all, breaks = deep_levels, drop = TRUE, na.value = "#BDBDBD") + 
        ggplot2::guides(fill = ggplot2::guide_legend(title = "Haplogroup")) + ggplot2::theme_bw(base_size = 12) + 
        ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 10), plot.margin = ggplot2::margin(10, 
            40, 10, 10)) + ggplot2::coord_cartesian(clip = "off") + ggplot2::labs(y = "Frequency", x = "", title = "Multilayer Alluvial (Frequency)")
    base <- tools::file_path_sans_ext(out_pdf)
    suppressWarnings(save_plot(paste0(base, ".pdf"), p, width = plot_base_width * 2, height = plot_base_height * 
        1.25, limitsize = FALSE))
}
plot_alluvial_multilayer_freq_singlehap_raw <- function(mdat, root_hap, out_pdf, label_min_freq = 0.01) {
    long <- mdat$long
    axes <- mdat$axes
    deepest <- mdat$deepest
    pal_all <- mdat$pal_all
    pop_levels <- mdat$pop_levels
    root_col <- axes[1]
    sub_long <- long %>% dplyr::filter(.data[[root_col]] == root_hap)
    if (nrow(sub_long) == 0) {
        message(sprintf("Root %s has no data under this HH chain; skip raw plot", root_hap))
        return(invisible(NULL))
    }
    hap_axes <- axes[-length(axes)]
    label_keep <- .compute_label_keep(sub_long, hap_axes, label_min_freq)
    mapping <- do.call(ggplot2::aes, c(list(y = quote(freq)), setNames(lapply(axes, as.name), paste0("axis", seq_along(axes)))))
    p <- ggplot2::ggplot(sub_long, mapping) + ggalluvial::geom_alluvium(ggplot2::aes(fill = .data[[deepest]]), 
        width = 1/14, knot.pos = 0.38, alpha = 0.9) + ggalluvial::geom_stratum(ggplot2::aes(fill = as.character(after_stat(stratum))), 
        width = 1/14, color = "grey30", show.legend = FALSE) + ggplot2::geom_text(stat = "stratum", ggplot2::aes(label = ifelse(as.character(after_stat(stratum)) %in% 
        !!pop_levels, as.character(after_stat(stratum)), ifelse(as.character(after_stat(stratum)) %in% !!label_keep, 
        as.character(after_stat(stratum)), NA_character_)), hjust = 0.5), size = 3, vjust = 0.5) + ggplot2::scale_x_discrete(limits = c(paste0("HH=", 
        sub("hap", "", hap_axes)), "Population"), expand = c(0.1, 0.05)) + ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
        ggplot2::scale_fill_manual(values = pal_all, breaks = unique(sub_long[[deepest]]), drop = TRUE, na.value = "#BDBDBD") + 
        ggplot2::guides(fill = ggplot2::guide_legend(title = "Haplogroup")) + ggplot2::theme_bw(base_size = 12) + 
        ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 10), plot.margin = ggplot2::margin(10, 
            40, 10, 10)) + ggplot2::coord_cartesian(clip = "off") + ggplot2::labs(y = "Frequency", x = "", title = paste0("Multilayer Alluvial (", 
        root_hap, ", raw)"))
    base <- tools::file_path_sans_ext(out_pdf)
    suppressWarnings(save_plot(paste0(base, ".pdf"), p, width = plot_base_width * 2, height = plot_base_height * 
        1.25, limitsize = FALSE))
}
plot_alluvial_multilayer_freq_singlehap_scaled <- function(mdat, root_hap, out_pdf, label_min_freq = 0.01) {
    long <- mdat$long
    axes <- mdat$axes
    deepest <- mdat$deepest
    pal_all <- mdat$pal_all
    pop_levels <- mdat$pop_levels
    root_col <- axes[1]
    sub_long <- long %>% dplyr::filter(.data[[root_col]] == root_hap)
    if (nrow(sub_long) == 0) {
        message(sprintf("Root %s has no data under this HH chain; skip scaled plot", root_hap))
        return(invisible(NULL))
    }
    sub_long <- sub_long %>% dplyr::group_by(dplyr::across(dplyr::all_of(axes))) %>% dplyr::mutate(freq_scaled = freq/sum(freq)) %>% 
        dplyr::ungroup()
    hap_axes <- axes[-length(axes)]
    sub_long_lbl <- sub_long %>% dplyr::mutate(freq = freq_scaled)
    label_keep <- .compute_label_keep(sub_long_lbl, hap_axes, label_min_freq)
    mapping <- do.call(ggplot2::aes, c(list(y = quote(freq_scaled)), setNames(lapply(axes, as.name), paste0("axis", 
        seq_along(axes)))))
    p <- ggplot2::ggplot(sub_long, mapping) + ggalluvial::geom_alluvium(ggplot2::aes(fill = .data[[deepest]]), 
        width = 1/14, knot.pos = 0.38, alpha = 0.9) + ggalluvial::geom_stratum(ggplot2::aes(fill = as.character(after_stat(stratum))), 
        width = 1/14, color = "grey30", show.legend = FALSE) + ggplot2::geom_text(stat = "stratum", ggplot2::aes(label = ifelse(as.character(after_stat(stratum)) %in% 
        !!pop_levels, as.character(after_stat(stratum)), ifelse(as.character(after_stat(stratum)) %in% !!label_keep, 
        as.character(after_stat(stratum)), NA_character_)), hjust = 0.5), size = 3, vjust = 0.5) + ggplot2::scale_x_discrete(limits = c(paste0("HH=", 
        sub("hap", "", hap_axes)), "Population"), expand = c(0.1, 0.05)) + ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
        ggplot2::scale_fill_manual(values = pal_all, breaks = unique(sub_long[[deepest]]), drop = TRUE, na.value = "#BDBDBD") + 
        ggplot2::guides(fill = ggplot2::guide_legend(title = "Haplogroup")) + ggplot2::theme_bw(base_size = 12) + 
        ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 10), plot.margin = ggplot2::margin(10, 
            40, 10, 10)) + ggplot2::coord_cartesian(clip = "off") + ggplot2::labs(y = "Scaled frequency", x = "", 
        title = paste0("Multilayer Alluvial (", root_hap, ", scaled)"))
    base <- tools::file_path_sans_ext(out_pdf)
    suppressWarnings(save_plot(paste0(base, ".pdf"), p, width = plot_base_width * 2, height = plot_base_height * 
        1.25, limitsize = FALSE))
}
plot_sankey_pop_series_withbars <- function(freq_mat_df, out_pdf, pop_order = NULL, label_threshold = 0.01, bar_halfwidth = 0.33, 
    plot_title = NULL) {
    long <- freq_mat_df %>% rename(group = 1) %>% pivot_longer(-group, names_to = "component", values_to = "value")
    if (is.null(pop_order)) 
        pop_order <- order_pop_by_dominant_hap(freq_mat_df)
    long$group <- factor(long$group, levels = pop_order)
    long <- long %>% group_by(group) %>% mutate(value = value/sum(value)) %>% ungroup()
    components <- sort(unique(long$component))
    pal <- get_palette_for_levels(components, min_n = max(256, length(components)))
    stack_df <- long %>% arrange(group, component) %>% group_by(group) %>% mutate(ymin = cumsum(lag(value, default = 0)), 
        ymax = ymin + value) %>% ungroup()
    group_pos <- data.frame(group = levels(long$group), x = seq_along(levels(long$group)), row.names = NULL)
    stack_df <- left_join(stack_df, group_pos, by = "group")
    rect_df <- stack_df %>% mutate(xmin = x - bar_halfwidth, xmax = x + bar_halfwidth, fill_col = pal[component])
    flows <- list()
    for (h in components) {
        s <- rect_df %>% filter(component == h) %>% arrange(x)
        if (nrow(s) < 2) 
            next
        for (i in 1:(nrow(s) - 1)) {
            left <- s[i, ]
            right <- s[i + 1, ]
            flows[[length(flows) + 1]] <- data.frame(component = h, x = c(left$xmax, left$xmax, right$xmin, right$xmin), 
                y = c(left$ymin, left$ymax, right$ymax, right$ymin), seg = paste0(h, "_", i), fill_col = pal[h], 
                row.names = NULL)
        }
    }
    flow_df <- if (length(flows)) 
        do.call(rbind, flows)
    else data.frame(component = character(), x = numeric(), y = numeric(), seg = character(), fill_col = character(), 
        row.names = NULL)
    label_df <- rect_df %>% filter(value >= label_threshold) %>% mutate(ymid = (ymin + ymax)/2, label_txt = component)
    p <- ggplot() + geom_polygon(data = flow_df, aes(x = x, y = y, group = seg), fill = flow_df$fill_col, alpha = 0.45, 
        color = NA) + geom_rect(data = rect_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = rect_df$fill_col, 
        color = "grey25", linewidth = 0.25, alpha = 0.98) + geom_text(data = label_df, aes(x = xmax - 0.01, y = ymid, 
        label = label_txt), hjust = 1, size = 2.8, fontface = "bold") + scale_x_continuous(breaks = group_pos$x, 
        labels = group_pos$group, position = "bottom") + scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0.02), 
        labels = scales::percent_format(accuracy = 1)) + labs(x = NULL, y = "Proportion", title = plot_title %||% 
        "Population-series haplogroup flows (with bars)") + theme_minimal(base_size = 13) + theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1), panel.background = element_rect(fill = "white", 
            color = NA), plot.background = element_rect(fill = "white", color = NA), plot.title = element_text(face = "bold", 
            size = 12, hjust = 0)) + coord_cartesian(clip = "off")
    base <- tools::file_path_sans_ext(out_pdf)
    save_plot(paste0(base, ".pdf"), p, width = nature_sankey_width, height = nature_sankey_height, limitsize = FALSE)
}
plot_sankey_pop_series_pure <- function(freq_mat_df, out_pdf, pop_order = NULL, label_threshold = 0.01, bar_halfwidth = 0.32, 
    plot_title = NULL) {
    long <- freq_mat_df %>% rename(group = 1) %>% pivot_longer(-group, names_to = "component", values_to = "value")
    if (is.null(pop_order)) 
        pop_order <- order_pop_by_dominant_hap(freq_mat_df)
    long$group <- factor(long$group, levels = pop_order)
    long <- long %>% group_by(group) %>% mutate(value = value/sum(value)) %>% ungroup()
    components <- sort(unique(long$component))
    pal <- get_palette_for_levels(components, min_n = max(256, length(components)))
    stack_df <- long %>% arrange(group, component) %>% group_by(group) %>% mutate(ymin = cumsum(lag(value, default = 0)), 
        ymax = ymin + value) %>% ungroup()
    group_pos <- data.frame(group = levels(long$group), x = seq_along(levels(long$group)), stringsAsFactors = FALSE, 
        row.names = NULL)
    stack_df <- stack_df %>% left_join(group_pos, by = "group")
    flow_polys <- list()
    label_rows <- list()
    for (comp in components) {
        sub_comp <- stack_df %>% filter(component == comp) %>% arrange(x)
        if (nrow(sub_comp) < 2) 
            next
        for (i in 1:(nrow(sub_comp) - 1)) {
            x0 <- sub_comp$x[i]
            x1 <- sub_comp$x[i + 1]
            ymin0 <- sub_comp$ymin[i]
            ymax0 <- sub_comp$ymax[i]
            ymin1 <- sub_comp$ymin[i + 1]
            ymax1 <- sub_comp$ymax[i + 1]
            v0 <- sub_comp$value[i]
            v1 <- sub_comp$value[i + 1]
            seg_id <- paste0(comp, "_", x0, "_", x1)
            poly_df <- data.frame(component = as.character(comp), x = c(x0, x0, x1, x1), y = c(ymin0, ymax0, ymax1, 
                ymin1), seg_id = seg_id, stringsAsFactors = FALSE, row.names = NULL)
            flow_polys[[length(flow_polys) + 1]] <- poly_df
            mid_x <- (x0 + x1)/2
            mid_y0 <- (ymin0 + ymax0)/2
            mid_y1 <- (ymin1 + ymax1)/2
            mid_y <- (mid_y0 + mid_y1)/2
            label_rows[[length(label_rows) + 1]] <- data.frame(component = as.character(comp), seg_id = seg_id, 
                x = mid_x, y = mid_y, v0 = v0, v1 = v1, stringsAsFactors = FALSE, row.names = NULL)
        }
    }
    flow_polys_df <- if (length(flow_polys)) 
        do.call(rbind, flow_polys)
    else data.frame(component = character(), x = numeric(), y = numeric(), seg_id = character(), row.names = NULL)
    flow_polys_df$fill_col <- pal[flow_polys_df$component]
    label_df <- if (length(label_rows)) 
        do.call(rbind, label_rows)
    else data.frame(component = character(), seg_id = character(), x = numeric(), y = numeric(), v0 = numeric(), 
        v1 = numeric(), row.names = NULL)
    if (nrow(label_df)) {
        label_df <- label_df %>% mutate(max_v = pmax(v0, v1), show = max_v >= label_threshold, label_txt = component) %>% 
            filter(show)
    }
    p <- ggplot() + geom_polygon(data = flow_polys_df, aes(x = x, y = y, group = seg_id), fill = flow_polys_df$fill_col, 
        alpha = 0.55, color = NA) + geom_text(data = label_df, aes(x = x, y = y, label = label_txt), size = 2.8, 
        color = "black", fontface = "bold", lineheight = 0.9) + scale_x_continuous(breaks = group_pos$x, labels = group_pos$group, 
        position = "bottom", expand = expansion(mult = 0.08)) + scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0.02), 
        labels = scales::percent_format(accuracy = 1)) + labs(x = NULL, y = "Proportion", title = plot_title %||% 
        "Population-series haplogroup flows (pure)") + theme_minimal(base_size = 13) + theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1), plot.margin = margin(10, 20, 10, 20), panel.background = element_rect(fill = "white", 
            color = NA), plot.background = element_rect(fill = "white", color = NA), plot.title = element_text(face = "bold", 
            hjust = 0))
    base <- tools::file_path_sans_ext(out_pdf)
    save_plot(paste0(base, ".pdf"), p, width = nature_sankey_width, height = nature_sankey_height, limitsize = FALSE)
}
plot_diversity_metric <- function(diversity_long, metric = c("shannon", "simpson"), pop_levels, pop_colors) {
    metric <- match.arg(metric)
    pop_order_by_mean <- diversity_long %>% group_by(population) %>% summarise(mv = mean(.data[[metric]], na.rm = TRUE), 
        .groups = "drop") %>% arrange(desc(mv)) %>% pull(population)
    diversity_long$population <- factor(diversity_long$population, levels = pop_order_by_mean, ordered = TRUE)
    y_lower <- min(diversity_long[[metric]], na.rm = TRUE)
    y_upper <- max(diversity_long[[metric]], na.rm = TRUE)
    p <- ggplot(diversity_long, aes(x = population, y = .data[[metric]])) + ggbeeswarm::geom_quasirandom(aes(color = population), 
        method = qr_method, width = qr_width, alpha = qr_alpha, size = qr_size, show.legend = FALSE) + geom_boxplot(aes(color = population), 
        width = 0.5, outlier.shape = NA, fill = NA, size = 0.6, show.legend = FALSE) + stat_summary(fun = median, 
        geom = "crossbar", width = 0.5, fatten = 0, color = "red", size = 0.6) + scale_color_manual(values = pop_colors[pop_order_by_mean]) + 
        labs(x = "Population", y = metric, title = paste0("Population diversity (", metric, ")")) + theme_bw() + 
        theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 8), axis.text.y = element_text(size = 8), 
            axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9), plot.title = element_text(hjust = 0.5, 
                size = 10)) + coord_cartesian(ylim = c(y_lower, y_upper))
    p
}
