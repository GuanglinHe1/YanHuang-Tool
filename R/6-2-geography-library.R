auto_legend_ncol <- function(n_items, max_rows = sankey_legend_max_rows, max_cols = sankey_legend_max_cols) {
    if (!is.finite(n_items) || n_items <= 0) 
        return(1L)
    ncol <- as.integer(ceiling(n_items/max_rows))
    ncol <- max(1L, min(ncol, as.integer(max_cols)))
    ncol
}
scale_color_rb <- function(..., limits = NULL, name = waiver()) {
    if (is.null(limits) || length(limits) != 2 || any(!is.finite(limits))) {
        scale_color_gradient(low = "blue", high = "red", name = name, oob = scales::squish, ...)
    }
    else {
        scale_color_gradient(low = "blue", high = "red", limits = limits, name = name, oob = scales::squish, ...)
    }
}
scale_fill_rb <- function(..., limits = NULL, name = waiver()) {
    if (is.null(limits) || length(limits) != 2 || any(!is.finite(limits))) {
        scale_fill_gradient(low = "blue", high = "red", name = name, oob = scales::squish, ...)
    }
    else {
        scale_fill_gradient(low = "blue", high = "red", limits = limits, name = name, oob = scales::squish, ...)
    }
}
apply_pop_order <- function(pops, pop_order_from_file) {
    pops <- unique(as.character(pops))
    if (is.null(pop_order_from_file) || length(pop_order_from_file) == 0) 
        return(pops)
    in_file <- pop_order_from_file[pop_order_from_file %in% pops]
    rest <- setdiff(pops, in_file)
    c(in_file, rest)
}
bezier_path <- function(x0, y0, x1, y1, x2, y2, n = 16) {
    t <- seq(0, 1, length.out = n)
    bx <- (1 - t)^2 * x0 + 2 * (1 - t) * t * x1 + t^2 * x2
    by <- (1 - t)^2 * y0 + 2 * (1 - t) * t * y1 + t^2 * y2
    data.frame(x = bx, y = by, row.names = NULL)
}
adjust_bar_positions <- function(pop_coords, bar_halfwidth_deg = default_bar_halfwidth_deg, bar_height_deg = default_bar_height_deg, 
    max_iter = 200, eps = 1e-04) {
    n <- nrow(pop_coords)
    if (n <= 1) 
        return(pop_coords)
    centers <- as.matrix(pop_coords[, c("longitude", "latitude")])
    width <- 2 * bar_halfwidth_deg
    height <- bar_height_deg
    min_dx <- width * 1.05
    min_dy <- height * 1.05
    for (iter in seq_len(max_iter)) {
        moved <- FALSE
        for (i in 1:(n - 1)) {
            for (j in (i + 1):n) {
                dx <- centers[j, 1] - centers[i, 1]
                dy <- centers[j, 2] - centers[i, 2]
                overlap_x <- abs(dx) < min_dx
                overlap_y <- abs(dy) < min_dy
                if (overlap_x && overlap_y) {
                  need_x <- min_dx - abs(dx)
                  need_y <- min_dy - abs(dy)
                  if (need_x <= need_y) {
                    dir <- ifelse(dx >= 0, 1, -1)
                    shift <- need_x/2 + eps
                    centers[i, 1] <- centers[i, 1] - dir * shift
                    centers[j, 1] <- centers[j, 1] + dir * shift
                  }
                  else {
                    dir <- ifelse(dy >= 0, 1, -1)
                    shift <- need_y/2 + eps
                    centers[i, 2] <- centers[i, 2] - dir * shift
                    centers[j, 2] <- centers[j, 2] + dir * shift
                  }
                  moved <- TRUE
                }
            }
        }
        if (!moved) 
            break
    }
    pop_coords$longitude <- centers[, 1]
    pop_coords$latitude <- centers[, 2]
    pop_coords
}
separate_close_longitudes <- function(pop_coords, min_diff = 1, frac_shift = 0.05, shift_deg = NA_real_, passes = 2) {
    if (nrow(pop_coords) < 2) 
        return(pop_coords)
    lon_range <- range(pop_coords$longitude, na.rm = TRUE)
    total_span <- lon_range[2] - lon_range[1]
    if (!is.finite(total_span) || total_span <= 0) 
        total_span <- 1
    shift_each <- if (is.finite(shift_deg)) 
        shift_deg
    else frac_shift * total_span
    if (!is.finite(shift_each) || shift_each <= 0) 
        shift_each <- 0.1
    for (pp in seq_len(passes)) {
        moved <- FALSE
        for (i in 1:(nrow(pop_coords) - 1)) {
            for (j in (i + 1):nrow(pop_coords)) {
                dlon <- abs(pop_coords$longitude[i] - pop_coords$longitude[j])
                if (is.finite(dlon) && dlon < min_diff) {
                  mid <- (pop_coords$longitude[i] + pop_coords$longitude[j])/2
                  pop_coords$longitude[i] <- mid - shift_each
                  pop_coords$longitude[j] <- mid + shift_each
                  moved <- TRUE
                }
            }
        }
        if (!moved) 
            break
    }
    pop_coords
}
max_geo_dist_for_pops <- function(pop_coords, pops) {
    coords <- pop_coords %>% dplyr::filter(population %in% pops)
    if (nrow(coords) < 2) 
        return(1)
    mat <- as.matrix(coords[, c("longitude", "latitude")])
    comb <- utils::combn(nrow(mat), 2)
    d <- sqrt(colSums((mat[comb[1, ], , drop = FALSE] - mat[comb[2, ], , drop = FALSE])^2))
    max(d, na.rm = TRUE)
}
plot_geo_sankey_allpairs <- function(freq_mat_df, pop_coords, out_pdf, bar_halfwidth_deg = default_bar_halfwidth_deg, 
    bar_height_deg = default_bar_height_deg, plot_title = NULL, flow_min_freq = 0.03, curve_factor = 0.35, bezier_points = 16, 
    component_min_freq = sankey_component_min_freq, label_min_freq = sankey_label_min_freq) {
    long <- freq_mat_df %>% dplyr::rename(population = 1) %>% tidyr::pivot_longer(-population, names_to = "component", 
        values_to = "value") %>% dplyr::left_join(pop_coords, by = "population") %>% dplyr::filter(!is.na(longitude), 
        !is.na(latitude))
    if (nrow(long) == 0) {
        message2("No populations with coordinates for geo-Sankey; skip all-pairs.")
        return(invisible(NULL))
    }
    long <- long %>% dplyr::group_by(population) %>% dplyr::mutate(value = value/sum(value)) %>% dplyr::mutate(keep = ((value >= 
        component_min_freq) & (value > 0)) | (value == max(value, na.rm = TRUE))) %>% dplyr::filter(keep) %>% dplyr::mutate(value = value/sum(value)) %>% 
        dplyr::ungroup()
    components_all <- sort(unique(long$component))
    pal <- get_palette_for_levels(components_all, min_n = max(256, length(components_all)))
    comp_stats <- long %>% dplyr::group_by(component) %>% dplyr::summarise(maxp = max(value, na.rm = TRUE), .groups = "drop")
    legend_components <- comp_stats %>% dplyr::filter(is.finite(maxp) & maxp > sankey_legend_min_show_freq) %>% 
        dplyr::pull(component) %>% sort()
    if (length(legend_components) == 0) 
        legend_components <- components_all[1]
    legend_ncol <- auto_legend_ncol(length(legend_components))
    stack_df <- long %>% dplyr::arrange(population, component) %>% dplyr::group_by(population) %>% dplyr::mutate(ymin_prop = cumsum(dplyr::lag(value, 
        default = 0)), ymax_prop = ymin_prop + value) %>% dplyr::ungroup() %>% dplyr::mutate(bar_ymin = latitude - 
        bar_height_deg/2, ymin = bar_ymin + ymin_prop * bar_height_deg, ymax = bar_ymin + ymax_prop * bar_height_deg, 
        xmin = longitude - bar_halfwidth_deg, xmax = longitude + bar_halfwidth_deg)
    label_pad <- max(0.02, bar_halfwidth_deg * sankey_label_left_pad_frac)
    min_h_deg <- bar_height_deg * sankey_label_min_height_prop
    label_df <- stack_df %>% dplyr::filter(value >= label_min_freq) %>% dplyr::filter((ymax - ymin) >= min_h_deg) %>% 
        dplyr::mutate(x_label = xmin + label_pad, y_label = (ymin + ymax)/2)
    pops <- unique(stack_df$population)
    max_dist <- max_geo_dist_for_pops(pop_coords, pops)
    if (!is.finite(max_dist) || max_dist <= 0) 
        max_dist <- 1
    flows <- list()
    for (h in components_all) {
        s <- stack_df %>% dplyr::filter(component == h)
        if (nrow(s) < 2) 
            next
        s <- s %>% dplyr::arrange(longitude)
        for (i in 1:(nrow(s) - 1)) {
            for (j in (i + 1):nrow(s)) {
                a <- s[i, ]
                b <- s[j, ]
                max_v <- max(a$value, b$value)
                if (max_v < flow_min_freq) 
                  next
                if (b$longitude >= a$longitude) {
                  left <- a
                  right <- b
                }
                else {
                  left <- b
                  right <- a
                }
                x_start <- left$xmax
                y_top_s <- left$ymax
                y_bot_s <- left$ymin
                x_end <- right$xmin
                y_top_e <- right$ymax
                y_bot_e <- right$ymin
                dx <- right$longitude - left$longitude
                dy <- right$latitude - left$latitude
                dist_ij <- sqrt(dx^2 + dy^2)
                w_norm <- (max_dist - dist_ij)/max_dist
                if (!is.finite(w_norm)) 
                  w_norm <- 0
                w_norm <- max(0, min(1, w_norm))
                center_left <- (left$ymin + left$ymax)/2
                center_right <- (right$ymin + right$ymax)/2
                sign_dir <- ifelse(center_right >= center_left, 1, -1)
                amp <- curve_factor * bar_height_deg * (1 + w_norm) * sign_dir
                mid_x <- (x_start + x_end)/2
                mid_top_y <- (y_top_s + y_top_e)/2 + amp
                mid_bot_y <- (y_bot_s + y_bot_e)/2 + amp
                top_path <- bezier_path(x_start, y_top_s, mid_x, mid_top_y, x_end, y_top_e, n = bezier_points)
                bottom_path <- bezier_path(x_end, y_bot_e, mid_x, mid_bot_y, x_start, y_bot_s, n = bezier_points)
                seg_id <- paste0(h, "_", left$population, "_", right$population)
                flows[[length(flows) + 1]] <- rbind(data.frame(component = h, x = top_path$x, y = top_path$y, seg = seg_id, 
                  stringsAsFactors = FALSE, row.names = NULL), data.frame(component = h, x = bottom_path$x, y = bottom_path$y, 
                  seg = seg_id, stringsAsFactors = FALSE, row.names = NULL))
            }
        }
    }
    flow_df <- if (length(flows)) 
        do.call(rbind, flows)
    else data.frame(component = character(), x = numeric(), y = numeric(), seg = character(), stringsAsFactors = FALSE, 
        row.names = NULL)
    pop_label_df <- pop_coords %>% dplyr::filter(population %in% pops) %>% dplyr::mutate(y_label = latitude - bar_height_deg/2 - 
        0.3)
    world <- ggplot2::map_data("world")
    lon_range <- range(pop_coords$longitude, na.rm = TRUE)
    lat_range <- range(pop_coords$latitude, na.rm = TRUE)
    xlim <- c(lon_range[1] - map_x_pad, lon_range[2] + map_x_pad)
    ylim <- c(lat_range[1] - bar_height_deg/2 - map_y_pad, lat_range[2] + bar_height_deg/2 + map_y_pad)
    asp <- diff(xlim)/diff(ylim)
    if (!is.finite(asp) || asp <= 0) 
        asp <- 1
    height_in <- map_base_height_in
    width_in <- height_in * asp * (1.1 + 0.35 * (legend_ncol - 1))
    p <- ggplot() + geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "grey95", color = "grey80", 
        linewidth = 0.2) + geom_polygon(data = flow_df, aes(x = x, y = y, group = seg, fill = component), alpha = 0.25, 
        color = NA) + geom_rect(data = stack_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = component), 
        color = "grey25", linewidth = 0.25, alpha = 0.98) + geom_text(data = label_df, aes(x = x_label, y = y_label, 
        label = component), size = 2.8, color = "black", hjust = 0, vjust = 0.5, check_overlap = TRUE) + geom_text(data = pop_label_df, 
        aes(x = longitude, y = y_label, label = population), size = 3, vjust = 1) + coord_quickmap(xlim = xlim, 
        ylim = ylim, expand = FALSE) + scale_fill_manual(values = pal, limits = components_all, breaks = legend_components, 
        guide = guide_legend(title = "Haplogroup", ncol = legend_ncol, byrow = TRUE, override.aes = list(alpha = 1))) + 
        labs(x = "Longitude", y = "Latitude", title = plot_title %||% sprintf("Geo-Sankey (all-pairs, HH) (freq>=%.2f)", 
            component_min_freq)) + theme_minimal(base_size = 12) + theme(panel.grid = element_line(color = "grey90", 
        linewidth = 0.3), legend.position = "right", legend.title = element_text(face = "bold"), legend.text = element_text(size = 9), 
        legend.key.height = grid::unit(9, "pt"), legend.key.width = grid::unit(10, "pt"), plot.title = element_text(face = "bold", 
            hjust = 0))
    base <- tools::file_path_sans_ext(out_pdf)
    save_plot(paste0(base, ".pdf"), p, width = width_in, height = height_in, limitsize = FALSE)
}
plot_geo_sankey_star <- function(freq_mat_df, pop_coords, target_population, out_pdf, bar_halfwidth_deg = default_bar_halfwidth_deg, 
    bar_height_deg = default_bar_height_deg, plot_title = NULL, flow_min_freq = 0.03, curve_factor = 0.35, bezier_points = 16, 
    component_min_freq = sankey_component_min_freq, label_min_freq = sankey_label_min_freq) {
    long <- freq_mat_df %>% dplyr::rename(population = 1) %>% tidyr::pivot_longer(-population, names_to = "component", 
        values_to = "value") %>% dplyr::left_join(pop_coords, by = "population") %>% dplyr::filter(!is.na(longitude), 
        !is.na(latitude))
    if (nrow(long) == 0) {
        message2("No coordinates; skip star.")
        return(invisible(NULL))
    }
    if (!target_population %in% long$population) {
        message2("Target '%s' not in data; skip.", target_population)
        return(invisible(NULL))
    }
    long <- long %>% dplyr::group_by(population) %>% dplyr::mutate(value = value/sum(value)) %>% dplyr::mutate(keep = ((value >= 
        component_min_freq) & (value > 0)) | (value == max(value, na.rm = TRUE))) %>% dplyr::filter(keep) %>% dplyr::mutate(value = value/sum(value)) %>% 
        dplyr::ungroup()
    components_all <- sort(unique(long$component))
    pal <- get_palette_for_levels(components_all, min_n = max(256, length(components_all)))
    comp_stats <- long %>% dplyr::group_by(component) %>% dplyr::summarise(maxp = max(value, na.rm = TRUE), .groups = "drop")
    legend_components <- comp_stats %>% dplyr::filter(is.finite(maxp) & maxp > sankey_legend_min_show_freq) %>% 
        dplyr::pull(component) %>% sort()
    if (length(legend_components) == 0) 
        legend_components <- components_all[1]
    legend_ncol <- auto_legend_ncol(length(legend_components))
    stack_df <- long %>% dplyr::arrange(population, component) %>% dplyr::group_by(population) %>% dplyr::mutate(ymin_prop = cumsum(dplyr::lag(value, 
        default = 0)), ymax_prop = ymin_prop + value) %>% dplyr::ungroup() %>% dplyr::mutate(bar_ymin = latitude - 
        bar_height_deg/2, ymin = bar_ymin + ymin_prop * bar_height_deg, ymax = bar_ymin + ymax_prop * bar_height_deg, 
        xmin = longitude - bar_halfwidth_deg, xmax = longitude + bar_halfwidth_deg)
    label_pad <- max(0.02, bar_halfwidth_deg * sankey_label_left_pad_frac)
    min_h_deg <- bar_height_deg * sankey_label_min_height_prop
    label_df <- stack_df %>% dplyr::filter(value >= label_min_freq) %>% dplyr::filter((ymax - ymin) >= min_h_deg) %>% 
        dplyr::mutate(x_label = xmin + label_pad, y_label = (ymin + ymax)/2)
    pops <- unique(stack_df$population)
    max_dist <- max_geo_dist_for_pops(pop_coords, pops)
    if (!is.finite(max_dist) || max_dist <= 0) 
        max_dist <- 1
    flows <- list()
    for (h in components_all) {
        t_row <- stack_df %>% dplyr::filter(population == target_population, component == h)
        if (nrow(t_row) == 0) 
            next
        t_row <- t_row[1, ]
        s_rows <- stack_df %>% dplyr::filter(component == h, population != target_population)
        if (nrow(s_rows) == 0) 
            next
        for (k in seq_len(nrow(s_rows))) {
            s <- s_rows[k, ]
            max_v <- max(s$value, t_row$value)
            if (max_v < flow_min_freq) 
                next
            if (t_row$longitude >= s$longitude) {
                x_start <- s$xmax
                x_end <- t_row$xmin
            }
            else {
                x_start <- s$xmin
                x_end <- t_row$xmax
            }
            y_top_s <- s$ymax
            y_bot_s <- s$ymin
            y_top_e <- t_row$ymax
            y_bot_e <- t_row$ymin
            dx <- t_row$longitude - s$longitude
            dy <- t_row$latitude - s$latitude
            dist_ij <- sqrt(dx^2 + dy^2)
            w_norm <- (max_dist - dist_ij)/max_dist
            if (!is.finite(w_norm)) 
                w_norm <- 0
            w_norm <- max(0, min(1, w_norm))
            center_s <- (s$ymin + s$ymax)/2
            center_t <- (t_row$ymin + t_row$ymax)/2
            sign_dir <- ifelse(center_t >= center_s, 1, -1)
            amp <- curve_factor * bar_height_deg * (1 + w_norm) * sign_dir
            mid_x <- (x_start + x_end)/2
            mid_top_y <- (y_top_s + y_top_e)/2 + amp
            mid_bot_y <- (y_bot_s + y_bot_e)/2 + amp
            top_path <- bezier_path(x_start, y_top_s, mid_x, mid_top_y, x_end, y_top_e, n = bezier_points)
            bottom_path <- bezier_path(x_end, y_bot_e, mid_x, mid_bot_y, x_start, y_bot_s, n = bezier_points)
            seg_id <- paste0(h, "_", s$population, "_to_", target_population)
            flows[[length(flows) + 1]] <- rbind(data.frame(component = h, x = top_path$x, y = top_path$y, seg = seg_id, 
                stringsAsFactors = FALSE, row.names = NULL), data.frame(component = h, x = bottom_path$x, y = bottom_path$y, 
                seg = seg_id, stringsAsFactors = FALSE, row.names = NULL))
        }
    }
    flow_df <- if (length(flows)) 
        do.call(rbind, flows)
    else data.frame(component = character(), x = numeric(), y = numeric(), seg = character(), stringsAsFactors = FALSE, 
        row.names = NULL)
    pop_label_df <- pop_coords %>% dplyr::filter(population %in% pops) %>% dplyr::mutate(y_label = latitude - bar_height_deg/2 - 
        0.3)
    world <- ggplot2::map_data("world")
    lon_range <- range(pop_coords$longitude, na.rm = TRUE)
    lat_range <- range(pop_coords$latitude, na.rm = TRUE)
    xlim <- c(lon_range[1] - map_x_pad, lon_range[2] + map_x_pad)
    ylim <- c(lat_range[1] - bar_height_deg/2 - map_y_pad, lat_range[2] + bar_height_deg/2 + map_y_pad)
    asp <- diff(xlim)/diff(ylim)
    if (!is.finite(asp) || asp <= 0) 
        asp <- 1
    height_in <- map_base_height_in
    width_in <- height_in * asp * (1.1 + 0.35 * (legend_ncol - 1))
    p <- ggplot() + geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "grey95", color = "grey80", 
        linewidth = 0.2) + geom_polygon(data = flow_df, aes(x = x, y = y, group = seg, fill = component), alpha = 0.25, 
        color = NA) + geom_rect(data = stack_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = component), 
        color = "grey25", linewidth = 0.25, alpha = 0.98) + geom_text(data = label_df, aes(x = x_label, y = y_label, 
        label = component), size = 2.8, color = "black", hjust = 0, vjust = 0.5, check_overlap = TRUE) + geom_text(data = pop_label_df, 
        aes(x = longitude, y = y_label, label = population), size = 3, vjust = 1) + coord_quickmap(xlim = xlim, 
        ylim = ylim, expand = FALSE) + scale_fill_manual(values = pal, limits = components_all, breaks = legend_components, 
        guide = guide_legend(title = "Haplogroup", ncol = legend_ncol, byrow = TRUE, override.aes = list(alpha = 1))) + 
        labs(x = "Longitude", y = "Latitude", title = plot_title %||% sprintf("Geo-Sankey (star, target=%s) (freq>=%.2f)", 
            target_population, component_min_freq)) + theme_minimal(base_size = 12) + theme(panel.grid = element_line(color = "grey90", 
        linewidth = 0.3), legend.position = "right", legend.title = element_text(face = "bold"), legend.text = element_text(size = 9), 
        legend.key.height = grid::unit(9, "pt"), legend.key.width = grid::unit(10, "pt"), plot.title = element_text(face = "bold", 
            hjust = 0))
    base <- tools::file_path_sans_ext(out_pdf)
    save_plot(paste0(base, ".pdf"), p, width = width_in, height = height_in, limitsize = FALSE)
}
as_square_matrix <- function(x, expected_names = NULL) {
    if (is.data.frame(x)) {
        df <- x
        if (ncol(df) >= 2 && is.character(df[[1]]) && all(df[[1]] %in% colnames(df))) {
            rn <- df[[1]]
            df <- df[, -1, drop = FALSE]
            rownames(df) <- rn
        }
        mat_chr <- as.matrix(df)
        mat_num <- suppressWarnings(apply(mat_chr, 2, as.numeric))
        if (is.null(dim(mat_num))) 
            mat_num <- matrix(mat_num, nrow = nrow(mat_chr))
        rownames(mat_num) <- rownames(df)
        colnames(mat_num) <- colnames(df)
        x <- mat_num
    }
    if (is.matrix(x) && nrow(x) > 1 && ncol(x) == 1) 
        x <- as.vector(x[, 1])
    if (is.matrix(x) && nrow(x) == 1 && ncol(x) > 1) 
        x <- as.vector(x[1, ])
    if (is.vector(x) && !is.list(x)) {
        L <- length(x)
        n <- sqrt(L)
        if (!is.finite(n) || abs(n - round(n)) > 1e-08) 
            stop("as_square_matrix: vector length not perfect square: ", L)
        n <- as.integer(round(n))
        mat <- matrix(suppressWarnings(as.numeric(x)), nrow = n, ncol = n, byrow = FALSE)
        if (!is.null(expected_names) && length(expected_names) == n) {
            rownames(mat) <- expected_names
            colnames(mat) <- expected_names
        }
        return(mat)
    }
    mat <- as.matrix(x)
    storage.mode(mat) <- "numeric"
    if (!is.null(expected_names) && length(expected_names) > 0 && !is.null(rownames(mat)) && !is.null(colnames(mat))) {
        common <- expected_names[expected_names %in% rownames(mat)]
        common <- common[common %in% colnames(mat)]
        mat <- mat[common, common, drop = FALSE]
    }
    else if (!is.null(rownames(mat)) && !is.null(colnames(mat))) {
        common <- intersect(rownames(mat), colnames(mat))
        mat <- mat[common, common, drop = FALSE]
    }
    if (nrow(mat) != ncol(mat)) 
        stop("as_square_matrix: not square; dims=", nrow(mat), "x", ncol(mat))
    mat
}
symmetrize_matrix <- function(mat, diag_value = 1, expected_names = NULL) {
    mat <- as_square_matrix(mat, expected_names = expected_names)
    mat[!is.finite(mat)] <- NA_real_
    mat <- (mat + t(mat))/2
    diag(mat) <- diag_value
    mat
}
clamp01 <- function(x) {
    x[!is.finite(x)] <- 0
    x[x < 0] <- 0
    x[x > 1] <- 1
    x
}
row_normalize <- function(M) {
    M <- as.matrix(M)
    rs <- rowSums(M)
    rs[!is.finite(rs) | rs == 0] <- 1
    M/rs
}
cosine_sim_matrix <- function(M) {
    M <- as.matrix(M)
    num <- tcrossprod(M)
    norms <- sqrt(rowSums(M^2))
    denom <- outer(norms, norms, "*")
    S <- num/denom
    S[!is.finite(S)] <- 0
    diag(S) <- 1
    S
}
js_similarity_matrix <- function(M, eps = 1e-12) {
    M <- as.matrix(M)
    n <- nrow(M)
    out <- matrix(0, n, n)
    for (i in 1:n) {
        out[i, i] <- 1
        if (i < n) {
            for (j in (i + 1):n) {
                d <- JSD_pair(M[i, ], M[j, ], eps = eps)^2
                d01 <- sqrt(d/log(2))
                s <- 1 - d01
                if (!is.finite(s)) 
                  s <- 0
                out[i, j] <- s
                out[j, i] <- s
            }
        }
    }
    out
}
within_pop_consistency <- function(M, eps = 1e-12) {
    M <- row_normalize(as.matrix(M) + eps)
    K <- ncol(M)
    HHI <- rowSums(M^2)
    dom <- apply(M, 1, max)
    eff <- 1/HHI
    if (K <= 1) {
        ent <- rep(0, nrow(M))
        ent_norm <- rep(0, nrow(M))
        cons_1mE <- rep(1, nrow(M))
    }
    else {
        ent <- -rowSums(M * log(M))
        ent_norm <- ent/log(K)
        cons_1mE <- 1 - ent_norm
    }
    data.frame(HHI = HHI, Consistency_1minusEntropy = cons_1mE, Dominance_maxp = dom, Effective_num = eff, Shannon = ent, 
        Shannon_norm = ent_norm, stringsAsFactors = FALSE, row.names = NULL)
}
haversine_km <- function(lon1, lat1, lon2, lat2) {
    R <- 6371.0088
    toRad <- pi/180
    lon1 <- lon1 * toRad
    lat1 <- lat1 * toRad
    lon2 <- lon2 * toRad
    lat2 <- lat2 * toRad
    dlon <- lon2 - lon1
    dlat <- lat2 - lat1
    a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
    sqrt_a <- sqrt(a)
    if (!is.finite(sqrt_a)) 
        sqrt_a <- 0
    if (sqrt_a > 1) 
        sqrt_a <- 1
    c <- 2 * asin(sqrt_a)
    R * c
}
geo_dist_matrix_km <- function(lon_deg, lat_deg) {
    n <- length(lon_deg)
    out <- matrix(0, n, n)
    for (i in 1:n) {
        if (i < n) {
            for (j in (i + 1):n) {
                d <- haversine_km(lon_deg[i], lat_deg[i], lon_deg[j], lat_deg[j])
                out[i, j] <- d
                out[j, i] <- d
            }
        }
    }
    out
}
knn_mean_similarity <- function(sim_mat, geo_dist_mat, k = 5) {
    n <- nrow(sim_mat)
    if (n <= 1) 
        return(rep(NA_real_, n))
    k <- min(k, n - 1)
    out <- rep(NA_real_, n)
    for (i in 1:n) {
        ord <- order(geo_dist_mat[i, ], na.last = NA)
        ord <- ord[ord != i]
        nn <- head(ord, k)
        out[i] <- mean(sim_mat[i, nn], na.rm = TRUE)
    }
    out
}
plot_heatmap_matrix_clustered <- function(mat, out_prefix, title, is_similarity = TRUE, cluster_method = "ward.D2", 
    show_dendrogram = TRUE, expected_names = NULL) {
    mat <- as_square_matrix(mat, expected_names = expected_names)
    if (nrow(mat) < 2) {
        message2("Heatmap skipped (n<2): %s", title)
        return(invisible(NULL))
    }
    if (is_similarity) {
        mat <- symmetrize_matrix(mat, diag_value = 1, expected_names = expected_names)
        mat <- clamp01(mat)
        diag(mat) <- 1
        D <- 1 - mat
        D[!is.finite(D)] <- 1
        D[D < 0] <- 0
        D[D > 1] <- 1
        diag(D) <- 0
        d <- stats::as.dist(D)
    }
    else {
        mat <- symmetrize_matrix(mat, diag_value = 0, expected_names = expected_names)
        mat[!is.finite(mat)] <- 0
        mat[mat < 0] <- 0
        diag(mat) <- 0
        d <- stats::as.dist(mat)
    }
    hc <- stats::hclust(d, method = cluster_method)
    labs <- rownames(mat)[hc$order]
    mat2 <- mat[labs, labs, drop = FALSE]
    if (show_dendrogram) {
        open_pdf(paste0(out_prefix, "_dendrogram.pdf"), width = 7, height = 5)
        plot(hc, main = paste0(title, " (", cluster_method, ")"), xlab = "", sub = "")
        dev.off()
    }
    df <- as.data.frame(as.table(mat2))
    colnames(df) <- c("Pop1", "Pop2", "Value")
    df$Pop1 <- factor(df$Pop1, levels = labs)
    df$Pop2 <- factor(df$Pop2, levels = rev(labs))
    show_text <- length(labs) <= max_heatmap_labels
    p <- ggplot(df, aes(x = Pop1, y = Pop2, fill = Value)) + geom_tile() + scale_fill_rb(limits = c(0, 1), name = "Value") + 
        labs(title = title, x = NULL, y = NULL) + theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold", 
        hjust = 0), panel.grid = element_blank())
    if (show_text) {
        p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    }
    else {
        p <- p + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())
    }
    save_plot(paste0(out_prefix, ".pdf"), p, width = 8.5, height = 7.5, limitsize = FALSE)
}
plot_metric_map <- function(pop_df, metric_col, out_prefix, title) {
    world <- ggplot2::map_data("world")
    lon_range <- range(pop_df$longitude, na.rm = TRUE)
    lat_range <- range(pop_df$latitude, na.rm = TRUE)
    xlim <- c(lon_range[1] - map_x_pad, lon_range[2] + map_x_pad)
    ylim <- c(lat_range[1] - map_y_pad, lat_range[2] + map_y_pad)
    asp <- diff(xlim)/diff(ylim)
    if (!is.finite(asp) || asp <= 0) 
        asp <- 1
    height_in <- map_base_height_in
    width_in <- height_in * asp
    show_labels <- nrow(pop_df) <= label_pop_if_nle
    p <- ggplot() + geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "grey95", color = "grey80", 
        linewidth = 0.2) + geom_point(data = pop_df, aes(x = longitude, y = latitude, color = .data[[metric_col]], 
        size = n_samples), alpha = 0.95) + scale_color_rb(name = metric_col) + scale_size_continuous(range = c(1.8, 
        5.2)) + coord_quickmap(xlim = xlim, ylim = ylim, expand = FALSE) + labs(title = title, x = "Longitude", 
        y = "Latitude", size = "N") + theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold", 
        hjust = 0), panel.grid = element_line(color = "grey90", linewidth = 0.3))
    if (show_labels) {
        p <- p + geom_text(data = pop_df, aes(x = longitude, y = latitude, label = population), size = 2.8, vjust = -0.8)
    }
    save_plot(paste0(out_prefix, ".pdf"), p, width = width_in, height = height_in, limitsize = FALSE)
}
make_edge_table <- function(pop_df, sim_mat, edge_mode = c("topn", "threshold"), top_n = 1, min_sim = 0.7, max_edges = 1200, 
    cons_col = "Consistency_1minusEntropy") {
    edge_mode <- match.arg(edge_mode)
    pop_df <- pop_df %>% dplyr::distinct(population, .keep_all = TRUE)
    pops <- pop_df$population
    sim_mat <- as_square_matrix(sim_mat, expected_names = pops)
    sim_mat <- symmetrize_matrix(sim_mat, diag_value = 1, expected_names = pops)
    sim_mat <- clamp01(sim_mat)
    if (!(cons_col %in% names(pop_df))) 
        stop("cons_col not found: ", cons_col)
    cons <- pop_df[[cons_col]]
    names(cons) <- pop_df$population
    edges <- list()
    n <- length(pops)
    if (edge_mode == "topn") {
        for (i in 1:n) {
            s <- sim_mat[i, ]
            s[i] <- -Inf
            ord <- order(s, decreasing = TRUE)
            pick <- head(ord, min(top_n, sum(is.finite(s) & s > -Inf)))
            for (j in pick) {
                if (!is.finite(s[j])) 
                  next
                pi <- pops[i]
                pj <- pops[j]
                ci <- cons[pi]
                cj <- cons[pj]
                if (!is.finite(ci)) 
                  ci <- 0
                if (!is.finite(cj)) 
                  cj <- 0
                if (ci < 0) 
                  ci <- 0
                if (cj < 0) 
                  cj <- 0
                pc <- sqrt(ci * cj)
                edges[[length(edges) + 1]] <- data.frame(from = pi, to = pj, sim = as.numeric(s[j]), pair_cons = as.numeric(pc), 
                  stringsAsFactors = FALSE, row.names = NULL)
            }
        }
    }
    else {
        for (i in 1:n) {
            for (j in 1:n) {
                if (i == j) 
                  next
                s <- sim_mat[i, j]
                if (!is.finite(s) || s < min_sim) 
                  next
                pi <- pops[i]
                pj <- pops[j]
                ci <- cons[pi]
                cj <- cons[pj]
                if (!is.finite(ci)) 
                  ci <- 0
                if (!is.finite(cj)) 
                  cj <- 0
                if (ci < 0) 
                  ci <- 0
                if (cj < 0) 
                  cj <- 0
                pc <- sqrt(ci * cj)
                edges[[length(edges) + 1]] <- data.frame(from = pi, to = pj, sim = as.numeric(s), pair_cons = as.numeric(pc), 
                  stringsAsFactors = FALSE, row.names = NULL)
            }
        }
    }
    ed <- if (length(edges)) 
        do.call(rbind, edges)
    else data.frame(from = character(), to = character(), sim = numeric(), pair_cons = numeric(), stringsAsFactors = FALSE, 
        row.names = NULL)
    if (nrow(ed) == 0) 
        return(ed)
    ed <- ed %>% dplyr::filter(from %in% pops, to %in% pops)
    ed <- ed %>% dplyr::arrange(dplyr::desc(sim), dplyr::desc(pair_cons))
    if (is.finite(max_edges) && nrow(ed) > max_edges) 
        ed <- ed[1:max_edges, , drop = FALSE]
    ed
}
empty_path_df <- function() {
    data.frame(edge_id = character(), from = character(), to = character(), x = numeric(), y = numeric(), sim = numeric(), 
        pair_cons = numeric(), dist_deg = numeric(), stringsAsFactors = FALSE, row.names = NULL)
}
bezier_edge_points <- function(ed, pop_df, curv_base = curv_base, curv_alpha = curv_alpha, curv_clip = curv_clip, 
    n_points = bezier_points) {
    if (nrow(ed) == 0) 
        return(empty_path_df())
    ed2 <- ed %>% dplyr::left_join(pop_df %>% dplyr::select(population, longitude, latitude), by = c(from = "population")) %>% 
        dplyr::rename(x0 = longitude, y0 = latitude) %>% dplyr::left_join(pop_df %>% dplyr::select(population, 
        longitude, latitude), by = c(to = "population")) %>% dplyr::rename(x2 = longitude, y2 = latitude) %>% dplyr::filter(is.finite(x0), 
        is.finite(y0), is.finite(x2), is.finite(y2))
    if (nrow(ed2) == 0) 
        return(empty_path_df())
    ed2$dist_deg <- sqrt((ed2$x2 - ed2$x0)^2 + (ed2$y2 - ed2$y0)^2)
    max_deg <- max(ed2$dist_deg, na.rm = TRUE)
    if (!is.finite(max_deg) || max_deg <= 0) 
        max_deg <- 1
    ed2$dist_norm <- ed2$dist_deg/max_deg
    ed2$theta <- atan2(ed2$y2 - ed2$y0, ed2$x2 - ed2$x0)
    ed2$curv_raw <- curv_base * (ed2$dist_norm^curv_alpha) * sin(ed2$theta)
    ed2$curv <- ed2$curv_raw
    ed2$curv[!is.finite(ed2$curv)] <- 0
    ed2$curv[ed2$curv > curv_clip] <- curv_clip
    ed2$curv[ed2$curv < -curv_clip] <- -curv_clip
    max_amp <- 10
    ed2$amp <- ed2$curv * ed2$dist_deg
    ed2$amp[!is.finite(ed2$amp)] <- 0
    ed2$amp[ed2$amp > max_amp] <- max_amp
    ed2$amp[ed2$amp < -max_amp] <- -max_amp
    res_list <- vector("list", nrow(ed2))
    t <- seq(0, 1, length.out = n_points)
    for (k in seq_len(nrow(ed2))) {
        x0 <- ed2$x0[k]
        y0 <- ed2$y0[k]
        x2 <- ed2$x2[k]
        y2 <- ed2$y2[k]
        mx <- (x0 + x2)/2
        my <- (y0 + y2)/2
        vx <- x2 - x0
        vy <- y2 - y0
        vlen <- sqrt(vx^2 + vy^2)
        if (!is.finite(vlen) || vlen < 1e-06) {
            x <- x0 + (x2 - x0) * t
            y <- y0 + (y2 - y0) * t
        }
        else {
            ux <- -vy/vlen
            uy <- vx/vlen
            A <- ed2$amp[k]
            x1 <- mx + ux * A
            y1 <- my + uy * A
            x <- (1 - t)^2 * x0 + 2 * (1 - t) * t * x1 + t^2 * x2
            y <- (1 - t)^2 * y0 + 2 * (1 - t) * t * y1 + t^2 * y2
        }
        res_list[[k]] <- data.frame(edge_id = paste0(ed2$from[k], "->", ed2$to[k]), from = ed2$from[k], to = ed2$to[k], 
            x = as.numeric(x), y = as.numeric(y), sim = as.numeric(ed2$sim[k]), pair_cons = as.numeric(ed2$pair_cons[k]), 
            dist_deg = as.numeric(ed2$dist_deg[k]), stringsAsFactors = FALSE, row.names = NULL)
    }
    do.call(rbind, res_list)
}
plot_edges_fierce <- function(pop_df, sim_mat, out_prefix, title, edge_mode = edge_mode, top_n = top_edge_per_pop, 
    min_sim = edge_min_sim, max_edges = max_edges, cons_col = cons_col_for_edges) {
    pop_df <- pop_df %>% dplyr::distinct(population, .keep_all = TRUE)
    world <- ggplot2::map_data("world")
    ed <- make_edge_table(pop_df, sim_mat, edge_mode = edge_mode, top_n = top_n, min_sim = min_sim, max_edges = max_edges, 
        cons_col = cons_col)
    write_result(ed, paste0(out_prefix, "_edges.tsv"), row.names = FALSE, fileEncoding = "UTF-8")
    lon_range <- range(pop_df$longitude, na.rm = TRUE)
    lat_range <- range(pop_df$latitude, na.rm = TRUE)
    xlim <- c(lon_range[1] - map_x_pad, lon_range[2] + map_x_pad)
    ylim <- c(lat_range[1] - map_y_pad, lat_range[2] + map_y_pad)
    asp <- diff(xlim)/diff(ylim)
    if (!is.finite(asp) || asp <= 0) 
        asp <- 1
    height_in <- map_base_height_in
    width_in <- height_in * asp
    show_labels <- nrow(pop_df) <= label_pop_if_nle
    if (nrow(ed) == 0) {
        message2("Edges: none (mode=%s). Plot points only: %s", edge_mode, basename(out_prefix))
        p0 <- ggplot() + geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "grey95", color = "grey80", 
            linewidth = 0.2) + geom_point(data = pop_df, aes(x = longitude, y = latitude, fill = .data[[cons_col]], 
            size = n_samples), shape = 21, color = "grey15", stroke = 0.25, alpha = 0.95) + scale_fill_rb(name = cons_col, 
            guide = guide_colorbar(order = 1)) + scale_size_continuous(range = c(1.8, 5.2), name = "N") + coord_quickmap(xlim = xlim, 
            ylim = ylim, expand = FALSE) + labs(title = title, x = "Longitude", y = "Latitude") + theme_minimal(base_size = 12) + 
            theme(plot.title = element_text(face = "bold", hjust = 0), panel.grid = element_line(color = "grey90", 
                linewidth = 0.3))
        if (show_labels) {
            p0 <- p0 + geom_text(data = pop_df, aes(x = longitude, y = latitude, label = population), size = 2.8, 
                vjust = -0.8)
        }
        save_plot(paste0(out_prefix, ".pdf"), p0, width = width_in, height = height_in, limitsize = FALSE)
        return(invisible(NULL))
    }
    path_df <- bezier_edge_points(ed, pop_df, curv_base = curv_base, curv_alpha = curv_alpha, curv_clip = curv_clip, 
        n_points = bezier_points)
    has_paths <- nrow(path_df) > 0 && all(c("x", "y", "edge_id") %in% names(path_df))
    if (!has_paths) 
        message2("Edges exist but no drawable paths; plot points only: %s", basename(out_prefix))
    if (has_paths) {
        sim_limits <- range(path_df$sim, na.rm = TRUE)
        if (length(sim_limits) != 2 || any(!is.finite(sim_limits)) || diff(sim_limits) <= 1e-06) 
            sim_limits <- c(0, 1)
    }
    else {
        sim_limits <- c(0, 1)
    }
    p <- ggplot() + geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "grey95", color = "grey80", 
        linewidth = 0.2)
    if (has_paths) {
        p <- p + geom_path(data = path_df, aes(x = x, y = y, group = edge_id, color = sim, linewidth = pair_cons), 
            alpha = 0.9, lineend = "round") + scale_color_rb(limits = sim_limits, name = "Similarity", guide = guide_colorbar(order = 1)) + 
            scale_linewidth_continuous(range = c(0.3, 3), limits = c(0, 1), name = "Pair consistency", guide = guide_legend(order = 3))
    }
    p <- p + geom_point(data = pop_df, aes(x = longitude, y = latitude, fill = .data[[cons_col]], size = n_samples), 
        shape = 21, color = "grey15", stroke = 0.25, alpha = 0.95) + scale_fill_rb(name = cons_col, guide = guide_colorbar(order = 2)) + 
        scale_size_continuous(range = c(1.8, 5.2), name = "N") + coord_quickmap(xlim = xlim, ylim = ylim, expand = FALSE) + 
        labs(title = title, x = "Longitude", y = "Latitude") + theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold", 
        hjust = 0), panel.grid = element_line(color = "grey90", linewidth = 0.3))
    if (show_labels) {
        p <- p + geom_text(data = pop_df, aes(x = longitude, y = latitude, label = population), size = 2.8, vjust = -0.8)
    }
    save_plot(paste0(out_prefix, ".pdf"), p, width = width_in, height = height_in, limitsize = FALSE)
}
