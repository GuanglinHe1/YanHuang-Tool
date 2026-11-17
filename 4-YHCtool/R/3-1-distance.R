source("R/1-2-runtime.R")
for (HH in HH_range) {
    f <- file.path(frequency, sprintf("freq_matrix_HH%02d.tsv", HH))
    if (!file.exists(f)) 
        next
    tab <- read_table(f)
    mat <- as.matrix(tab[, -1, drop = FALSE])
    rownames(mat) <- tab$population
    matrices <- list(TVD = pairwise_matrix(mat, TVD), JSDist = pairwise_matrix(mat, JSD_pair), NeiD = pairwise_neiD(mat), 
        DA = pairwise_DA(mat), BrayCurtis = as.matrix(vegan::vegdist(mat, method = "bray")), Euclidean = as.matrix(dist(mat)), 
        Manhattan = as.matrix(dist(mat, method = "manhattan")))
    for (metric in names(matrices)) write_result(matrices[[metric]], file.path(output, sprintf("pairwise_%s_HH%02d.tsv", 
        metric, HH)))
    if (nzchar(metadata)) {
        meta <- read_table(metadata)
        if (anyDuplicated(meta$population)) 
            stop("Duplicate metadata populations")
        meta <- meta[match(rownames(mat), meta$population), , drop = FALSE]
        if (anyNA(meta$population)) 
            stop("Missing metadata populations")
        D <- matrices$TVD
        if (all(c("latitude", "longitude") %in% names(meta))) {
            if (anyNA(meta[c("latitude", "longitude")])) 
                stop("Invalid coordinates")
            G <- as.matrix(dist(meta[c("latitude", "longitude")]))
            test <- vegan::mantel(as.dist(D), as.dist(G), permutations = permutations)
            write_result(data.frame(statistic = unname(test$statistic), p = test$signif), file.path(output, sprintf("Mantel_HH%02d.tsv", 
                HH)))
        }
        if ("region" %in% names(meta)) {
            test <- vegan::adonis2(as.dist(D) ~ region, data = meta, permutations = permutations)
            write_result(data.frame(term = rownames(test), as.data.frame(test)), file.path(output, sprintf("PERMANOVA_HH%02d.tsv", 
                HH)))
        }
    }
}
if (nzchar(ystr)) {
    z <- read_table(ystr)
    if (!"population" %in% names(z)) 
        stop("Missing Y-STR population column")
    if (!any(vapply(z[setdiff(names(z), c("sample_id", "population"))], is.numeric, logical(1)))) 
        stop("No numeric Y-STR loci")
    write_result(compute_pairwise_Rst(z), file.path(output, "pairwise_Rst.tsv"))
} else message("Y-STR input omitted; skipping Rst")
