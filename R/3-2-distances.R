TVD <- function(p, q) {
    0.5 * sum(abs(p - q))
}
JSD_pair <- function(p, q, eps = jsd_eps) {
    p <- p + eps
    p <- p/sum(p)
    q <- q + eps
    q <- q/sum(q)
    m <- 0.5 * (p + q)
    KL <- function(a, b) sum(a * log(a/b))
    sqrt(0.5 * KL(p, m) + 0.5 * KL(q, m))
}
pairwise_matrix <- function(mat, FUN) {
    n <- nrow(mat)
    D <- matrix(0, n, n)
    for (i in seq_len(n)) {
        D[i, i] <- 0
        if (i < n) {
            for (j in (i + 1):n) {
                v <- FUN(mat[i, ], mat[j, ])
                D[i, j] <- v
                D[j, i] <- v
            }
        }
    }
    rownames(D) <- rownames(mat)
    colnames(D) <- rownames(mat)
    D
}
.nei_I <- function(p, q) {
    num <- sum(p * q)
    den <- sqrt(sum(p^2) * sum(q^2))
    if (den == 0) 
        return(NA_real_)
    num/den
}
.nei_D <- function(p, q) {
    I <- .nei_I(p, q)
    if (is.na(I) || I <= 0) 
        return(NA_real_)
    -log(I)
}
.DA_pair <- function(p, q) {
    1 - sum(sqrt(p * q))
}
pairwise_neiD <- function(mat) {
    pairwise_matrix(mat, .nei_D)
}
pairwise_DA <- function(mat) {
    pairwise_matrix(mat, .DA_pair)
}
compute_pairwise_Rst <- function(ystr_df, pop_col = "population") {
    stopifnot(pop_col %in% names(ystr_df))
    loci <- setdiff(colnames(ystr_df), c("sample_id", pop_col))
    loci <- loci[sapply(ystr_df[loci], function(x) is.numeric(x) || is.integer(x))]
    pops <- unique(ystr_df[[pop_col]])
    P <- length(pops)
    R <- matrix(NA_real_, P, P, dimnames = list(pops, pops))
    rst_two <- function(xA, xB) {
        xA <- xA[is.finite(xA)]
        xB <- xB[is.finite(xB)]
        if (length(xA) < 2 || length(xB) < 2) 
            return(NA_real_)
        nA <- length(xA)
        nB <- length(xB)
        nT <- nA + nB
        mA <- mean(xA)
        mB <- mean(xB)
        mT <- (nA * mA + nB * mB)/nT
        SSB <- nA * (mA - mT)^2 + nB * (mB - mT)^2
        SSW <- sum((xA - mA)^2) + sum((xB - mB)^2)
        if ((SSB + SSW) <= 0) 
            return(NA_real_)
        SSB/(SSB + SSW)
    }
    for (i in seq_len(P)) {
        for (j in seq_len(P)) {
            if (j < i) {
                R[i, j] <- R[j, i]
                next
            }
            if (i == j) {
                R[i, j] <- 0
                next
            }
            rloc <- c()
            for (L in loci) {
                xi <- ystr_df[ystr_df[[pop_col]] == pops[i], L][[1]]
                xj <- ystr_df[ystr_df[[pop_col]] == pops[j], L][[1]]
                rloc <- c(rloc, rst_two(xi, xj))
            }
            if (length(rloc)) 
                R[i, j] <- mean(rloc, na.rm = TRUE)
            if (is.nan(R[i, j])) 
                R[i, j] <- NA_real_
        }
    }
    R
}
sanitize_dist_matrix <- function(D, name = "D", eps = 1e-06) {
    if (is.null(D)) 
        stop("distance matrix is NULL: ", name)
    D <- as.matrix(D)
    if (!isSymmetric(D, tol = 1e-08)) 
        D <- (D + t(D))/2
    if (any(!is.finite(D))) {
        finite_vals <- D[is.finite(D)]
        repl <- if (length(finite_vals)) 
            max(finite_vals) + eps
        else 1
        D[!is.finite(D)] <- repl
    }
    diag(D) <- 0
    up <- D[upper.tri(D)]
    if (length(up) > 1) {
        if (all(abs(up - up[1]) < 1e-12)) {
            n <- nrow(D)
            if (n >= 2) {
                for (i in seq_len(n)) for (j in seq_len(n)) if (i != j) 
                  D[i, j] <- abs(i - j) * eps
                diag(D) <- 0
            }
        }
    }
    return(D)
}
