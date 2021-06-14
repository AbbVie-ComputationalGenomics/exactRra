#' Functions for computing exact RRA P values


#' Compute a matrix of per-rank RRA statistics by approximating relative ranks
#' as uniform distribution variables.
#'
#' Compute a lookup table F of shape (M, N), where M is the total number of
#' ranks to be drawn from ranks {1, 2, ..., N} without replacement and F[m, n]
#' indicates probability that out of a set of M randomly selected ranks from, at
#' least m of the selected ranks have a rank of n or lower. This probability is
#' computed by approximating the relative ranks n / N using independent uniform
#' distributions.
#'
#' In other words,
#' F[m, n] = pbeta(n / N, 1, M + 1 - m), where
#' 1 <= m <= M and 1 <= m <= N.
#'
#' @param size Number of ranks to be drawn, i.e. number of rows in the output
#'     matrix.
#' @param max_rank Total number of available ranks to be drawn from, i.e. the
#'     number of columns in the output matrix.
#' @param log_p Whether to return log-probabilities.
#'
#' @return A M * N probability lookup matrix.
#'
#' @export
compute_pbeta_table = function(size, max_rank, log_p = F) {
    pbeta_mat = t(sapply(
        seq(size),
        function(i) {
            pbeta(seq(max_rank) / max_rank, i, size + 1 - i,
                  log.p = log_p)
        }
    ))
    pbeta_mat
}


#' Compute a matrix of exact per-rank RRA statistics using combinatorics.
#'
#' Helper function for computing a table F of shape (M, N), where F[m, n] is the
#' number of ways in which M ranks can be chosen from {1, 2, ..., N} such that
#' at least m of them have a rank equal or lower to n.
compute_combinations_table = function(size, max_rank) {
    # combos_mat[m, n] corresponds to the number of ways in which M ranks can be
    # chosen such that exactly m ranks are less than or equal to n and exactly
    # M - m ranks are greater than n.
    combos_mat = t(sapply(
        seq(size),
        function(m) {
            choose(seq(max_rank), m) * choose(max_rank - seq(max_rank),
                                              size - m)
        }
    ))

    # cum_combos_mat[m, n] corresponds to the number of ways in which M ranks
    # can be chosen such that at least m out of M ranks are less than or equal
    # to n. This is obtained by cum_combos_mat[m, n] = sum(combos_mat[m:M, n]).
    #
    # Steps:
    # 1. Reverse the rows in combos_mat.
    # 2. Calculate cumulative sums along each column.
    # 3. Reverse the rows again, so that the cumulative sums are from the end
    #    to the beginning of each column.
    cum_combos_mat = apply(
        combos_mat[rev(seq(size)), ],  # Reverse the rows.
        2,
        cumsum
    )[rev(seq(size)), ]  # Re-reverse the rows.

    if (any(cum_combos_mat == Inf)) {
        msg = paste0("Got infinite values in the resulting matrix. Please use ",
                     "compute_pbeta_table() instead.")
        stop(msg)
    }

    cum_combos_mat
}


#' Helper functions for generating a memoized order statistic table of shape
#' (size, max_rank).
#'
#' .stat_table_memoized() is a memoised version of .stat_table().
#'
#' @export
.stat_table = function(size, max_rank, stat = c("beta", "comb")) {
    stat = match.arg(stat)
    if (stat == "beta") {
        out_mat = compute_pbeta_table(size, max_rank)
    } else {
        out_mat = compute_combinations_table(size, max_rank)
    }
    out_mat
}
#' @describeIn dot-stat_table
#'
#' @export
.stat_table_memoized = memoise::memoise(.stat_table)


#' Helper function for computing RRA P value from a score matrix.
#'
#' Computes the exact probability of randomly choosing M of ordered ranks
#' (r_1, ..., r_M) that satisfy the following conditions. score_mat is a matrix
#' of shape (M, N).
#' 1. r_1 <= r_2 <= ... <= r_M.
#' 2. min(score_mat[m, r_m]) > cutoff for all m in {1, ..., M}.
#'
#' @param score_mat A M * N score matrix.
#' @param cutoffs Significance cutoffs to be used for the ranks.
#' @return P values corresponding to each provided cutoff.
#'
#' @md
.exact_rra_pval_from_score_mat = function(score_mat, cutoffs) {
    score_mat = as.matrix(score_mat)
    M = nrow(score_mat)
    N = ncol(score_mat)

    # max_j[m] is the value such that
    #     all(score_mat[m, 1:(max_j[m])] <= cutoff) == TRUE, and
    #     all(score_mat[m, -(1:(max_j[m]))] > cutoff) == TRUE.
    pvals = sapply(
        cutoffs,
        function(cutoff) {
            max_j = exactRra::matrix_binary_search(score_mat, cutoff)
            less_signif_log_combos = exactRra:::trunc_it_log_cumsum_cpp(
                max_j, N)
            total_log_combos = lchoose(N, M)
            1 - exp(less_signif_log_combos - total_log_combos)
        }
    )

    pvals
}


#' Compute the exact RRA P value for a rank vector (r_1, ..., r_m).
#'
#' The robust rank aggregation statistic (RRA) is defined as min(stat(i, r_i)).
#' This function computes the exact P value for an observed ranks vector r =
#' (r_1, ..., r_m) under the null hypothesis that the ranks are drawn uniformly
#' from ranks 1 ... n (without replacement).
#'
#' @param ranks A vector of potentially unordered ranks. NAs are not allowed.
#' @param max_rank The maximum value a rank can have, i.e. the total number of
#'     ranked values in the experiments.
#' @param stat_method If "beta", then the test statistic is obtained by via
#'     k'th order statistics by approximating normalised ranks as the uniform
#'     distribution, which follows beta distribution. If "comb", then the test
#'     statistic is obtained exactly through combinatorics.
#' @param memoize Whether to memoize the computed probability lookup tables.
#'
#' @return A vector of the test statistic, the P value and the index in ranks
#' that produced the final test statistic.
#'
#' @export
exact_rra_pval = function(ranks, max_rank, stat_method = c("beta", "comb"),
                          memoize = T) {
    # Sort each set of ranks (rows in ranks) and assign NAs last in each vector.
    if (any(is.na(ranks))) {
        stop("NAs are not allowed in argument <ranks>.")
    }
    ranks = sort(ranks)
    n_ranks = length(ranks)

    if (memoize) {
        stat_table = .stat_table_memoized(n_ranks, max_rank, stat_method)
    } else {
        stat_table = .stat_table(n_ranks, max_rank, stat_method)
    }

    # Calculate the overall RRA statistic for the current set of ranks.
    rho = stat_table[cbind(1:n_ranks, ranks)]
    rho_argmin = which.min(rho)
    rho_min = rho[rho_argmin]

    pval = .exact_rra_pval_from_score_mat(stat_table, rho_min)
    c(
        rho = rho_min,
        pval = pval,
        rho_argmin = rho_argmin
    )
}
