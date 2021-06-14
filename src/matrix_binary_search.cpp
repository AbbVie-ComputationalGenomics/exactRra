
#include <Rcpp.h>
using namespace Rcpp;

//' Binary search for each row in two-way monotonous matrix score_mat
//'
//' Matrix `score_mat` is a M * N matrix, which is monotonous in two ways in the
//' sense that score_mat[m, n] < score_mat[m, n + 1] and score_mat[m, n] >
//' score_mat[m + 1, n] for each 1 <= m <= M and 1 <= n <= N.
//'
//' @param score_mat A M * N matrix of values to be searched.
//' @param cutoff The cutoff value.
//' @return A vector J, where value J[m] is such that score_mat[m, n] <= cutoff
//' for n <= J[m] and score_mat[m, n] > cutoff for n > J[m]. Note that values in
//' J are 1-based.
//'
//' @md
//' @export
// [[Rcpp::export]]
NumericVector matrix_binary_search(NumericMatrix score_mat, double cutoff) {
    int M = score_mat.nrow();
    int N = score_mat.ncol();
    NumericVector J(M);
    int a = 0;
    int b = N - 1;
    int mid;

    // Index a is reused in each iteration m thanks thanks to the monotonicity
    // of score_mat along each column.
    for (int m = 0; m < M; m++) {
        // Edge cases.
        if (score_mat(m, 0) > cutoff) {
            J[m] = 0;  // Values in J are 1-based.
            continue;
        }
        if (score_mat(m, N - 1) <= cutoff) {
            J[m] = N;  // Values in J are 1-based.
            continue;
        }

        while (b > a + 1) {
            mid = (a + b) / 2;
            if (score_mat(m, mid) > cutoff) {
                b = mid;
            } else if (score_mat(m, mid) < cutoff) {
                a = mid;
            } else {
                // score_mat(m, mid) == cutoff
                break;
            }
        }

        // Now either score_mat(m, mid) == cutoff or a >= b - 1. Move a forward
        // to find the supremum a that score_mat(m, a) <= cutoff.
        while (a + 1 < N && score_mat(m, a + 1) <= cutoff) {
            a++;
        }

        J[m] = a + 1;  // Values in J are 1-based.
    }

    return J;
}
