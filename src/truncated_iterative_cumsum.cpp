
#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;


//' Helper function for computing a truncated iterative cumulative sum.
//'
//' Let max_j = J be natural numbers 0 <= J[0] <= ... <= J[m] <= n.
//' The truncated iterative cumulative sum is defined as
//' sum_{i_{n-1} = J_{n-1}}^{n-1} sum_{i_{n-2} = J_{n-2}}^{J_{n-1}-1} ... sum_{i_0 = J_0}^{J_1-1} 1.
//'
//' @param max_j Natural numbers 0 <= max_j[1] <= ... <= max_j[m] <= n.
//' @param n The number of values to be cumulatively summed. See description.
//'
//' @return The calculated sum or log-sum.
//'
// [[Rcpp::export]]
int trunc_it_cumsum_cpp(NumericVector max_j, int n) {
    // Special case.
    if (any(max_j >= n).is_true()) {
        return 0;
    }

    // Create the cumulative sums matrix. Pad with a constant value to make
    // the matrix indexing 1-based.
    NumericMatrix cumsums(max_j.length() + 1, n + 1);
    std::fill(cumsums.begin(), cumsums.end(), 0.0);
    cumsums(0, _) = NumericVector(n + 1, 1.0);

    // As mentioned above, i and j are now 1-based - same as max_j.
    for (int i = 1; i <= max_j.length(); i++) {
        int start = (max_j[i - 1] < 0 ? 1 : max_j[i - 1] + 1);
        for (int j = start; j <= n; j++) {
            if (j == start) {
                cumsums(i, j) = cumsums(i - 1, j - 1);
            } else {
                cumsums(i, j) = cumsums(i - 1, j - 1) + cumsums(i, j - 1);
            }
        }
    }

    return(cumsums(max_j.length(), n));
}


//' @describeIn trunc_it_cumsum_cpp Helper function for
//' computing a truncated iterative cumulative log-sum.
//'
// [[Rcpp::export]]
double trunc_it_log_cumsum_cpp(NumericVector max_j, int n) {
    // Special case.
    if (any(max_j >= n).is_true()) {
        return(R_NegInf);
    }

    // Create the cumulative sums matrix. Pad with a constant value to make
    // the matrix indexing 1-based.
    NumericMatrix cumsums(max_j.length() + 1, n + 1);
    std::fill(cumsums.begin(), cumsums.end(), 0.0);
    cumsums(0, _) = NumericVector(n + 1, 1.0);

    // Store the order of magnitude here to keep the significant digits close to
    // zero.
    double cur_order_of_magnitude = 0.0;
    double total_order_of_magnitude = 0.0;

    // As mentioned above, i and j are now 1-based - same as max_j.
    for (int i = 1; i <= max_j.length(); i++) {
        int start = (max_j[i - 1] < 0 ? 1 : max_j[i - 1] + 1);
        for (int j = start; j <= n; j++) {
            if (j == start) {
                cumsums(i, j) = (cumsums(i - 1, j - 1)
                                    / exp(cur_order_of_magnitude));
            } else {
                cumsums(i, j) = cumsums(i - 1, j - 1)
                                    / exp(cur_order_of_magnitude)
                                + cumsums(i, j - 1);
            }
        }

        // Keep the order of magnitude of the sum close to 0 and update
        // cur_order_of_magnitude to be used in the next row.
        total_order_of_magnitude += cur_order_of_magnitude;
        cur_order_of_magnitude = floor(log(cumsums(i, n)));
    }

    return(log(cumsums(max_j.length(), n)) + total_order_of_magnitude);
}
