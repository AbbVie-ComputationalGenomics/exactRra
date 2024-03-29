% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{trunc_it_cumsum_cpp}
\alias{trunc_it_cumsum_cpp}
\alias{trunc_it_log_cumsum_cpp}
\title{Helper function for computing a truncated iterative cumulative sum.}
\usage{
trunc_it_cumsum_cpp(max_j, n)

trunc_it_log_cumsum_cpp(max_j, n)
}
\arguments{
\item{max_j}{Natural numbers 0 <= max_j[1] <= ... <= max_j[m] <= n.}

\item{n}{The number of values to be cumulatively summed. See description.}
}
\value{
The calculated sum or log-sum.
}
\description{
Let max_j = J be natural numbers 0 <= J[0] <= ... <= J[m] <= n.
The truncated iterative cumulative sum is defined as
sum_{i_{n-1} = J_{n-1}}^{n-1} sum_{i_{n-2} = J_{n-2}}^{J_{n-1}-1} ... sum_{i_0 = J_0}^{J_1-1} 1.
}
\section{Functions}{
\itemize{
\item \code{trunc_it_log_cumsum_cpp}: Helper function for
computing a truncated iterative cumulative log-sum.
}}

