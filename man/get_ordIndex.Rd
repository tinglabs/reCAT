% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_ordIndex.R
\name{get_ordIndex}
\alias{get_ordIndex}
\title{Perform reCAT}
\usage{
get_ordIndex(
  test_exp,
  threadnum,
  step_size = 2,
  base_cycle_range = c(6:9),
  max_loop = NULL,
  clust_method = NULL
)
}
\arguments{
\item{test_exp}{The gene expression input obtained from get_test_exp}

\item{threadnum}{The number of threads to use}

\item{step_size}{The interval of cycles to merge in the algorithm. Defaults to 2}

\item{base_cycle_range}{Sets the cluster number used to formulate the base cycle, can be set as c(6:9) or c(7:10). Defaults to c(6:9)}

\item{max_loop}{This sets the maximum number for not applying step_size in reCAT, i.e. if N = 100, step_size = 5, and max_loop = 50. Then reCAT merges consecutive cycles from 10 to 50, then merges cycles by intervals of 5 from 55 to 100.}

\item{clust_method}{There are two clustering approaches: GMM ("GMM") and correlation-based ("Corr"). This defaults to NULL, where reCAT uses "GMM" for datasets with less than 300 cells and "Corr" for > 300. We recommend "GMM" for better accuracy but only for small datasets (as it fails to cluster on larger datasets), and "Corr" for efficiency and larger datasets.}
}
\description{
This function performs reCAT on your data, generating the order index of each cell as output.
}
\examples{
get_ordIndex(data, threadnum = 2, step_size = 2, base_cycle_range = c(6:9))
}
