#source("get_hmm.R")
library(doParallel)

#' Determining the starting cell
#'
#' This function is used to find the starting cell of the cycle
#'@param bayes_score Obtained from get_score function
#'@param mean_score Obtained from get_score function
#'@param rdata This sets the region of each cycle stage manually by observing the scores. Usage: rdata = t(data.frame(c(2,5), c(10,15), c(20,25))) corresponding to G1: 2-5, S: 10-15 and G2/M: 20-25. If NULL, it will prompt you to enter the region values such as '2 5', please enter twice to confirm each line of input.
#'@param ordIndex Obtained from get_ordIndex function
#'@param cls_num The number of cell cycle stages that you want to segment, ie 3 - G1, S, G2/M.
#'@param nthreads The number of threads. Defaults to 3
#'@export
#'@examples
#'test_exp <- get_test_exp(data)
#'data_ordIndex <- get_ordIndex(test_exp, 2)
#'score_results <- get_score(t(test_exp))
#'rdata = t(data.frame(c(2,5), c(10,15), c(20,25)))
#'data_start <- get_start(score_results$bayes_score, score_results$mean_score, data_ordIndex, cls_num = 3, rdata,  nthread = 3)

get_start <- function(bayes_score, mean_score, ordIndex, cls_num, rdata = NULL, nthread = 3)
{
	cl <- makeCluster(nthread)
	registerDoParallel(cl)
	le = length(ordIndex)
	start = 0
	p = -Inf
	for(i in 1:le)
	{
		myord = c (i:1, le:(i+1))
		if (i == le)
			myord = c(le:1)

	    re = get_hmm_order(bayes_score, mean_score, ordIndex, cls_num, myord, rdata)
		#if (max(re$p) > p)
	  if (max(re) > p)
		{
			start = i
			#p = max(re$p)
			p = max(re)
		}
	}

	#p = apply(log_lk, 2, max)
	#start = which(p == max(p))

	return(start)
}
