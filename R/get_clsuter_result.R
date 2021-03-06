#source("get_score.R")

#' Obtain the clustering result
#'
#' This function obtains the clustering result from reCAT, bayes and mean scores, and the time series ordering.
#'@param test_exp The pre-processed input data
#'@param ensembleResultLst This can be obtained from the automatically saved .RData file in your working directory.
#'@param resultLst This can be obtained from the automatically saved .RData file in your working directory.
#'@param cls_num The desired number of clusters to view (value > 10)
#'@export
#'@examples
#'test_exp <- get_test_exp(data)
#'data_ordIndex <- get_ordIndex(test_exp, 2)
#'# Load the following .RData file from your local working directory, it is automatically generated when you perform the function 'get_ordIndex'
#'load("../data/bestEnsembleComplexTSP 10 - N  .RData")
#'cls_result = get_cluster_result(test_exp = test_exp, ensembleResultLst = ensembleResultLst,
#'                                resultLst = resultLst, cls_num = 20)

get_cluster_result <- function(test_exp, ensembleResultLst, resultLst, cls_num)
{
	if (cls_num <= 10)
	{
		print("the cluster number is too small!")
		return(0)
	}
	else
	{
		rownum = cls_num-10+1

		ordIndex2 = ensembleResultLst[rownum, ]
		EM_result = resultLst[rownum, ]
		cls_mean = c()
		for(i in c(1:cls_num))
		{
		  id = which(EM_result == i)
		  tmp = sum(complex(argument = ordIndex2[id] * 2 * pi - pi))
		  cls_mean = cbind(cls_mean, tmp)
		}
		meanResults <- Arg(cls_mean) / 2 / pi
		ordIndex3 = order(meanResults)

		cls_mean = c()
		for(i in c(1:cls_num))
		{
		  id = which(EM_result == i)

		  if (length(id) == 1)
		  {
		    tmp = test_exp[id, ]
		  }
		  else
		  {
		    tmp = apply(test_exp[id, ], 2, mean)
		  }
		  cls_mean = cbind(cls_mean, tmp)
		}

		score_r = get_score(cls_mean)

		return(list(cls = EM_result, bayes_score = score_r$bayes_score, mean_score = score_r$mean_score, ordIndex = ordIndex3))
	}
}
