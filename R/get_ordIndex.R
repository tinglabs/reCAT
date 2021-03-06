#' Perform reCAT
#'
#' This function performs reCAT on your data, generating the order index of each cell as output.
#'@param test_exp The gene expression input obtained from get_test_exp
#'@param threadnum The number of threads to use
#'@param step_size The interval of cycles to merge in the algorithm. Defaults to 2
#'@param base_cycle_range Sets the cluster number used to formulate the base cycle, can be set as c(6:9) or c(7:10). Defaults to c(6:9)
#'@param max_loop This sets the maximum number for not applying step_size in reCAT, i.e. if N = 100, step_size = 5, and max_loop = 50. Then reCAT merges consecutive cycles from 10 to 50, then merges cycles by intervals of 5 from 55 to 100.
#'@param clust_method There are two clustering approaches: GMM ("GMM") and correlation-based ("Corr"). This defaults to NULL, where reCAT uses "GMM" for datasets with less than 300 cells and "Corr" for > 300. We recommend "GMM" for better accuracy but only for small datasets (as it fails to cluster on larger datasets), and "Corr" for efficiency and larger datasets.
#'@export
#'@examples
#'get_ordIndex(data, threadnum = 2, step_size = 2, base_cycle_range = c(6:9))


# get the order using reCAT
#
# you can choose thread you want to run the order (threadnum)


#source("bestEnsembleComplexTSP.r")

get_ordIndex <- function(test_exp, threadnum, step_size = 2, base_cycle_range = c(6:9), max_loop = NULL, clust_method = NULL)
                        
{
  if(!is.null(clust_method)){
    clust_method=clust_method
  }else{
    if(nrow(test_exp) < 300){
      clust_method = "GMM"
    }else{
      clust_method = "Corr"
    }
  }
  result <- bestEnsembleComplexTSP(test_exp = test_exp, TSPFold = 2, beginNum = 10, endNum = dim(test_exp)[1], threads = threadnum, 
                                   base_cycle_range = base_cycle_range, step_size = step_size, max_loop = max_loop, clustMethod = clust_method, max_num = 300)
  ord <- result$ensembleResultLst[dim(result$ensembleResultLst)[1], ]
  #print(ord)
  ordIndex <- order(ord)
  return(ordIndex)
}
