#source("forceTSP.r")
library(scales)
correlation_distance <- function(x, y = NULL, method = NULL, use = "everything")
{
  #method = match.arg(method)
  x <- t(as.matrix(x))
  if(!is.null(y))
    y <- t(as.matrix(y))
  1 - (stats::cor(x, y) + 1) / 2
}
rotateResult <- function(clust_num_EM, result, rotate = 0, back = FALSE, normalize = 0, type = "points")
{
  #pause(0.1)
  sum <- 0
  newresult <- result
  if (back)
  {
    result <- clust_num_EM + 1 - result
  }
  if (type == "points")
  {
    for (i in 1: clust_num_EM)
    {
      j <- i - rotate
      if (j <= 0)
        j <- j + clust_num_EM
      id <- which(result == j)
      newresult[id] <- sum + length(id) / 2
      sum <- sum + length(id)
    }
  }
  else
  {
    for (i in 1: clust_num_EM)
    {
      j <- i - rotate
      if (j <= 0)
        j <- j + clust_num_EM
      id <- which(result == j)
      newresult[id] <- i
    }
  }
  if (normalize == 0)
  {
    return(newresult)
  }else
  {
    if (type == "points")
    {
      return((newresult-0.5) * normalize / length(result))
    }
    else
    {
      return((newresult-0.5) * normalize / clust_num_EM)
    }
  }
}


highCorFind <- function(clust_num_EM,test_exp,stageIdxSmp,TSPFold = 2, TSPEvaluteTime = 100, TSPMethod = "arbitrary_insertion", clustMethod = "GMM", rankType = "class", debug = FALSE)
{
  stageTest <- !is.null(stageIdxSmp)
  if (stageTest)
  {
    G1IdxSmp <- stageIdxSmp[["G1"]]
    SIdxSmp <- stageIdxSmp[["S"]]
    G2MIdxSmp <- stageIdxSmp[["G2M"]]
  }
  #pause(0.1)
  #EM GMM Algorithm  Mclust()  mclust package
  if (clust_num_EM >= dim(test_exp)[1])
  {
    EM_result <- c(1: dim(test_exp)[1])
  }
  
  else if (clustMethod == "Corr2")
  {
    dist_mat <- correlation_distance(test_exp)
    # 3-dim pcoa
    pcoa_mat <- stats::cmdscale(dist_mat, k = 3) 
    scaling <- scale(pcoa_mat)#scale_uniform(pcoa_mat)#scale(pcoa_mat)
    
    fit_EM <- kmeans(scaling, clust_num_EM, iter.max = 2500, nstart = 100)
    EM_result <- fit_EM$cluster
    
  }
  
  else if (clustMethod == "Corr2")#2
  {
    dist_mat <- correlation_distance(test_exp)
    # 3-dim pcoa
    pcoa_mat <- stats::cmdscale(dist_mat, x.ret = T) 
    #scaling <- scale(pcoa_mat)#scale_uniform(pcoa_mat)#scale(pcoa_mat)
    
    fit_EM <- hclust(dist(pcoa_mat$x))
    EM_result <- cutree(fit_EM, clust_num_EM)
    #fit_EM <- kmeans(scaling, clust_num_EM, iter.max = 2500, nstart = 100)
    #EM_result <- fit_EM$cluster
    
  }
  else if (clustMethod == "Corr")
  {
    require(mclust)
    dist_mat <- correlation_distance(test_exp)
    # 20-dim pcoa
    kn <- 10#dim(test_exp)[1] - 1
    pcoa_mat <- stats::cmdscale(dist_mat, k = kn) 
    scaling <- pcoa_mat#scale_uniform(pcoa_mat)#scale(pcoa_mat)
    
    #fit_EM <- kmeans(scaling, clust_num_EM, iter.max = 2500, nstart = 100)
    #EM_result <- fit_EM$cluster
    fit_EM <- Mclust(as.matrix(scaling), G = clust_num_EM, verbose = F)
    # Avoids empty clusters from GMM
    if(length(unique(fit_EM$classification)!=clust_num_EM)){
      fit_EM <- kmeans(scaling, clust_num_EM, iter.max = 2500, nstart = 100)
      EM_result <- fit_EM$cluster
    }else{
      EM_result <- fit_EM$classification
    }
    
    
  }
  else if (clustMethod == "GMM")
  {
    require(mclust)
    fit_EM = Mclust(test_exp,G=clust_num_EM)#, modelNames = "EII")#G can be a vector, then use BIC to calculate the best components number
    EM_result <- fit_EM$classification
  }
  else if (clustMethod == "Pam")
  {
    require(cluster)
    set.seed(10)
    pamx = pam(test_exp,clust_num_EM)#G can be a vector, then use BIC to calculate the best components number
    EM_result <- pamx$clustering
    table(EM_result)
  }
  else if (clustMethod == "Kmeans")
  {
    fitKMeans <- kmeans(test_exp,centers=clust_num_EM)
    KMeansResult <- fitKMeans$cluster
    EM_result <- KMeansResult
  }
  
  #find samples in each cluster, then calculate the expression mean of cell cycle genes in each class
  obj_mean_result <- EM_result
  obj_mean_data <- test_exp
  cls_means <- c()
  for(clst_idx in 1:clust_num_EM)
  {
    clsti = which(obj_mean_result==clst_idx)
    if(length(clsti)==1)
    {
      currClsMeans <- obj_mean_data[clsti,]
    }
    else
    {
      currClsMeans <- apply(obj_mean_data[clsti,], 2, mean)
    }
    cls_means <- cbind(cls_means, currClsMeans)#mean(mean can be used in sd,sum,median and so on ),2 is to use column
  }
  
  #Get the distance matrix 
  distance_obj <- cls_means
  distance_result <- as.matrix(dist(t(distance_obj)),upper=TRUE,diag=TRUE)
  clust_num_TSP <- dim(distance_obj)[2]# <- clust_num_EM
  
  #solve TSP problem
  require(TSP)
  rownames(distance_result) <- paste("Clust",1:clust_num_TSP,sep="")
  colnames(distance_result) <- paste("Clust",1:clust_num_TSP,sep="")
  #Multi-types of TSP result to get correlation
  tourLengthLst <- c()
  highCorLst <- c()
  tourLength <- 0
  resultLst <- matrix(0, nrow = (TSPFold*clust_num_EM), ncol = length(EM_result))
  #pause(0.1)
  #####
  # assign("dm1", dist_mat, .GlobalEnv)
  # assign("mr1", EM_result, .GlobalEnv)
  # assign("rl1", resultLst, .GlobalEnv)
  # assign("dr1", distance_result, .GlobalEnv)
  # print(range(distance_result))
  # print(dim(distance_result))
  #####
  if (TSPMethod != "force")
  {
    #resultLst <- matrix(0, nrow = (TSPFold*clust_num_EM), ncol = length(EM_result))
    
    for (t in 1:(TSPFold*clust_num_EM)){
      stClsNm <- paste("Clust",(t%%clust_num_TSP)+(clust_num_TSP)*(t%%clust_num_TSP==0),sep="")
      #print(stClsNm)
      distTSP <- as.ATSP(distance_result)
      s1 <- which(labels(distTSP) == stClsNm)
      e1 <- which(labels(distTSP) == stClsNm)
      atsp <- ATSP(distance_result[-c(s1,e1),-c(s1,e1)])
      atsp <- insert_dummy(atsp, label = paste(stClsNm,stClsNm,sep="/"))
      s1_e1 <- which(labels(atsp)== paste(stClsNm,stClsNm,sep="/"))
      atsp[s1_e1, ] <- c(distance_result[-c(s1,e1), s1], 0)
      atsp[, s1_e1] <- c(distance_result[e1, -c(s1,e1)], 0)
      tour <- solve_TSP(atsp, method ="arbitrary_insertion")
      path_labels <- c(stClsNm,labels(cut_tour(tour, s1_e1)), stClsNm)
      path_idx <- match(path_labels, labels(distTSP))# Match Indices
      tourLength <- tour_length(tour)
      
      ordIndex <- head(path_idx, -1)
      ordResult <- (match(c(1:clust_num_EM), ordIndex))[EM_result]
      timeSequence <- rotateResult(clust_num_EM, ordResult, rotate = 0, type = rankType)
      tourLengthLst <- c(tourLengthLst, tourLength)
      resultLst[t, ] <- timeSequence
    }}
  
  else
  {
    #resultLst <- matrix(0, nrow = 1, ncol = length(EM_result))
    
    tmp <- forceTSP(distance_result)
    tourLength <- tmp[1]
    path_idx <- tmp[2: length(tmp)]
    path_labels <- paste("C",1:clust_num_TSP,sep="")[path_idx]
    #print(path_labels)
    #print(path_idx)
    ordIndex <- head(path_idx, -1)
    #print(ordIndex)
    ordResult <- (match(c(1:clust_num_EM), ordIndex))[EM_result]
    #print(ordResult)
    timeSequence <- rotateResult(clust_num_EM, ordResult, rotate = 0, type = rankType)
    #print(timeSequence)
    tourLengthLst <- c(tourLengthLst, tourLength)
    #print(dim(resultLst))
    resultLst[1, ] <- timeSequence
    
  }
  
  if (stageTest)
  {
    stageToNum <- rep(c(1,2,3),c(length(G1IdxSmp),length(SIdxSmp),length(G2MIdxSmp)))
    seq1 <- ordIndex
    markedTimeSequence <- stageToNum
    corLst <- c()
    
    #print(ordResult)
    for (i in 0:(length(seq1)-1))
    {
      timeSequence <- rotateResult(clust_num_EM, ordResult, rotate = i, type = rankType)
      corLst <- c(corLst,cor(markedTimeSequence,timeSequence))
      timeSequence <- rotateResult(clust_num_EM, ordResult, rotate = i, back = TRUE)
      corLst <- c(cor(markedTimeSequence,timeSequence),corLst)
    }
    
    max(corLst)
    which.max(corLst)
    if(which.max(corLst) > clust_num_EM)
    {
      i <- which.max(corLst)-clust_num_EM-1
      timeSequence <- rotateResult(clust_num_EM, ordResult, rotate = i, type = rankType)
      highCor <- cor(markedTimeSequence,timeSequence)
      timeSequence <- rotateResult(clust_num_EM, ordResult, rotate = i, type = "class")
    }else
    {
      i <- (clust_num_EM - which.max(corLst))
      timeSequence <- rotateResult(clust_num_EM, ordResult, rotate = i, back = TRUE, type = rankType)
      highCor <- cor(markedTimeSequence,timeSequence)
      timeSequence <- rotateResult(clust_num_EM, ordResult, rotate = i, back = TRUE, type = "class")
    }
    highCorLst <- c(highCorLst,highCor)
  }
  #tourLengthLst <- c(tourLengthLst, tourLength)
  #resultLst[t, ] <- timeSequence
  if (TSPMethod == "force" & debug)
  {
    message("force")
    break
  }
  
  
  if (debug)
  {
    message(Sys.time())
  }
  minLengthID <- which(tourLengthLst == min(tourLengthLst))[1]
  if (stageTest)
  {
    bestTour <- list(highCor = highCorLst[minLengthID], result = resultLst[minLengthID, ])
  }
  else
  {
    #print(minLengthID)
    bestTour <- list(highCor = NULL, result = resultLst[minLengthID, ])
  }
  
  return(bestTour)
}
