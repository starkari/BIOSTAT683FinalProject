library(foreach)
library(doParallel)
B = 1000
n <- nrow(ObsData)
registerDoParallel(cores = detectCores())
estimates <-
  foreach(b = 1:B, .combine = rbind) %dopar% {
    # sample the indices 1 to n with replacement
    bootIndices<- sample(1:n, replace=T)
    bootData <- ObsData[bootIndices,]
    # calling the above function
    output <- run.tmle(ObsData=bootData, SL.library=SL.library)$estimates
    return(output)
  }

colnames(estimates)<-c("SimpSubs", "IPTW", "TMLE")
