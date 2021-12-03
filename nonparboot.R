library(foreach)
library(doParallel)
library(tidyverse)
library(SuperLearner)
library(ltmle)
library(dplyr)

data <-read.csv("fertility_sperm.csv")
colnames(data)

data <- data %>% 
  mutate(W11 = age,
         W12 = sitting_hour,
         W13 = trauma,
         W14 = case_when(alcohol<0.8 ~ '2',
                         alcohol==0.8 ~ '1',
                         alcohol==1 ~ '0'),
         W2 = case_when(season==-1 ~ 'Fall/Winter',
                        season==-0.33 ~'Spring/Summer',
                        season==0.33 ~ 'Spring/Summer',
                        season==1 ~ 'Fall/Winter'),
         A = case_when(smoking==-1 ~ 0,
                       smoking==0 ~ 1,
                       smoking==1 ~ 1),
         Y = ifelse(diagnosis=='N',1,0))

ObsData <- data %>% dplyr::select(W11, W12, W13, W14, W2, A, Y)

set.seed(123)


SL.library<- c('SL.glm', 'SL.glm.interaction', "SL.step",
               "SL.randomForest","SL.step.forward","SL.stepAIC","SL.mean")


run.tmle <- function(ObsData, SL.library){
  
  #------------------------------------------
  # Simple substitution estimator 
  #------------------------------------------
  
  # dataframe X with baseline covariates and exposure
  X <- subset(ObsData, select=c(A, W11, W12, W13, W14,W2))
  
  # set the exposure=1 in X1 and the exposure=0 in X0
  X1 <- X0 <- X
  X1$A <- 1 
  X0$A <- 0
  
  # Estimate E_0(Y|A,W) with Super Learner
  SL.outcome <- SuperLearner(Y=ObsData$Y, X=X, SL.library=SL.library,
                             family="binomial", cvControl=list(V=10,stratifyCV=TRUE))
  
  # get the expected outcome, given the observed exposure and covariates
  expY.givenAW <- predict(SL.outcome, newdata=ObsData)$pred
  # expected outcome, given A=1 and covariates 
  expY.given1W <- predict(SL.outcome, newdata=X1)$pred
  # expected outcome, given A=0 and covariates
  expY.given0W <- predict(SL.outcome, newdata=X0)$pred
  
  # simple substitution estimator would be 
  PsiHat.SS <- mean(expY.given1W - expY.given0W)
  
  #------------------------------------------
  # Inverse probability of txt weighting
  #------------------------------------------
  
  #  Super Learner for the exposure mechanism  P_0(A=1|W)
  SL.exposure <- SuperLearner(Y=ObsData$A, 
                              X=subset(ObsData, select= -c(A,Y,W2)),
                              SL.library=SL.library, family="binomial",
                              cvControl=list(V=10))
  
  # generate the predicted prob of being exposed, given baseline cov
  probA1.givenW <- SL.exposure$SL.predict
  # generate the predicted prob of not being exposed, given baseline cov
  probA0.givenW <- 1- probA1.givenW
  
  # clever covariate
  H.AW <- as.numeric(ObsData$A==1)/probA1.givenW - as.numeric(ObsData$A==0)/probA0.givenW
  
  # also want to evaluate the clever covariate at A=1 and A=0 for all participants
  H.1W <- 1/probA1.givenW
  H.0W <- -1/probA0.givenW
  
  # IPTW estimate
  PsiHat.IPTW <- mean(H.AW*ObsData$Y)
  
  #------------------------------------------
  # Targeting & TMLE
  #------------------------------------------
  
  # Update the initial estimator of E_0(Y|A,W)
  # run logistic regression of Y on H.AW using the logit of the esimates as offset
  logitUpdate<- glm( ObsData$Y ~ -1 +offset(qlogis(expY.givenAW)) + 
                       H.AW, family='binomial')
  epsilon <- logitUpdate$coef
  
  # obtain the targeted estimates
  expY.givenAW.star<- plogis( qlogis(expY.givenAW)+ epsilon*H.AW )  
  expY.given1W.star<- plogis( qlogis(expY.given1W)+ epsilon*H.1W )	
  expY.given0W.star<- plogis( qlogis(expY.given0W)+ epsilon*H.0W )
  
  # TMLE point estimate
  PsiHat.TMLE<- mean(expY.given1W.star - expY.given0W.star)
  
  #------------------------------------------
  # Return a list withthe point estimates, targeted estimates of E_0(Y|A,W), 
  # and the vector of clever covariates
  #------------------------------------------
  
  estimates <- data.frame(cbind(PsiHat.SS=PsiHat.SS, PsiHat.IPTW, PsiHat.TMLE))
  predictions <- data.frame(cbind(expY.givenAW.star, expY.given1W.star, expY.given0W.star))
  colnames(predictions) <- c('givenAW', 'given1W', 'given0W')
  list(estimates=estimates, predictions=predictions, H.AW=H.AW)
}



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

save(estimates, file='boot_par.Rdata')
