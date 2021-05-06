################################################################################
##100 times permutaions
################################################################################

rfPermutation <- function(data.cv,cate.cv,fold,return.list=T){
  message('100 runs of CV.')
  error.cv <- vector(mode = 'numeric')
  pb <- txtProgressBar(style=3)
  result=data.frame(row.names = rownames(data.cv))
  for (run in 1:100) {
    setTxtProgressBar(pb, run/100)
    sig <- -1
    n <- floor(nrow(data.cv)/fold)
    sam <- sample(1:nrow(data.cv),n,replace = F,)
    data.train <- data.cv[-sam,]
    data.test <- data.cv[sam,]
    cate.train <- cate.cv[-sam]
    cate.test <- cate.cv[sam]
    min.cluster.num <- min(table(cate.train))
    num.of.cluster <- length(unique(cate.train))
    if (min.cluster.num<5|num.of.cluster!=length(unique(cate.cv))) {
      sig <- 50
      while (sig>0) {
        sam <- sample(1:nrow(data.cv),n,replace = F)
        data.train <- data.cv[-sam,]
        data.test <- data.cv[sam,]
        cate.train <- cate.cv[-sam]
        cate.test <- cate.cv[sam]
        sig <- sig-1
        min.cluster.num <- min(table(cate.train))
        if (min.cluster.num>=5) {
          sig <- -1
        }
      }
    }
    if (sig==0) {
      message('Data set too small!')
      message('Rfcv stops at:')
      message(run)
      break
    }
    if (sig==-1) {
      model.cv <- randomForest(x = data.train,
                               y = cate.train,
                               proximity = F,
                               importance = T,
      )
      pred <- predict(object = model.cv,newdata = data.cv)
      result[,run] <- pred
      error.cv <- c(error.cv,mean(pred!=cate.cv))
      
    }
  }
  close(pb)
  if (return.list) {
    return(list(data=data.cv,cate=cate.cv,error=error.cv,result=result))
  }else{
    return(result)
  }
  
}
