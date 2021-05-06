################################################################################
##Plain cross validation
################################################################################

myrfcv <- function(data.cv,cate.cv,fold,return.list=T,e=F){
  message('100 runs of CV.')
  error.cv <- vector(mode = 'numeric')
  pb <- txtProgressBar(style=3)
  for (run in 1:100) {
    setTxtProgressBar(pb, run/100)
    sig <- -1
    n <- floor(nrow(data.cv)/fold)
    sam <- sample(1:nrow(data.cv),n,replace = F)
    data.train <- data.cv[-sam,]
    data.test <- data.cv[sam,]
    cate.train <- cate.cv[-sam]
    cate.test <- cate.cv[sam]
    min.cluster.num <- min(table(cate.train))
    num.of.cluster <- length(unique(cate.train))
    if (min.cluster.num<5|num.of.cluster==1) {
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
                               sampsize = rep(min.cluster.num,num.of.cluster),
      )
      pred <- predict(object = model.cv,newdata = data.test)
      error.cv <- c(error.cv,mean(pred!=cate.test))
      data.cv <- rbind(data.train,data.test)
      cate.cv <- c(cate.train%>%as.character(),pred%>%as.character())%>%as.factor()
      set <- sample(1:nrow(data.cv),nrow(data.cv),replace = F)
      data.cv <- data.cv[set,]
      cate.cv <- cate.cv[set]
      if(e){
        if(error.cv[run]<e){
          break
        }
      }
    }
  }
  close(pb)
  if (return.list) {
    return(list(data=data.cv,cate=cate.cv,error=error.cv))
  }else{
    return(cate.cv)
  }
  
}
