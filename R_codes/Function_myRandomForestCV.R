###############################################################
##Use random forest to clarify clusters
###############################################################
myRandomForestCV <- function(obj,
                             do.SCT=T,
                             clu.lab,
                             gene.ini=T,
                             only.Seurat=T,
                             assay='SCT'){
  ##Require packages & functions----------
  #packages
  library(Seurat)
  library(parallel)
  library(foreach)
  library(iterators)
  library(doParallel)
  library(dplyr)
  library(randomForest)
  library(FNN)
  #library(caret)
  
  #Functions
  pickcell=function(data=data,p=p){
    prob=max(data)
    if(prob>=p) return(1)
    else return(0)
  }
  
  TransGeneNames=function(gene=gene){
    t=gene[nchar(gene)>10]
    t=which(gene%in%t)
    t=union(t,grep(pattern = '+-+',x = gene))
    for (j in t) {
      gene[j]=paste0('a',paste(strsplit(gene[j],split = '-')[[1]],collapse = '_'))
      print(gene[j])
    }
    return(gene)
  }
  
  pickRFCV <- function(rfcv,cate,mode='avg'){
    if(mode=='avg'){
      pred <- rfcv$predicted
      error.rate <- vector(mode = 'numeric',length = length(pred))
      for (i in 1:length(pred)) {
        tb <- table(cate,pred[[i]])
        er <- foreach(j=1:nrow(tb),.combine = c)%do%{
          tb[j,j]/sum(tb[j,])
        }
        error.rate[i] <- mean(er)
      }
      return(names(pred)[which(error.rate==min(error.rate))])
    }else{
      return(rfcv$n.var[which(rfcv$error.cv==min(rfcv$error.cv))])
    }
  }
  
  ##Process data--------
  #Label
  if (is.na(clu.lab)) {
    stop('Must have a label!')
  } else{
    Idents(obj) <- obj@meta.data[,clu.lab]
  }
  
  #SCT
  if (do.SCT) {
    obj <- SCTransform(object = obj,
                       assay = 'RNA',
                       do.scale = T,
                       do.center = T,
                       seed.use = 1,
                       )
  }
  
  #gene
  message('Finding genes!')
  if (gene.ini) {
    markers <- FindAllMarkers(object = obj,assay = assay,only.pos = T)
    if (nrow(markers)<20) {
      warning('Not enough feature! Set thred.hold = .15!')
      markers <- FindAllMarkers(object = obj,assay = assay,only.pos = T,logfc.threshold = .15)
      if(nrow(markers<20)){
        warning('These two group don\'t have enough features.')
        print(markers$gene)
      }
    }
    if (nrow(markers)>200) {
      markers <- markers%>%group_by(cluster)%>%top_n(n = 100,wt = avg_logFC)
    }
    gene <- markers$gene%>%unique()
  } else gene <- gene.ini
  print(length(gene))
  if (assay=='SCT') {
    obj <- GetResidual(object = obj,features = gene,assay = assay)
  }
  
  #Get RF inputdata
  message('Prepare input data!')
  DefaultAssay(obj) <- assay
  data <- FetchData(object = obj,vars = gene,slot = 'data')%>%as.matrix()
  colnames(data)=TransGeneNames(gene = colnames(data))
  cate <- Idents(obj)%>%as.character()%>%as.factor()
  names(cate)=rownames(obj@meta.data)
  
  ##Build primary RF model & pick out features----------
  #Primary RF model
  message('Build First model.')
  model.st <- randomForest(x = data,
                           y = cate,
                           importance=TRUE,
                           proximity=F,
                           replace=T,
  )
  imp <- importance(model.st,type = '1')
  
  #Pick out features
  message('Select features.')
  model.st.cv <- rfcv(trainx = data,
                      trainy = cate,
                      ntree=300,
                      cv.fold = 10,
                      scale = 'log',
                      step = 0.75,
                      recursive = F)
  num <- pickRFCV(model.st.cv,cate = cate)
  num <- max(num,min(10,length(gene)))
  message(num)
  ind <- order(-imp)[1:num]
  data <- data[,rownames(imp)[ind]]
  
  ##Rebuild RF model pick out training cells-------
  #Rebuild RF model
  message('Build Second model.')
  model.nd <- randomForest(x = data,
                           y = cate,
                           ntree=500,
                           importance=TRUE,
                           proximity=F,
                           replace=T,
                           )
  
  #Pick out votes>=0.6 as training set.
  votes <- model.nd$votes%>%as.data.frame()
  a.cells <- votes%>%filter(A>=0.6)%>%rownames()
  b.cells <- votes%>%filter(B>=0.6)%>%rownames()
  a.data <- data[a.cells,]
  b.data <- data[b.cells,]
  cate.cv <- c(rep('A',length(a.cells)),rep('B',length(b.cells)))%>%as.factor()
  data.cv <- rbind(a.data,b.data)
  
  ##100 runs of random forest 10-fold cross validation---------
  message('100 runs of CV.')
  error.cv <- vector(mode = 'numeric')
  pb <- txtProgressBar(style=3)
  for (run in 1:100) {
    setTxtProgressBar(pb, run/100)
    sig <- -1
    n <- floor(nrow(data.cv)/10)
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
      )
      pred <- predict(object = model.cv,newdata = data.test)
      error.cv <- c(error.cv,mean(pred!=cate.test))
      data.cv <- rbind(data.train,data.test)
      cate.cv <- c(cate.train%>%as.character(),pred%>%as.character())%>%as.factor()
      set <- sample(1:nrow(data.cv),nrow(data.cv),replace = F)
      data.cv <- data.cv[set,]
      cate.cv <- cate.cv[set]
      
    }
  }
  close(pb)
  
  ##Build Final RF model and predict the rest---------
  message('Build the Final model.')
  if (cate.cv %>% unique() %>% length()==1) {
    model.fn <- cate.cv %>% unique()
    cate.fn <- rep(cate.cv %>% unique(),ncol(obj))
    error.cv <- NA
  }else{
    model.fn <- randomForest(x = data.cv,
                             y = cate.cv,
                             proximity = F,
                             importance = T
    )
    pred.fn <- predict(object = model.fn,
                       newdata = data,
                       type = 'prob')
    cate.fin <- vector(mode = 'character',length = nrow(pred.fn))
    for (i in 1:nrow(pred.fn)) {
      mx <- max(pred.fn[i,])
      lb <- c('A','B')
      if(mx>.50){
        cate.fin[i] <-lb[pred.fn[i,]==mx] 
      }else cate.fin[i] <- 'C'
    }
  }
  
  
  ##Output the result------------
  message('Out put the result.')
  obj$rf <- cate.fin
  if (only.Seurat) {
    return(obj)
  }else{
    # rf.list=list(First=model.st,
    #              Firstcv=model.st.cv,
    #              Second=model.nd,
    #              Final=model.fn)
    return(list(Seurat=obj,error.cv=error.cv,rfmodel=model.fn,gene=gene))
  }
}
