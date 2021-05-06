cca_knn=function(obj1,obj2,SCT=F,k,ucomp,assay1='SCT',assay2='SCT',mode='int',gene.list='F'){
  if(SCT){
    message('SCTing!')
    obj1=SCTransform(obj1,
                     assay = 'RNA',
                     variable.features.n = 5000,
                     do.scale = T,
                     do.center = T)
    obj2=SCTransform(obj2,
                     assay = 'RNA',
                     variable.features.n = 5000,
                     do.scale = T,
                     do.center = T)
  }
  
  #RunCCA
  if (mode=='cca') {
    message('Running CCA!')
    obj1$orig.ident='object1'
    obj2$orig.ident='object2'
    if(length(gene.list)==1){
      gene.list <- union(x = VariableFeatures(object = obj1,assay = assay1), 
                         y = VariableFeatures(object = obj2,assay = assay2))
    }
    test=RunCCA(object1 = obj1,
                object2 = obj2,
                assay1 = assay1,
                assay2 = assay2,
                features = gene.list,
                num.cc = 50,
                rescale = T)
    test=RunUMAP(object = test,reduction = 'cca',dims = 1:30,n.components = ucomp)
  }else{
    message('Running integration!')
    obj1$orig.ident='object1'
    obj2$orig.ident='object2'
    Exi.list <- list(obj1,obj2)
    if(length(gene.list)==1){
      Exi.features <- SelectIntegrationFeatures(object.list = Exi.list, 
                                                nfeatures = 3000,
                                                assay = c(assay1,assay2))
    }else{
      Exi.features <- gene.list
    }
    
    Exi.list <- PrepSCTIntegration(object.list = Exi.list, 
                                   anchor.features = Exi.features,
                                   assay = c(assay1,assay2))
    
    #Identify anchors and integrate with normalization.method = 'SCT'
    Exi.anchors <- FindIntegrationAnchors(object.list = Exi.list, 
                                          assay = c(assay1,assay2),
                                          normalization.method = "SCT", 
                                          anchor.features = Exi.features, 
                                          verbose = T)
    rm(Exi.list)
    gc()
    test <- IntegrateData(anchorset = Exi.anchors, 
                                    new.assay.name = 'int',
                                    normalization.method = "SCT", 
                                    verbose = T)
    #Confirm the integrated data by umap
    test=RunPCA(object = test,npcs = 100)
    test=RunUMAP(test,dims=1:25,reduction = 'pca',n.components = ucomp)
    rm(Exi.anchors)
    gc()
  }
  
  #knncv
  message('Cross validation for training data.')
  library(FNN)
  Idents(test)=test$orig.ident
  cell1=WhichCells(test,idents = 'object1')
  cell2=WhichCells(test,idents = 'object2')
  trdata=test@reductions$umap@cell.embeddings[cell1,]
  trcl=test@meta.data[cell1,'subclass_label']%>%as.character()%>%as.factor()
  trcv=knn.cv(train = trdata,cl = trcl,k = k,prob = T)
  trcv=list(result=trcv,
            dist=attr(trcv,'nn.dist'),
            index=attr(trcv,'nn.index'),
            prob=attr(trcv,'prob'))
  
  filt=function(dist,label,out='result'){
    label=as.character(label)
    ind=which(label!='unclear')
    dist=dist[ind]
    label=label[ind]%>%as.factor()
    
    avg=mean(dist)
    sd=sd(dist)
    ind=which(dist<(avg+sd))
    dist=dist[ind]
    label=label[ind]
    
    tb=table(label)
    fin=names(tb[tb==max(tb)])
    
    prob=vector(mode = 'numeric',length = 9)
    names(prob)=c('L2/3 IT','L4/5 IT','L5 IT','L5 PT','L5/6 NP','L6 CT','L6 IT','L6 IT Car3','L6b')
    for (i in 1:length(levels(label))) {
      prob[levels(label)[i]]=tb[levels(label)[i]]/sum(tb)
    }
    
    if(length(fin)>1){
      if(out=='result')return('unclear')
      if(out=='prob')return(rep(0,length(levels(label))))
    }
    else {
      if(out=='result')return(fin)
      if(out=='prob')return(prob)
      }
  }
  
  library(foreach)
  trcvresult=foreach(i=1:nrow(trcv$dist),.combine = c)%do%{filt(dist = trcv$dist[i,],label = trcl[trcv$index[i,]],out = 'result')}
  test=subset(test,cells = c(cell1[trcvresult==trcl],cell2))
  cell1=WhichCells(test,idents = 'object1')
  cell2=WhichCells(test,idents = 'object2')
  trdata=test@reductions$umap@cell.embeddings[cell1,]
  trcl=test@meta.data[cell1,'subclass_label']%>%as.character()%>%as.factor()
  
  #knn
  message('Runing knn classification!')
  tedata=test@reductions$umap@cell.embeddings[cell2,]
  v2mnknn=knn(train = trdata,test = tedata,cl = trcl,k = k,prob = T)
  v2mnknn=list(result=v2mnknn,
               dist=attr(v2mnknn,'nn.dist'),
               index=attr(v2mnknn,'nn.index'),
               prob=attr(v2mnknn,'prob'))
  nv2mnknn=foreach(i=1:nrow(v2mnknn$dist),.combine = c)%do%{
     filt(dist = v2mnknn$dist[i,],label = trcl[v2mnknn$index[i,]],out = 'result')
  }
  # v2mnknn$result=nv2mnknn
  
  #mnknncv
  message('Cross validation for test data.')
  k2=mean(c(length(nv2mnknn)/100,min(table(nv2mnknn))))%>%floor()
  tecv=knn.cv(train = tedata,cl = nv2mnknn,k = k2,prob = T)
  tecv=list(result=tecv,
            dist=attr(tecv,'nn.dist'),
            index=attr(tecv,'nn.index'),
            prob=attr(tecv,'prob'))
  ntecv=foreach(i=1:nrow(tecv$dist),.combine = rbind)%do%{
    filt(dist = tecv$dist[i,],label = nv2mnknn[tecv$index[i,]],out = 'prob')
    }
  rownames(ntecv)=cell2
  
  #Result
  knn=list(train_cv=trcv,test_knn=v2mnknn,test_cv=tecv)
  result=list(Seurat=test,knn=knn,testprob=ntecv)
  
  return(result)
}
