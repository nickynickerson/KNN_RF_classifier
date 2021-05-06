################################################################################
##This is the workflow for hclust guided binary RF semi-supervised classification.
##First, do the integration. 
##Then, rename the result as 'en.sub'. 
##Finally, go through the rest of the code, including hclust, RF, and query the result.
################################################################################

####Initialize----------------------------
#Require packages
library(Seurat)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(dplyr)
library(randomForest)
library(FNN)
setwd('/active_data/share/Hyk/final.result.codes/')

#Source functions
source('Function_IntegrationBasedCCA.R')
source('Function_myRandomForestCV.R')
source('Function_myrfcv.R')

#Find all daughter nodes for a parent node.
findAllNode <- function(nd=nd,mg=mg){
  if(nd<0){
    return(-nd)
  }
  nodes <- vector(mode = 'numeric')
  nds <- mg[nd,]
  for (i in 1:2) {
    if (nds[i]<0) {
      nodes <- c(-nds[i],nodes)
    }else{
      nodes <- c(findAllNode(nd=nds[i],mg=mg),nodes)
    }
  }
  return(nodes)
}

#Change gene names so they can be used in formula.
TransGeneNames=function(gene){
  t=gene[nchar(gene)>10]
  t=which(gene%in%t)
  t=union(t,grep(pattern = '+-+',x = gene))
  for (j in t) {
    gene[j]=paste0('a',paste(strsplit(gene[j],split = '-')[[1]],collapse = '_'))
    print(gene[j])
  }
  return(gene)
}

#Reverse the changes done to gene names
reverseGeneNames <- function(gene){
  t=gene[nchar(gene)>10]
  t=which(gene%in%t)
  t=union(t,grep(pattern = '+_+',x = gene))
  for (j in t) {
    gene[j]=paste0(paste(strsplit(gene[j],split = '_')[[1]],collapse = '-'))
    gene[j]=substr(gene[j],2,100)
    print(gene[j])
  }
  return(gene)
}

#Predict binary classification RF result with a criteria 'mx'.
myPredict.rf <- function(model,data.blur,gene){
  gene <- TransGeneNames(gene = gene)
  data <- data.blur[,gene] %>% t() %>% scale() %>% t()
  result <- predict(object = model,newdata = data,type = 'prob')
  cate <- vector(mode = 'character',length = nrow(result))
  for (i in 1:nrow(result)) {
    mx <- max(result[i,])
    lb <- c('A','B')
    if(mx>.55){
      cate[i] <-lb[result[i,]==mx] 
    }else cate[i] <- 'C'
  }
  return(cate)
}

#Query the final node each cell land on. 
queryTree <- function(result,tree,node.p){
  if(is.null(node.p)){
    node.p <- nrow(tree)
  }
  if(node.p<0){
    return(node.p)
  }
  label <- result[node.p]
  node.d <- tree[node.p,]%>%sort(decreasing = T)
  lb <- c('A','B')
  if (label=='C') {
    return(node.p)
  }
  return(queryTree(result,tree,node.d[lb==label]))
}

####Run CCA-KNN---------------------------
#Prepare Data
load('en.sub.Rdata')
load('../../diffusionMap/IntegrationAdultInputData.Rdata')
Idents(exi_sc_v2) <- exi_sc_v2$subclass_label
v2.4000 <- subset(exi_sc_v2,cells = WhichCells(exi_sc_v2,downsample = 4000))
v2.6000 <- subset(exi_sc_v2,cells = WhichCells(exi_sc_v2,downsample = 6000))
Idents(exi_sc_v2) <- exi_sc_v2$orig.ident
v2p.4000 <- subset(exi_sc_v2,cells = WhichCells(exi_sc_v2,downsample = 36000))
rm(exi_sc_v3,exi_sc_v2)
gc()

#Run CCA&KNN
v2.4000 <- cca_knn(obj1 = v2.4000, obj2 = en.sub,
                   SCT = F,k = 100,ucomp = 2,
                   assay1 = 'SCT',assay2 = 'integrated')
v2en.4 <- cca_knn(obj1 = v2.4, obj2 = en.sub,
                   SCT = F,k = 100,ucomp = 2,
                   assay1 = 'SCT',assay2 = 'integrated',mode = 'int')
v2.6000 <- cca_knn(obj1 = v2.6000, obj2 = en.sub,
                   SCT = F,k = 100,ucomp = 2,
                   assay1 = 'SCT',assay2 = 'integrated')
v2.p2 <- cca_knn(obj1 = v2p.4000, obj2 = en.sub,
                SCT = F,k = 100,ucomp = 2,
                assay1 = 'SCT',assay2 = 'integrated')
Idents(v2.p$Seurat) <- v2.p$Seurat$orig.ident
#Check result
DimPlot(v2en.4$Seurat,group.by = c('subclass_label','date'),split.by = 'orig.ident')
FeaturePlot(v2en.4$Seurat,features = c('Ptn','Cux1','Fezf2','Bcl11b','Tshz2'),slot = 'data')

DimPlot(v2.p$Seurat,group.by = c('subclass_label','date'),split.by = 'orig.ident')
FeaturePlot(v2.p$Seurat,features = c('Ptn','Cux1','Fezf2','Bcl11b','Tshz2'),slot = 'data')

mean(v2.p$knn$test_cv$result==v2en.4$knn$test_cv$result)

####hcluster------------------------------
##Input Seurat object, out put Seurat object and hclust.
en.sub@meta.data[,'knn1'] <- en.sub$knn
en.sub@meta.data[,'knn'] <- v2.p2$knn$test_cv$result
Idents(en.sub) <- en.sub$knn

en.sub.split <- SplitObject(en.sub)
pivot <- foreach(i=1:length(en.sub.split),.combine = rbind)%do%{
  apply(X = en.sub.split[[i]]@reductions$pca@cell.embeddings[,1:20], 2, mean)
}

rownames(pivot) <- names(en.sub.split)
dist.ma <- dist(x = pivot)
ed.hc <- hclust(dist.ma)
plot(ed.hc)
en.sub$knn <- v2en.4$knn$test_cv$result

####Main loop-----------------------------
##Input Seurat object and hclust, output clarified Seurat object.
i <- 1
obj.list <- list()
obj.list <- c(obj.list,list(en.sub))
blur.cells <- vector(mode = 'character')
model.list <- list()
gene.list <- list()
while (T) {
  if (i>length(obj.list)|i>30) {
    break
  }
  #Find nodes and make labels.
  if (!'nodes'%in%colnames(obj.list[[i]]@meta.data)) {
    obj.list[[i]]$nodes <- nrow(ed.hc$merge)
    nodes.p <- nrow(ed.hc$merge)
  }else{
    nodes.p <- obj.list[[i]]$nodes[1]
  }
  if (nodes.p<0) {
    i <- i+1
    next
  }
  nodes.d <- ed.hc$merge[as.numeric(nodes.p),]
  nodes.d <- sort(nodes.d,decreasing = T)
  a.nodes <- findAllNode(nd = nodes.d[1],mg = ed.hc$merge)
  a.lab <- ed.hc$labels[a.nodes]
  b.nodes <- findAllNode(nd = nodes.d[2],mg = ed.hc$merge)
  b.lab <- ed.hc$labels[b.nodes]
  
  obj.list[[i]]$hclust <- obj.list[[i]]$knn%>%as.character()
  obj.list[[i]]@meta.data[obj.list[[i]]$knn%in%a.lab,"hclust"] <- 'A'
  obj.list[[i]]@meta.data[obj.list[[i]]$knn%in%b.lab,"hclust"] <- 'B'
  obj.list[[i]]$nodes <- obj.list[[i]]$hclust
  obj.list[[i]]@meta.data[obj.list[[i]]$hclust=='A',"nodes"] <- nodes.d[1]
  obj.list[[i]]@meta.data[obj.list[[i]]$hclust=='B',"nodes"] <- nodes.d[2]
  
  #Run rfcv
  if(min(table(obj.list[[i]]$hclust))<5){
    message('Cell group too small!')
    print(a.lab)
    print(b.lab)
  }
  rfcv <- myRandomForestCV(obj = obj.list[[i]],
                           do.SCT = F,
                           clu.lab = 'hclust',
                           only.Seurat = F,
                           assay = 'integrated')
  
  #Split object
  rfcv.cons <- subset(x = rfcv$Seurat,
                      cells = colnames(rfcv$Seurat)[rfcv$Seurat$rf==rfcv$Seurat$hclust])
  Idents(rfcv.cons) <- rfcv.cons$rf
  obj.split <- SplitObject(object = rfcv.cons)
  for (j in 1:length(obj.split)) {
    obj.list <- c(obj.list,list(obj.split[[j]]))
    }
  blur.cells <- c(colnames(rfcv$Seurat)[rfcv$Seurat$rf!=rfcv$Seurat$hclust],blur.cells)
  
  i <- i+1
  model.list[[as.numeric(nodes.p)]] <- rfcv$rfmodel
  gene.list[[as.numeric(nodes.p)]] <- rfcv$gene
}

####Result--------------------------------
#Predict all cells with the final model list.
all.gene <- foreach(i=1:length(model.list),.combine = c)%do%{names(model.list[[i]]$forest$ncat)}
all.gene <- all.gene%>%unique()%>%reverseGeneNames()
data.blur <- en.sub@assays$integrated@scale.data%>%as.matrix()%>%t()
colnames(data.blur) <- TransGeneNames(colnames(data.blur))
blur.cell.idents <- foreach(i=1:length(model.list),.combine = cbind)%do%{
  tryCatch(myPredict.rf(model.list[[i]],data.blur,gene.list[[i]]),error=function(e){print(model.list[[i]])})
}
rownames(blur.cell.idents) <- rownames(data.blur)
colnames(blur.cell.idents) <- 1:8

#query result
blur.nodes <- foreach(i=1:nrow(blur.cell.idents),.combine = c)%do%{queryTree(result = blur.cell.idents[i,],tree = ed.hc$merge,node.p = 8)}
blur.label <- vector(mode = 'character',length = length(blur.nodes))
for (i in 1:length(blur.nodes)) {
  if (blur.nodes[i]<0) {
    blur.label[i] <- ed.hc$labels[-blur.nodes[i]]
  }else{
    blur.label[i] <- 'unclear'
  }
}
obj.list[[1]]$rfcv <- blur.label
Idents(en.sub) <- en.sub$rfcv

FeaturePlot(en.sub,features = c('Nxph4'),reduction = 'umap')
DimPlot(en.sub,reduction = 'umap',group.by = c('rfcv','knn'))
