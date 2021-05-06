################################################################################
##This is the workflow for RF semi-supervised multi-classification.
##This includes segmentation, and classification. 
################################################################################

# Initialization ----------------------------------------------------------

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(readr)
library(pheatmap)
library(randomForest)
library(diffusionMap)
library(princurve)
setwd('/active_data/share/Hyk/final.result.codes/')
options(future.globals.maxSize = 400 * 1024^3)

####Functions####
source('Function_IntegrationBasedCCA.R')
source('Function_myRandomForestCV.R')
source('Function_myrfcv.R')
source('Function_rfPermutation.R')

pickcell=function(data=data,p=p){
  prob=max(data)
  if(prob>=p) return(1)
  else return(0)
}

#Change gene names so they can be used in formula.
TransGeneNames=function(gene=gene){
  t=gene[nchar(gene)>10]
  t=which(gene%in%t)
  t=union(t,grep(pattern = '+-+',x = gene))
  #Add 'a' in the front of long gene names and substitute '_' for '-'.
  for (j in t) {
    gene[j]=paste0('a',paste(strsplit(gene[j],split = '-')[[1]],collapse = '_'))
    print(gene[j])
  }
  return(gene)
}

#Calculate the center of a dataset. Didn't work out well. Use average instead.
centShift <- function(data=data,h=h,e=1e-4,maxit=Inf){
  pt=apply(data, 2, sum)/nrow(data)
  pt=pt%>%as.matrix()%>%t()
  
  ##Functions
  #dist of one point and a matrix
  myDist <- function(point=point,da=da){
    library(foreach)
    return(foreach(i = 1:nrow(da),.combine = c)%do%{
      dist(rbind(pt,da[i,]))
    }
    )
  }
  
  #Calculate shift vector
  shift <- function(fpoint=pt,fdata=data,h=h){
    d <- myDist(point = fpoint,da = fdata)
    knn.ind <- which(d<h)
    knn.mat <- fdata[knn.ind,]
    fk=nrow(knn.mat)
    vec <- apply(X = knn.mat,MARGIN = 2,FUN = sum)
    vec <- vec-fk*fpoint
    E <- norm(vec,type = '2')
    npoint <- vec+fpoint
    return(list(pt=npoint,error=E))
  }
  
  #Iterate
  h <- 2*dist(data)%>%mean()
  itr <- 1
  E <- 1
  while(itr<maxit&&E>e){
    itr <- itr+1
    temp <- shift(fpoint = pt,fdata = data,h = h)
    E <- temp$error
    pt <- temp$pt
  }
  
  #Output
  return(pt)
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

#Pick out high quality cells based on votes. This function doesn't return HQ cells directly. It returns figures.
pickRFHqCell <- function(rfmodel){
  votes<- rfmodel$votes
  dif <- vector(mode = 'numeric',length = nrow(votes))
  st <- vector(mode = 'numeric',length = nrow(votes))
  for (i in 1:nrow(votes)) {
    st[i] <- sort(votes[i,],decreasing = T)[1]
    nd <- sort(votes[i,],decreasing = T)[2]
    dif[i] <- st[i]-nd
  }
  dif <- dif %>% as.matrix() %>% as.data.frame()
  st <- st%>% as.matrix() %>% as.data.frame()
  gg<- ggplot(dif,aes(V1))+geom_histogram(bins = 60)
  gg2 <- ggplot(st,aes(V1))+geom_histogram(bins = 60)
  return(list(st=st,dif=dif,gg.dif=gg,gg.st=gg2))
}

####For real####
#seperate the neurons
load('/active_data/share/yuxy/obj.integrated_annotated_20210303.RData')
load('IntegrationAdultInputData.Rdata')

DimPlot(obj.integrated,group.by = 'integrated_snn_res.1.8',label = T,repel = T)+NoLegend()
Idents(obj.integrated) <- obj.integrated$celltype
lineage.en <- subset(x = obj.integrated,idents = c('RGP','IP','CP EN','SP EN'))
lineage.en <- FindNeighbors(lineage.en,reduction = 'umap',graph.name = 'umap.snn',dims = 1:2)
lineage.en <- FindClusters(lineage.en,graph.name = 'umap.snn',resolution = 2.8)

Idents(lineage.en) <- lineage.en$umap.snn_res.2.8
stage1 <- subset(lineage.en,idents = c(11,13,90,0,44,19,77,84,1,21,22,93,53,33,64,83,82,18,8,66))
stage2 <- subset(lineage.en,idents = c(25,16,71,46,6,88,78,67,17,62,63,55,48,68))
stage3 <- subset(lineage.en,idents = c(65,12,94,87,42,73,86,56,54,31))
cell.stage1 <- WhichCells(lineage.en,idents = c(11,13,90,0,44,19,77,84,1,21,22,93,53,33,64,83,82,18,8,66))
cell.stage2 <- WhichCells(lineage.en,idents = c(25,16,71,46,6,88,78,67,17,62,63,55,48,68))
cell.stage3 <- WhichCells(lineage.en,idents = c(65,12,94,87,42,73,86,56,54,31))
DimPlot(lineage.en,cells.highlight = list(cell.stage1,cell.stage2,cell.stage3),sizes.highlight = .1,cols.highlight = c('steelblue','green','red'))

lineage.neuron <- list(stage1=stage1,stage2=stage2,stage3=stage3)
cell.neuron <- list(stage1=cell.stage1,stage2=cell.stage2,stage3=cell.stage3)
save(lineage.en,file = 'lineage.en.Rdata')
save(lineage.neuron,cell.neuron,file = 'lineage.neuron.Rdata')

#individual umaps
for (i in 1:3) {
  lineage.neuron[[i]] <- RunPCA(lineage.neuron[[i]],
                                assay = 'integrated',
                                npcs = 80,
                                reduction.name = 'ipca')
  lineage.neuron[[i]] <- RunUMAP(lineage.neuron[[i]],
                                 reduction = 'ipca',
                                 dims = 1:20,
                                 reduction.name = 'iumap')
}
DimPlot(lineage.neuron[[2]],reduction = 'umap',group.by = 'umsp.snn_res.1.8',label = T,repel = T)
FeaturePlot(lineage.neuron[[2]],features = c('Pax6','Eomes','Neurog2'),ncol = 3)

#Integration
lineage.neuron$stage1$orig.ident <- 'stage1'
lineage.neuron$stage2$orig.ident <- 'stage2'
lineage.neuron$stage3$orig.ident <- 'stage3'
v2$orig.ident <- 'v2'
v2.stage1 <- cca_knn(obj1 = exi_sc_v2,obj2 = lineage.neuron[[1]],
                     k = 100,assay1 = 'SCT',assay2 = 'integrated',
                     SCT = F,ucomp = 2,mode = 'int')
v2.stage1$Seurat$knn <- v2.stage1$knn$test_knn$result
DimPlot(v2.stage1$Seurat,
        group.by = 'orig.ident',
        #split.by = 'orig.ident',
        #cells = WhichCells(v2.stage1$Seurat,idents = 'object2'),
        #label = T,repel = T, 
)
Idents(v2.stage1$Seurat) <- v2.stage1$Seurat$orig.ident
em <- subset(v2.stage1$Seurat,idents = 'object2')
em <- FindNeighbors(em,reduction = 'umap',dims = 1:2,graph.name = 'umap.snn')
em <- FindClusters(em,graph.name = 'umap.snn',resolution = 1.5)
em$knn <- v2.stage1$knn$test_cv$result
knn.hq.cell <- WhichCells(em,idents = c('12','35','11','16','41'),invert = T)
knn.lq.cell <- WhichCells(em,idents = c('12','35','11','16','41'),invert = F)
DimPlot(em,label=T,repel = T,cells = knn.hq.cell)
Idents(em) <- em$umap.snn_res.1
lineage.neuron[[1]]$knn <- v2.stage1$knn$test_cv$result
DimPlot(lineage.neuron[[1]],
        group.by = 'knn',
        ncol = 1,reduction = 'umap',
        repel = T,label = T,cells = knn.lq.cell)

####rfcv####

#Select feature
stage1.knn.hq <- subset(lineage.neuron$stage1,cells = knn.hq.cell)
Idents(stage1.knn.hq) <- stage1.knn.hq$knn
marker.knn.hq <- FindAllMarkers(object = stage1.knn.hq,assay = 'integrated',
                                only.pos = F)
marker.knn.hq.sct <- FindAllMarkers(object = stage1.knn.hq,assay = 'SCT',
                                    only.pos = T)
marker.knn.hq <- marker.knn.hq %>% filter(p_val_adj<.005)
gene <- marker.knn.hq$gene %>% unique()
data.s1 <- stage1.knn.hq@assays$integrated@data[gene,] %>% as.matrix() %>% t() %>% as.data.frame()
cate.s1 <- stage1.knn.hq$knn %>% as.character() %>% as.factor()
colnames(data.s1) <- TransGeneNames(colnames(data.s1))
rf.s1 <- randomForest(x = data.s1,y = cate.s1,
                      ntree = 400,mtry = 100,
                      replace = T,importance = T)
imp <- rf.s1$importance
rfcv.s1 <- rfcv(trainx = data.s1,cate.s1,cv.fold=10,step=.7,ntree=300)
gene.imp <- imp %>% as.data.frame %>% top_n(n = 120,wt = MeanDecreaseGini) %>% rownames()

#Build the refined model.
data.s1.f <- data.s1[,gene.imp]
rf.s1.f <- randomForest(x = data.s1.f,y = cate.s1,
                        ntree = 400,
                        importance = T,replace = T)
stage1.knn.hq$rf <- rf.s1.f$predicted
DimPlot(stage1.knn.hq,group.by = c('knn','rf'))

#Pickout HQ cells
s1.hq <- pickRFHqCell(rf.s1.f)
ind2 <- union(which(s1.hq$dif>.25),which(s1.hq$st>.5))
cell.rf.hq <- rownames(votes)[ind2]
stage1.knn.hq$rf <- rf.s1.f$predicted
DimPlot(stage1.knn.hq,cells = cell.rf.hq,group.by = 'rf')
FeaturePlot(stage1.knn.hq,cells = cell.rf.hq,features = 'Fezf2')

#Build final model
data.fin <- data.s1.f[cell.rf.hq,]
cate.fin <- stage1.knn.hq@meta.data[ind2,"rf"]
cate.fin <- cate.fin %>% as.character() %>% as.factor()
rf.fin <- randomForest(x = data.fin,y = cate.fin,
                       ntree = 600,
                       replace = T,importance = T,
                       sampsize = rep(10,6),
)
gene.imp.reverse <- reverseGeneNames(gene.imp)
data.s1 <- lineage.neuron$stage1@assays$integrated@data[gene.imp.reverse,] %>% 
  as.matrix() %>% 
  t() %>% 
  as.data.frame()
colnames(data.s1) <- TransGeneNames(colnames(data.s1))
result <- predict(rf.fin,newdata = data.s1)
lineage.neuron$stage1$rf <- result
DimPlot(lineage.neuron$stage1,group.by = 'rf')

#Do permutation
rf.fin.hq <- randomForest(x = data.s1,y = result,
                          ntree = 500,proximity = F,
                          sampsize = rep(200,6))


save(result,data.s1,file = 'permutation.Rdata')
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
cl.cores <- detectCores()
cl <-makeCluster(getOption("cl.cores", 10),type="FORK")
registerDoParallel(cl = cl)
per.result <- foreach(i=1:10,.combine = 'cbind')%dopar%{
  rfPermutation(data.cv = data.s2,cate.cv = result.s2,
                fold = 10,return.list = F)
}

votes <- data.frame(row.names = rownames(per.result))
cate <- per.result[,1] %>% as.character() %>%  unique()
for (i in 1:6) {
  votes[,i] <- foreach(j=1:nrow(votes),.combine = c)%do%{
    mean(per.result[j,]==cate[i])
  }
}
vmax <- vector(mode = 'numeric',length = nrow(votes))
dif <- vector(mode = 'numeric',length = nrow(votes))
for (i in 1:nrow(votes)) {
  vmax[i] <- sort(votes[i,],decreasing = T)[1,1]
  dif <- sort(votes[i,],decreasing = T)[1,1]-sort(votes[i,],decreasing = T)[2,1]
}
table(vmax)
ind4 <- which(vmax==1)
cell.stage1.fin <- rownames(votes)[ind4]
lineage.neuron$stage1$rf <- result
stage1.fin <- subset(lineage.neuron$stage1,cells = cell.stage1.fin)
DimPlot(stage1.fin,group.by = 'rf')
save(stage1.fin,file = 'stage1.fin.Rdata')

Idents(stage1.fin) <- stage1.fin$rf
markers.s1.int <- FindAllMarkers(object = stage1.fin,assay = 'integrated',
                                 only.pos = T)
markers.s1.sct <- FindAllMarkers(object = stage1.fin,assay = 'SCT',
                                 only.pos = F)

####predict stage2####

#int
gene.s2 <- FindVariableFeatures(lineage.neuron$stage2,assay = 'integrated',nfeatures = 2000)
gene.s2 <- gene.s2@assays$integrated@var.features
Idents(stage1.fin) <- stage1.fin$rf
markers.s1.int <- markers.s1.int %>% filter(p_val_adj<.01)
gene.int <- markers.s1.int$gene %>% unique()
stage1.fin$subclass_label <- stage1.fin$rf #To use function cca_knn object1 must have subclass_label as true label.
stage12.int <- cca_knn(obj1 = stage1.fin,
                       obj2 = lineage.neuron$stage2,
                       SCT = F,k = 15,ucomp = 2,
                       assay1 = 'integrated',assay2 = 'integrated',
                       mode = 'int',gene.list = intersect(gene.s2,gene.int))
Idents(stage12.int$Seurat) <- stage12.int$Seurat$orig.ident
em <- subset(stage12.int$Seurat,idents = 'object2')
em <- FindNeighbors(em,reduction = 'umap',dims = 1:2,graph.name = 'umap.snn')
em <- FindClusters(em,graph.name = 'umap.snn',resolution = 1.5)
em$knn <- v2.stage1$knn$test_cv$result
knn.hq.cell <- WhichCells(em,idents = c('15','12','5','27'),invert = T)
knn.lq.cell <- WhichCells(em,idents = c('15','12','5','27'),invert = F)
DimPlot(em,label=T,repel = T,cells = knn.hq.cell)

lineage.neuron$stage2$knn <- stage12.int$knn$test_knn$result
DimPlot(lineage.neuron$stage2,group.by = 'knn')
stage2.knn.hq <- subset(lineage.neuron$stage2,cells = knn.hq.cell)
Idents(stage2.knn.hq) <- stage2.knn.hq$knn
marker.s2 <- FindAllMarkers(stage2.knn.hq,only.pos = T,assay = 'integrated')
marker.s2 <- marker.s2 %>% filter(p_val_adj<.01)
gene.s2 <- marker.s2$gene %>% unique()
data.s2 <- stage2.knn.hq@assays$integrated@data %>% as.matrix() %>% t
data.s2 <- data.s2[,gene.s2]
cate.s2 <- stage2.knn.hq$knn %>% as.character() %>% as.factor()
rf.s2 <- randomForest(data.s2,cate.s2,
                      ntree = 800,importance = T,
                      proximity = F,mtry = 100,
                      sampsize = rep(20,6))
imp.s2 <- rf.s2$importance
gene.imp <- imp.s2 %>% as.data.frame %>% top_n(n = 200,wt = MeanDecreaseGini) %>% rownames()
data.s2 <- data.s2[,gene.imp]
rf.s2.f <- randomForest(data.s2,cate.s2,
                        ntree = 800,importance = T,
                        proximity = F,
                        sampsize = rep(20,6))
stage2.knn.hq$rf <- rf.s2.f$predicted

DimPlot(stage2.knn.hq,group.by = 'rf')
rf.s2.hq <- pickRFHqCell(rf.s2.f)
ind21 <- union(which(rf.s2.hq$st[,1]>0.45),which(rf.s2.hq$dif[,1]>0.18))
cell.s2.rf.hq <- rf.s2.f$votes[ind21,] %>% rownames()
stage2.rf.hq <- subset(stage2.knn.hq,cells = cell.s2.rf.hq)
DimPlot(stage2.rf.hq,group.by = 'rf')
stage2.rf.hq@meta.data[stage2.rf.hq$rf=='L5 IT','rf'] <- 'L4/5 IT'
data.s2.fin <- stage2.rf.hq@assays$integrated@data %>% as.matrix() %>% t
data.s2.fin <- data.s2.fin[,gene.imp]
cate.s2.fin <- stage2.rf.hq$rf %>% as.character() %>% as.factor()
rf.v2.fin <- randomForest(data.s2.fin,cate.s2.fin,
                          ntree = 800,importance = T,
                          proximity = F,
                          sampsize = rep(20,5))
data.s2 <- lineage.neuron$stage2@assays$integrated@data %>% as.matrix() %>% t()
data.s2 <- data.s2[,gene.imp]
result.s2 <- predict(rf.v2.fin,data.s2)
lineage.neuron$stage2$RF <- result.s2
DimPlot(lineage.neuron$stage2,group.by = 'RF')
save(lineage.neuron,file = 'lineage.neuron.Rdata')
per.result.s2 <- foreach(i=1:10,.combine = 'cbind')%dopar%{
  rfPermutation(data.cv = data.s2,cate.cv = result.s2,
                fold = 10,return.list = F)
}
save(per.result.s2,file = 'per.result.s2.Rdata')
vmax <- vector(mode = 'numeric',length = nrow(votes))
dif <- vector(mode = 'numeric',length = nrow(votes))
for (i in 1:nrow(votes)) {
  vmax[i] <- sort(votes[i,],decreasing = T)[1,1]
  dif <- sort(votes[i,],decreasing = T)[1,1]-sort(votes[i,],decreasing = T)[2,1]
}
table(vmax)
ind4 <- which(vmax==1)
cell.stage2.fin <- rownames(votes)[ind4]
stage2.fin <- subset(lineage.neuron$stage2,cells = cell.stage2.fin)
DimPlot(stage2.fin,group.by = 'RF')
Idents(stage2.fin) <- stage2.fin$RF
marker.s2.fin <- FindAllMarkers(stage2.fin,only.pos = T,assay = 'integrated')
marker.s2.fin.sct <- FindAllMarkers(stage2.fin,only.pos = T,assay = 'SCT')
marker.s2.fin <- marker.s2.fin %>% filter(p_val_adj<.005)
data.s2 <- stage2.fin@assays$integrated@data %>% as.matrix() %>% t
data.s2 <- data.s2[,marker.s2.fin$gene %>% unique()]
cate.s2 <- stage2.fin$RF %>% as.character() %>% as.factor()
rf.s2 <- randomForest(data.s2,cate.s2,
                      importance = T,proximity = F,
                      ntree = 1000,mtry = 200,
                      sampsize = rep(10,5))
imp.s2 <- rf.s2$importance
gene.s2 <- imp.s2 %>% as.data.frame() %>% top_n(n = 200,wt = MeanDecreaseGini) %>% rownames()
data.s2 <- data.s2[,gene.s2]
rf.s2 <- randomForest(data.s2,cate.s2,
                      importance = T,proximity = F,
                      ntree = 1000,
                      sampsize = rep(10,5))
stage2.fin$rf <- rf.s2$predicted
DimPlot(stage2.fin,group.by = 'rf')
save(stage2.fin,file = 'stage2.fin.Rdata')
save(marker.s2.fin,marker.s2.fin.sct,file = 'stage2.merkers.Rdata')

####Stage2->Stage3####
gene.s3 <- FindVariableFeatures(object = lineage.neuron$stage3,assay = 'SCT',nfeatures = 3000)
gene.s3 <- gene.s3@assays$SCT@var.features
stage2.fin$subclass_label <- stage2.fin$rf
stage23.int <- cca_knn(obj1 = stage2.fin,obj2 = lineage.neuron$stage3,
                       assay1 = 'SCT',assay2 = 'SCT',k = 20,
                       ucomp = 2,SCT = F,mode = 'int',gene.list = gene)
DimPlot(stage23.int$Seurat,group.by = 'umap.snn_res.2.8',split.by = 'orig.ident')
lineage.neuron$stage3$knn <- stage23.int$knn$test_cv$result
DimPlot(lineage.neuron$stage3,split.by = 'knn')
Idents(stage23.int$Seurat) <- stage23.int$Seurat$orig.ident
em <- subset(stage23.int$Seurat,idents = 'object2')
em$knn <- stage23.int$knn$test_cv$result
DimPlot(em,group.by = 'knn')
em <- FindNeighbors(em,reduction = 'umap',dims = 1:2,graph.name = 'umap.snn')
em <- FindClusters(em,graph.name = 'umap.snn',resolution = 2.5)
votes <- stage23.int$testprob
vmax <- vector(mode = 'numeric',length = nrow(votes))
dif <- vector(mode = 'numeric',length = nrow(votes))
for (i in 1:nrow(votes)) {
  vmax[i] <- sort(votes[i,],decreasing = T)[1]
  dif[i] <- sort(votes[i,],decreasing = T)[1]-sort(votes[i,],decreasing = T)[2]
}
dif <- dif %>% as.matrix() %>% as.data.frame()
vmax <- vmax%>% as.matrix() %>% as.data.frame()
gg<- ggplot(dif,aes(V1))+geom_histogram(bins = 60)
gg2 <- ggplot(vmax,aes(V1))+geom_histogram(bins = 60)
ind4 <- which(dif>=0.95)
knn.hq.cell <- rownames(votes)[ind4]
stage3.knn.hq <- subset(lineage.neuron$stage3,cells = knn.hq.cell)
DimPlot(stage3.knn.hq,split.by = 'knn')

#rfcv
data.s3 <- lineage.neuron$stage3@assays$integrated@data %>% as.matrix() %>% t()
cate.s3 <- lineage.neuron$stage3$knn %>% as.character() %>% as.factor()
Idents(stage3.knn.hq) <- stage3.knn.hq$knn
marker.s3.int <- FindAllMarkers(stage3.knn.hq,only.pos = T,
                                assay = 'integrated',)
marker.s3.int <- FindAllMarkers(stage3.knn.hq,only.pos = T,
                                assay = 'integrated',)
marker.s3.int <- marker.s3.int %>% filter(p_val_adj<.01)
data.s3 <- data.s3[,marker.s3.int$gene %>% unique()]
rf.s3 <- randomForest(data.s3,cate.s3,
                      ntree = 500,mtry = 100,
                      importance = T,proximity = F,
                      nodesize = 5)
imp.s3 <- rf.s3$importance
save(data.s3,cate.s3,file = 'rfcv.stage3.Rdata')
load('rfcv.stage3.Rdata')
gene.imp <- imp.s3 %>% as.data.frame %>% top_n(n = 85,wt = MeanDecreaseGini) %>% rownames()
data.s3 <- data.s3[,gene.imp]
rf.s3 <- randomForest(data.s3,cate.s3,
                      ntree = 500,
                      importance = T,proximity = F,
)
votes <- rf.s3$votes
vmax <- vector(mode = 'numeric',length = nrow(votes))
dif <- vector(mode = 'numeric',length = nrow(votes))
for (i in 1:nrow(votes)) {
  vmax[i] <- sort(votes[i,],decreasing = T)[1]
  dif[i] <- sort(votes[i,],decreasing = T)[1]-sort(votes[i,],decreasing = T)[2]
}
dif <- dif %>% as.matrix() %>% as.data.frame()
vmax <- vmax%>% as.matrix() %>% as.data.frame()
gg<- ggplot(dif,aes(V1))+geom_histogram(bins = 60)
gg2 <- ggplot(vmax,aes(V1))+geom_histogram(bins = 60)
ind31 <- which(vmax[,1]>.6)
cell.s3.rf.hq <- rownames(votes)[ind31]
stage3.knn.hq$rf <- rf.s3$predicted
DimPlot(stage3.knn.hq,cells = cell.s3.rf.hq,group.by = 'rf')
stage3.rf.hq <- subset(stage3.knn.hq,cells = cell.s3.rf.hq)
data.s3.fin <- stage3.rf.hq@assays$integrated@data[gene.imp,] %>% as.matrix() %>% t
cate.s3.fin <- stage3.rf.hq$rf %>% as.character() %>% as.factor()
rf.s3.fin <- randomForest(data.s3.fin,cate.s3.fin,
                          ntree = 1000,
                          importance = T,proximity = F,
                          sampsize = rep(80,3),nodesize = 5)
data.s3 <- lineage.neuron$stage3@assays$integrated@data[gene.imp,] %>% as.matrix() %>% t
result.s3 <- predict(rf.s3.fin,newdata = data.s3)
lineage.neuron$stage3$rf <- result.s3
save(data.s3,result.s3,file = 'permutation.Rdata')

votes <- data.frame(row.names = rownames(per.result))
cate <- per.result[,1] %>% as.character() %>%  unique()
for (i in 1:3) {
  votes[,i] <- foreach(j=1:nrow(votes),.combine = c)%do%{
    mean(per.result[j,]==cate[i])
  }
}
vmax <- vector(mode = 'numeric',length = nrow(votes))
dif <- vector(mode = 'numeric',length = nrow(votes))
for (i in 1:nrow(votes)) {
  vmax[i] <- sort(votes[i,],decreasing = T)[1,1]
  dif <- sort(votes[i,],decreasing = T)[1,1]-sort(votes[i,],decreasing = T)[1,2]
}
dif <- dif %>% as.matrix() %>% as.data.frame()
vmax <- vmax%>% as.matrix() %>% as.data.frame()
gg<- ggplot(dif,aes(V1))+geom_histogram(bins = 60)
gg2 <- ggplot(vmax,aes(V1))+geom_histogram(bins = 60)
ind32 <- which(vmax[,1]==1)
cell.s3.fin <- rownames(votes)[ind32]
stage3.fin <- subset(lineage.neuron$stage3,cells = cell.s3.fin)
DimPlot(stage3.fin,group.by = 'stage')
save(stage3.fin,file = 'stage3.fin.Rdata')

#Add stage2.fin's L6b into stage3 to create stage4.
Idents(stage2.fin) <- stage2.fin$rf
L6b <- subset(stage2.fin,idents = 'L6b')
DimPlot(L6b)
cell.L6b <- colnames(L6b)
cell.stage3.fin <- colnames(stage3.fin)
stage4 <- subset(lineage.en,cells = c(cell.stage3.fin,cell.L6b))
label.rf.s3 <- stage3.fin$rf %>% as.character()
label.rf.l6b <- L6b$rf %>% as.character()
stage4$rf <- c(label.rf.s3,label.rf.l6b)
DimPlot(stage4,group.by = 'rf')
save(stage4,file = 'stage4.Rdata')
Idents(stage4) <- stage4$rf
markers.s4.int <- FindAllMarkers(stage4,assay = 'integrated',only.pos = T)
markers.s4.int <- markers.s4.int %>% filter(p_val_adj<.005)
markers.s4.sct <- FindAllMarkers(stage4,assay = 'SCT',only.pos = T)


# Select IP ---------------------------------------------------------------

lineage.en <- CellCycleScoring(lineage.en,
                               cc.genes.updated.2019$s.genes,
                               cc.genes.updated.2019$g2m.genes)
Idents(lineage.en) <- lineage.en$umap.snn_res.2.8
IP <- subset(lineage.en,idents = c('57','58','38','41','35','47','7','10','80','95','75','28','37','26','34','91'))
Idents(IP) <- IP$Phase
IP.g1 <- subset(IP,idents = 'G1')
DimPlot(IP.g1)
IP.cc <- ScaleData(IP,vars.to.regress = c('S.Score','G2M.Score'),
                   assay = 'integrated',
                   features = IP@assays$integrated@var.features)
IP.cc <- RunPCA(IP.cc,assay = 'integrated',reduction.name = 'cc.pca')
IP.cc <- RunUMAP(IP.cc,reduction = 'cc.pca',dims = 1:28,reduction.name = 'cc.umap')
DimPlot(IP.g1,group.by = c('celltype','Phase','stage'))
IP.g1 <- FindNeighbors(IP.g1,reduction = 'umap',graph.name = 'snn',dims = 1:2)
IP.g1 <- FindClusters(IP.g1,graph.name = 'snn',resolution = .1)
ip.g1.1 <- subset(IP.g1,idents = c('1','2'))
ip.g1.1 <- RunPCA(ip.g1.1,assay = 'integrated')
ip.g1.1 <- RunUMAP(ip.g1.1,reduction = 'pca',reduction.name = 'iumap',dims = 1:30)

# IP+stage4 ---------------------------------------------------------------

#Select features
gene.s4 <- markers.s4.int$gene %>% unique()
gene.ip <- FindVariableFeatures(object = ip.g1.1,assay = 'integrated',nfeatures = 3000)
gene.ip <- gene.ip@assays$integrated@var.features
gene.ip <- intersect(gene.ip,gene.s4)

#int

IP.EN.g1 <- cca_knn(obj1 = stage4,obj2 = ip.g1.1,
                    SCT = F,
                    k = 20,ucomp = 2,
                    assay1 = 'integrated',assay2 = 'integrated',
                    mode = 'int',gene.list = gene.ip)
IP.EN <- cca_knn(obj1 = stage4,obj2 = IP,
                 SCT = F,
                 k = 20,ucomp = 2,
                 assay1 = 'integrated',assay2 = 'integrated',
                 mode = 'int',gene.list = gene.ip)
DimPlot(IP.EN.g1$Seurat,split.by = 'orig.ident',group.by = c('rf','Phase'))
ip.g1.1$knn <- IP.EN.g1$knn$test_cv$result

#rfcv
Idents(ip.g1.1) <- ip.g1.1$knn
marker.ip.g11 <- FindAllMarkers(ip.g1.1,assay = 'integrated',only.pos = T)
marker.ip.g11 <- marker.ip.g11 %>% filter(p_val_adj<.05)
gene.ip.g1<- marker.ip.g11$gene %>% unique()
data.ip <- ip.g1.1@assays$integrated@data[gene.ip.g1,] %>% as.matrix() %>% t
cate.ip <- ip.g1.1$knn %>% as.character() %>% as.factor()
colnames(data.ip) <- TransGeneNames(colnames(data.ip))
rf.ipg1 <- randomForest(data.ip,cate.ip,
                        importance = T,proximity = F)
rfcv.ipg1 <- rfcv(data.ip,cate.ip,ntree=200,
                  #importance = T,proximity = F,
                  cv.fold = 10,step = .75,scale = 'log')



