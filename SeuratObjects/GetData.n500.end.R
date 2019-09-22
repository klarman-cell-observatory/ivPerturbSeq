set.seed(1)
savedir="."
codedir="../CodeDirectory_Current"

source(paste(codedir,"/load_Seurat.R",sep=""),chdir=T)
source(paste(codedir,"/SharedVariable.R",sep=""),chdir=T)
source(paste(codedir,"/niceDraw.R",sep=""),chdir=T)


load(paste(savedir,"/key.500genes.Robj",sep=""))

##Iterate over all cell types to get initial subclustering, in a list of Seurat objects
seurs<-list()

 for(i in unique(key@ident)){print(i);cur<-SubsetData(key,WhichCells(key,i));cur<-SharedVariable(cur);cur<-ScaleData(cur,genes.use=cur@var.genes,vars.to.regress=c("nUMI","batch"));cur<-RunPCA(cur);cur<-RunTSNE(cur,dims.use=1:15,verbose=T);cur<-FindClusters(cur,dims.use=1:15);seurs[[i]]=cur;print(i)}


save(seurs,file="Clustered.Celltypes.Robj")

##Get Subset For Actual Analysis
PerturbCode=codedir
source(paste(PerturbCode,"/BatchSample.v2.R",sep=""),chdir=T)



key<-SubsetData(key,names(key@ident)[!is.na(key@meta.data$perturbation)])

key<-SubsetData(key,names(key@ident)[!(key@meta.data$perturbation %in% c("Suv420h1","Map1a","Ankrd11"))])

keep<-c()

for(celltype in unique(key@ident))
{
cells<-SubsetData(key,WhichCells(key,celltype))

cells<-SetAllIdent(cells,"perturbation")

lst<-table(as.character(cells@ident))


cells<-SubsetData(cells,names(cells@ident)[cells@ident %in% names(lst)[lst>10]])


lst<-table(as.character(cells@ident))

med=2*median(lst)

set.seed(123)
cells<-batchSample(cells,total=med,maxPerc=1.0)




lst<-table(as.character(cells@ident))


cells<-SubsetData(cells,names(cells@ident)[cells@ident %in% names(lst)[lst>10]])

keep<-c(keep,names(cells@ident))

}


key<-SubsetData(key,keep)






save(key,file=paste(savedir,"/key.500genes.for.analysis.Robj",sep=""))
