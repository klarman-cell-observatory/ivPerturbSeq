library(Seurat)
library(Matrix)
library(DESeq2)
library(edgeR)

##modified from seurat code
combData<-function(object,genes=c())
{
genes.use=rownames(object@data)

if(length(genes)>0){genes.use=genes}
data.all=data.frame(row.names = genes.use)
levs=c()
for(i in levels(object@ident)) {
temp.cells=WhichCells(object,i)
levs=c(levs,i)
if (length(temp.cells)==1) data.temp=(object@raw.data[genes.use,temp.cells])
if (length(temp.cells)>1) data.temp=apply(object@raw.data[genes.use,temp.cells],1,sum)
data.all=cbind(data.all,data.temp)
colnames(data.all)[ncol(data.all)]=i

}
colnames(data.all)=levs
return(data.all)
}


######
##Performs DE by combining cells in same batch
######
##Input:
##seur: seurat object, with count data stored in seur@raw.data, condition information (ko_vs_wt) is stored in seur@data.info["ko_vs_wt"]. Assume two conditions, one ko, one for wt
##id: the cluster (or list of clusters) of interest
##kowt: the name of the column in the meta.data table containing ko vs wt information
##combineOn: gives the name of the column in seur@meta.data with batch information (so one factor per batch)
##method: method used for DE. Options: EdgeR, DESeq,
##form: usually leave blank. Can use to include extra covariate infomation. NOT IMPLEMENTED YET.
##minCells: min number of cells per batch for batch to be included.
##minBatches: minimum number of batches per condition for analysis to run (if not met returns NULL).
##minReads: minimum number of reads in genes (requires this number of reads in at least 2 organoids)
##
##Note: ignore genes that do not have at least 10 reads mapped to them in at least two batches.
###########################
##Output:
##res: a data frame with DE information (pvalue, FDR corrected pvalue, logFC, etc). Exact formatting depends on method selected.
##########################################


combDE<-function(seur,id,kowt="ko_vs_wt",combineOn="orig.ident",method="DESeq",form="",minCells=20,minBatches=2,minReads=10,mapIt=c())
{
print("subsample")
seur<-SubsetData(seur,WhichCells(seur,id))

seur<-SetAllIdent(seur,combineOn)
print("combine data")
dat<-combData(seur)

print("Filter 10X Lanes")
keep=names(summary(seur@ident))[summary(seur@ident)>minCells]
numOrg=length(keep)
dat=dat[keep]
print(dim(dat))
print("get ko vs wt data")
print(head(dat))
cols=colnames(dat)
print(cols)
lst=c()
lst=cols


keep<-strsplit(form,"+",fixed=T)[[1]]

print(c(kowt,combineOn,keep))
print(head(seur@meta.data))

keep<-strsplit(form,"+",fixed=T)[[1]]
val=seur@meta.data[,c(kowt,combineOn,keep)]
val=val[!duplicated(val[,2]),]
rownames(val)=val[,2]
lst_backup=lst
lst=as.character(val[lst,1])



if(length(mapIt)>0){
lst2=c()
lst=mapIt[as.character(lst)];print(class(lst));lst=as.character(lst)}


print(lst)
print(summary(lst))
if(min(summary(factor(lst)))<minBatches)
{
print("Not enough batches!")
return(NULL)
}



if(length(summary(factor(lst)))<2)
{
print("Not enough batches!")
return(NULL)
}


print("Make col data")
colDat=cbind(factor(lst))
colDat=data.frame(colDat)
colnames(colDat)="condition"
rownames(colDat)=colnames(dat)
colDat[1]=factor(colDat[,1])
print(class(colDat[,1]))
if(length(keep)>0)
{
colDat[keep]=val[lst_backup,keep]
}


print(head(colDat))

print("Filter")

if(method=="DESeq")
{
design= ~ condition
if(nchar(form)>0)
{
design=as.formula(paste("~",form,"+ condition",sep=""))
}
print(design)

dds <- DESeqDataSetFromMatrix(countData = dat,colData = colDat,design = design)
dds <- dds[ rowSums(counts(dds) > minReads)>=2, ]
print("run!")
dds <- DESeq(dds)
out=data.frame(results(dds))
out=out[order(out$pvalue),]
return(out)
}

if(method=="EdgeR")
{
print("Filter")
group=colDat[,1]
dat=dat[rowSums(dat > minReads)>2,]
print("run")
y <- DGEList(counts=dat,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
 y <- estimateDisp(y,design)
 fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  res=data.frame(qlf$table)
  res["padj"]=p.adjust(res$PValue,"fdr")
res=res[order(res$PValue),]
res["Gene"]=rownames(res)

numCol=dim(res)[2]
print(head(res))
res=res[c(numCol,1:(numCol-1))]
print(head(res))
return(res)
}


return(NULL)

}








