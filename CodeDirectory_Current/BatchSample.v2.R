library(Seurat)


##samples so that evenly distributed between batches
batchSample<-function(seur,batchs="batch",total=100,perturbation="perturbation",maxPerc=.3)
{
seur<-SubsetData(seur,names(seur@ident)[!is.na(seur@meta.data[,perturbation])])

seur<-SetAllIdent(seur,perturbation)

seur2<-seur



ret<-c()

for(pert in unique(seur2@ident))
{

seur<-SubsetData(seur,WhichCells(seur,pert))

uniq=unique(seur@meta.data[,batchs])
seur@meta.data["samp"]=0

for(batch in uniq)
{
names=names(seur@ident)[seur@meta.data[,batchs]==batch]

seur@meta.data[names,"samp"]=sample(1:length(names))

}

vals<-seur@meta.data[order(seur@meta.data[,"samp"]),"samp"][1:total]

#print(tail(vals))
vals<-vals/(1:total)

#print(tail(vals))

#print(sum(vals<maxPerc))

cutoff=max((1:total)[vals<maxPerc & !is.na(vals)])

#if(cutoff<1){cutoff=total}

#print(cutoff)

if(cutoff>1)
{
names=names(seur@ident)[order(seur@meta.data[,"samp"])[1:cutoff]]

ret<-c(ret,names)
}

seur<-seur2


}

return(SubsetData(seur,ret))

}


