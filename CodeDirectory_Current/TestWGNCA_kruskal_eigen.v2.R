#source("NiceLM.R")
library(Seurat)
source("WGCNA.Get.Eigen.R")
source("BatchSample.v2.R")
library(coin)


WGCNA_easy<-function(seur,celltype)
{
#print("Subset")
#seur<-SubsetData(seur,WhichCells(seur,celltype))

print("Load")
load("../WGCNA_DE_Results/lists.tidy.current.Robj")
#load("/psych/genetics_data/ssimmons/Perturb/data/WGCNA.genes.clean.Robj")

tab=tab[tab$CellType==celltype,]

print("Run")

ret<-TestWGCNA_module(seur,tab)

return(ret)
}


TestWGCNA_module<-function(seur,tab,ref="GFP",getSeur=F,maxCells=0,seed=1,maxPerc=.3,useCoin=T)
{

print("Pert only")
print(length(seur@ident))
seur=SubsetData(seur,cells=names(seur@ident)[!is.na(seur@meta.data$perturbation)])
print(length(seur@ident))

mods<-unique(tab[,1])

ret=list()

for(i in mods)
{
print(i)

print("Subset")
tab_cur=tab[tab[,1]==i,]

cur<-SubsetData(seur,WhichCells(seur,tab_cur[1,"CellType"]))

print("Make list")

genes<-tab_cur[,"Gene"]

lst<-list()
lst[["Genes"]]=genes

print("Get Module")
#cur<-AddModuleScore(cur,lst)

cur<-getEigen(cur,genes)

if(getSeur)
{
return(cur)

}

print("Perform DE")

cur@meta.data["nUMI"]=scale(cur@meta.data[,"nUMI"])
cur<-SubsetData(cur,names(cur@ident)[!is.na(cur@meta.data$perturbation)])
cur<-SetAllIdent(cur,"perturbation")
if(maxCells>0)
{
print("Downsample!")
lst<-table(cur@meta.data$perturbation)
print(lst)
maxCells=2*median(as.numeric(lst))
print(maxCells)
#cur<-SubsetData(cur,max.cells.per.ident=maxCells)
set.seed(seed)
cur<-batchSample(cur,total=maxCells,maxPerc=maxPerc)
}
cur@meta.data["perturbation"]=factor(cur@meta.data[,"perturbation"])
cur@meta.data["batch"]=factor(cur@meta.data[,"batch"])

#mrk<-lm.pert(cur@meta.data,Cluster1~perturbation+batch+nUMI,Cluster1~batch+nUMI,useBatch=T,dep_var = "Cluster1",ref = ref)
mrk=normal_test(Cluster1~perturbation|batch,data=cur@meta.data,distribution="approximate")
print(mrk)

ret[[i]]=mrk

}

return(ret)

}



drawPlot<-function(out,cutoff=.1)
{
nice<-lapply(out,function(x){return(x$"Pr(>F)"[2])})

dat<-data.frame(as.numeric(nice))

dat["Module"]=names(out)
colnames(dat)[1]="pvalue"

dat["padj"]=p.adjust(dat[,"pvalue"],"fdr")

dat["Sig"]="No"
dat[dat[,"padj"]<cutoff,"Sig"]="Yes"

dat["Sig"]=factor(dat[,"Sig"],c("Yes","No"))


p=ggplot(dat,aes(x=Module,y=-log(pvalue,10),color=Sig))+geom_point()+coord_flip()+scale_color_manual(values=c("red","black"),name=paste("FDR < ",cutoff,sep=""))+geom_vline(xintercept=c(3.5,5.5,7.5,8.5),linetype="dotted")+xlab("")+ylab("-log pvalue, base 10")+ggtitle("Overall Perturbation Effect on WGCNA")

return(p)

}
