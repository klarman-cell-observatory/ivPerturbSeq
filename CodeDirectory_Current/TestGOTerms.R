library(dplyr);library(tidyr)
#library(Seurat)
library(Seurat)


##Tests for perturbation effects on different GO Terms
ScoreGOTerms<-function(seur,celltype,minPerc=.05,minNum=5,GOTerms="/stanley/levin_dr/ssimmons/SingleCell2/GOEnrich/GO.Clean.noFactors.Robj")
{
print("Load GO Terms")
load(GOTerms)
print("Subsample")
seur<-SubsetData(seur,WhichCells(seur,celltype))

genes=rownames(seur@data)
mn=Matrix::rowMeans(seur@data>0)
genes=genes[mn>minPerc]
seur@meta.data["nGene"]=scale(seur@meta.data[,"nGene"])

seur@meta.data["perturbation"]=factor(seur@meta.data[,"perturbation"])

seur@meta.data["perturbation"]=relevel(seur@meta.data[,"perturbation"],ref="GFP")

print("Apply Test!")
iter=0
##Gets average score(in log TPM) for each term
scores<-lapply(Terms,function(x){
iter<<-iter+1
if(iter %% 100 ==0){
print(iter)}
x=intersect(x,genes)
if(length(x)<minNum){return(NULL)}

scor=Matrix::colMeans(seur@data[x,])

meta=seur@meta.data

meta["Score"]=scor


fit<-lm(Score~nGene+perturbation+batch,data=meta)

coef=summary(fit)$coefficients

coef=data.frame(coef)

coef=coef[grep("pert",rownames(coef)),]

coef["Pert"]=rownames(coef)

return(coef)

})

print("Clean up!")

names(scores)=names(Terms)

scores=scores[-which(sapply(scores, is.null))]

for(i in names(scores)){scores[[i]]["Name"]=i}

scores=do.call(rbind,scores)

scores["CellType"]=celltype

colnames(scores)[4]="pval"
scores["padj"]=p.adjust(scores[,4])

scores=scores[order(scores[,"pval"]),]

return(scores)

}



##Test for all celltypes
ScoreGOTerms_All<-function(seur,minPerc=.05,minNum=5)
{

celltypes=as.character(unique(seur@ident))

out<-lapply(celltypes,function(celltype){print(celltype);return(ScoreGOTerms(key,celltype,minPerc=minPerc,minNum=minNum))})

return(out)

}



