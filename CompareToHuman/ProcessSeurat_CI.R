library(Seurat)
library(openxlsx)
library(Matrix)
library(boot)



##This code goes through each module and calculates the average correlation between genes in the associated cell types.
##Uses permutations to get p values and bootstrap to get CI
##also gets average/SD of permutations for normalization
averageCor<-function(cur,indices=c(),genes=c(),tot=1)
{

if(length(indices)>0)
{
cur=data.frame(as.matrix(cur))[indices,]
}



cur=t(cur)

if(length(genes)>0)
{
genes=intersect(genes,rownames(cur))
print(dim(cur))
cur=t(as.matrix(t(cur)[,genes]))

print(dim(cur))
}

dat=as.matrix(t(cur))

#tot=sum(exp(data.frame(as.matrix(cur))[,1])-1)
dat=log(1000000*(exp(dat)-1)/tot+1)

aveCOR=mean(cor(dat))

return(aveCOR)

}

testCoexpression<-function(cur,genes,numPerm=1000)
{

print("Set up for permuting")
tab=data.frame(Matrix::rowMeans(cur@data>0))
print(head(tab))
tab["Temp"]=1

rownames(tab)=rownames(cur@data)

tab=tab[tab[,1]>.05,]

tab=tab[order(tab[,1]),]


tab["Pos"]=1:dim(tab)[1]

tab["Bin"]=floor(100*tab[,"Pos"]/(dim(tab)[1]+1))

genes=intersect(genes,rownames(tab))

bins=tab[genes,"Bin"]



perms=c()

seur=cur

tot=sum(exp(cur@data[,1])-1)

cur=t(cur@data[genes,])

print("Run Boot Strap!")
modCor=boot(cur,averageCor,1000,tot=tot)
print("Get CI!")
modCor=boot.ci(modCor,type="perc")

print(str(modCor))

print(modCor$perc)

vals=c(modCor$t0,modCor$perc[4],modCor$perc[5])

modCor=averageCor(cur,tot=tot)

#return(vals)
perms=c()


print("Permute")
for(i in 1:numPerm)
{
genes_new=c()


while(length(genes_new)<length(genes))
{
genes_new=c()
for(b in bins)
{
genes_new<-c(genes_new,sample(rownames(tab)[tab[,"Bin"]==b],1))
}
genes_new=unique(genes_new)
}


cur=t(seur@data[genes_new,])
perms=c(perms,averageCor(cur,tot=tot))

}

return(c(modCor,((sum(perms>=modCor)+1)/(numPerm+1)),mean(perms) ,sd(perms),vals[1],vals[2],vals[3]))

}

##
##Gets results of different tests in modules
##
TestModules<-function(seurfile,moduleFile)
{
print("Load data")
seur=readRDS(seurfile)

if("CellType_new" %in% colnames(seur@meta.data)){seur@meta.data["CellType"]=seur@meta.data[,"CellType_new"]}

seur<-SetAllIdent(seur,"CellType")

seur<-SubsetData(seur,names(seur@ident)[!is.na(seur@ident)])

print("Load module")
load(moduleFile)

mods=unique(tab[,1])

minCells=100

res=list()

for(mod in mods)
{
print(paste("Current Module: ",mod,sep=""))
celltype=tab[tab[,1]==mod,"CellType"]
celltype=celltype[1]

genes=tab[tab[,1]==mod,"Gene"]

genes=intersect(genes,rownames(seur@data))

#if(length(genes)<3)
#{
#next
#}

if(sum(seur@ident==celltype)<minCells)
{
next
}

print("Subset data")
cur<-SubsetData(seur,names(seur@ident)[seur@ident==celltype])

print("Analysis!")
print("First, see if the genes are expressed!")
numExp<-Matrix::rowMeans(cur@data[genes,]>0)

numGenes=length(genes)
genes_temp=genes[numExp>.01]
numGenesExp=length(genes_temp)
genes=genes[numExp>.05]
numGenesExpHigh=length(genes)

if(numGenesExpHigh<3)
{
next
}

print("Next, coexpression information")
test<-testCoexpression(cur,genes)

results<-data.frame(c("Total Genes","Total Genes >1%","Total Genes >5%","Correlation","pval","Mean","SD","CorBoot","CI1","CI2"))
colnames(results)[1]="Name"
results["Values"]=c(numGenes,numGenesExp,numGenesExpHigh,test[1],test[2],test[3],test[4],test[5],test[6],test[7])
results["Module"]=mod
results["CellType"]=celltype

print(results)

res[[mod]]=results

}


return(res)

}



if(!interactive())
{
args=commandArgs(trailingOnly=TRUE)
seurfile=args[1]
outfile=args[2]
if(length(args)==2)
{
res<-TestModules(seurfile)
}
else
{
modfile=args[3]
res<-TestModules(seurfile,modfile)
}


write.xlsx(res,outfile)
}

