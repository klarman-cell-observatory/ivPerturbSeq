library(Seurat)



##
##Gets shared variable genes between all batches
##
SharedVariable<-function(seur,x.low.cutoff=1,x.high.cutoff=5,minNum=3,batch="orig.ident",minCells=10)
{
batchs<-unique(seur@meta.data[,batch])

seur<-SetAllIdent(seur,batch)

vars<-c()

for(bat in batchs)
{
print(bat)
temp<-SubsetData(seur,WhichCells(seur,bat))
if(length(temp@ident)>minCells)
{
temp<-FindVariableGenes(temp,x.low.cutoff=x.low.cutoff,x.high.cutoff=x.high.cutoff,do.plot=F)
vars<-c(vars,temp@var.genes)
}
print(" ")
}

print("Combine!")
vars<-table(vars)

genes<-names(vars)[vars>minNum]

print(length(genes))

seur@var.genes<-genes

return(seur)

}


