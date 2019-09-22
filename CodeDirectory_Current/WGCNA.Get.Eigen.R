library(WGCNA)
library(Seurat)

#seur is Seurat object
#genes, the array describing the current module
getEigen<-function(seur,genes,mod="No Name",nam="Cluster1")
{

print("Prep Data")
mn=Matrix::rowSums(seur@data)

dat=data.frame(t(as.matrix(seur@data)))

dat=dat[,mn>0]

cols<-rep("grey",dim(dat)[2])

results=list()

#modules<-as.character(unique(tab[,1]))

#for(mod in modules)
#{
print("Running on:")
print(mod)
cur<-cols

#genes<-tab[tab[,1]==mod,2]
genes<-intersect(genes,colnames(dat))
names(cur)=colnames(dat)
cur[genes]="red"
eigen<-moduleEigengenes(dat,cur,excludeGrey = T,verbose=T)

#results[[mod]]=eigen
#}


#res<-lapply(results,function(x){x$eigengenes})
#tab<-data.frame(res)
#colnames(tab)=names(res)

seur@meta.data[nam]=eigen$eigengenes

return(seur)

}
