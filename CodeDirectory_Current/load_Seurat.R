
library(Seurat)
library(stringr)
library(Matrix)

#loads all 10X lanes from a given directort
dir10X<-function(dir="",dat=NULL,lst=c(),makeSeurat=T,minGenes=500,regress=c("nUMI"),species="human",num=c())
{
if(length(lst)==0)
{
print(paste("ls ",dir,"/*/*/outs/filt*/* | grep : | sed 's/://g'",sep=""))
lst=system(paste("ls ",dir,"/*/*/outs/filt*/* | grep : | sed 's/://g'",sep=""),intern=T)

}
if(length(num)>0)
{
lst=lst[num]
}
print(lst)
if(is.null(dat))
{
print("Read in!")
dat=Read10X(lst)
print("Fix colnames!")
cols=colnames(dat)
cols_new=c()
for(col in cols)
{
start=str_sub(col,1,1)
cur=col
if(start %in% c("A","T","G","C")){cur=paste("1_",cur,sep="")}
cols_new<-c(cols_new,cur)
}

colnames(dat)=cols_new
#print("Return!")
}
print(paste("Dims: ",toString(dim(dat))))
if(!makeSeurat){return(dat)}

print("Make object!")
seur<-CreateSeuratObject(dat,"Seurat",min.genes=minGenes,normalization.method="LogNormalize",scale.factor=1000000)

print("Get variable genes!")
seur<-FindVariableGenes(seur,x.low.cutoff=1,do.plot =F)

print("Regress out!")
if(length(regress)>0)
{
seur<-ScaleData(seur,genes.use=seur@var.genes,vars.to.regress=regress)
}
else{
seur<-ScaleData(seur,genes.use=seur@var.genes)
}
print("Run PCA!")
seur<-RunPCA(seur,pcs.compute=60)
#if(species=="human")
#{
#print("Get ribosomal and mito")
#bot=colSums(seur@raw.data)
#mito=grep("^MT-",rownames(seur@raw.data))
#if(length(mito)>0)
#{
#top=colSums(seur@raw.data[mito,])
#mn=top/bot
#seur@data.info["mito"]=names(seur@ident)
#}

#mito=grep("^RP[S,L]",rownames(seur@raw.data))
#if(length(mito)>0)
#{
#top=colSums(seur@raw.data[mito,])
#mn=top/bot
#seur@data.info["ribo"]=names(seur@ident)
#}


#}

print("load!")


print("Return!")

return(seur)

}


library(dplyr)
library(tidyr)
##Returns the percent cells in each cluster expressing given genes
AveExpress<-function(seur,genes.use=c(),percent=T)
{
genes.use=intersect(genes.use,rownames(seur@data))
if(length(genes.use)<1){return(c())}

dat=data.frame(as.matrix(t(seur@data[genes.use,])))

dat["clust"]=seur@ident

if(percent)
{
dat<- dat %>% gather(Gene,Expression,-clust)%>% group_by(Gene,clust) %>% summarise(Percent=mean(Expression>0)) %>% spread(Gene,Percent) %>% as.data.frame()
}
else
{
dat<- dat %>% gather(Gene,Expression,-clust)%>% group_by(Gene,clust) %>% summarise(Average=log(mean(exp(Expression+1)))) %>% spread(Gene,Average) %>% as.data.frame()
}
rownames(dat)=dat[,"clust"]

dat=dat[,2:dim(dat)[2]]

return(dat)
}


