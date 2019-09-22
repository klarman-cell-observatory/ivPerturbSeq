library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

niceFeaturePlot<-function(seur,gene_name,datainfo=F,low="gray",high="red",size=.5,usePCA=F,useScale=F,dr=NULL)
{
if(is.null(dr))
{
dat=seur@dr$tsne@cell.embeddings

dat=data.frame(dat)
}
if(useScale)
{
seur<-RegressOut(seur,"nGene",genes.regress=c(gene_name,seur@var.genes[1:10]))
}

if(usePCA)
{
dat=seur@pca.rot
}
if(!is.null(dr)){dat=dr;colnames(dat)=c("tSNE_1","tSNE_2")}
print(head(dat))
if(!datainfo)
{print("hi");dat["Gene"]=as.numeric(seur@data[gene_name,]);print("bye")}
else{dat["Gene"]=as.numeric(seur@meta.data[,gene_name])}
if(useScale)
{
print("hi")
print(corner(seur@scale.data))
dat["Gene"]=as.numeric(seur@scale.data[gene_name,])
}
print(head(dat))
p=c()
dat=dat[order(dat$Gene),]
if(usePCA)
{

p=ggplot(dat,aes(x=PC1,y=PC2,color=Gene))
}
else{
p=ggplot(dat,aes(x=tSNE_1,y=tSNE_2,color=Gene))
}
p=p+geom_point(size=size)+scale_colour_gradient(low=low,high=high)+ggtitle(gene_name)
return(p)

}



niceBoxPlot<-function(seur,gene,celltype="CellType",kovswt="ko_vs_wt",vln=F)
{
dat=seur@meta.data[,c(celltype,kovswt)]
dat["Gene"]=as.numeric(seur@data[gene,])
colnames(dat)=c("CellType","Condition","Gene")

p=ggplot(dat,aes(x=CellType,y=Gene,fill=Condition))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))

if(vln)
{
p=ggplot(dat,aes(x=CellType,y=Gene,fill=Condition))+geom_violin(scale="width")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ggtitle(gene)
}


return(p)


}


orderedBoxPlot<-function(seur,gene,cluster,minGene=10)
{

seur<-SubsetData(seur,names(seur@ident)[!is.na(seur@meta.data[,cluster])])

seur<-SetAllIdent(seur,cluster)

tab=table(seur@ident)

seur<-SubsetData(seur,WhichCells(seur,names(tab)[tab>minGene]))


dat=seur@meta.data[,c(cluster,gene)]

colnames(dat)=c("Cluster","Value")

tab <- dat %>% group_by(Cluster) %>% summarise(Med=median(Value)) %>% as.data.frame()

dat["Cluster"]=factor(dat[,"Cluster"],tab$Cluster[order(tab[,"Med"],decreasing=T)])

p=ggplot(dat,aes(x=Cluster,y=Value,fill=Cluster))+geom_boxplot()+theme(legend.position="none")+theme(axis.text.x = element_text(angle = 90, hjust = 1))

return(p)

}




