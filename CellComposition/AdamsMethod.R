library(Seurat)
library(dplyr);
library(tidyr)



##
##Tests to see if differences in composition between perturbations. Fee in:
##seur, the Seurat object
##celltype, the column in meta.data containing celltype info
##perturbation, the column in meta.data containing celltype info
##batch, the column in meta.data containing batch info
##cutoff, the min number of cells to perform test
##
CellComp_Poisson<-function(seur,celltype="CellType",perturatbations="perturbation",batch="batch",cutoff=10)
{

print("Clean Data")
meta=seur@meta.data[,c(celltype,perturatbations,batch)]
colnames(meta)=c("CellType","Pert","Batch")

meta <- meta %>% group_by(CellType,Pert,Batch) %>% summarise(Num=length(CellType)) %>% as.data.frame()

meta2<-meta %>% group_by(Pert,Batch) %>% summarise(Tot=sum(Num)) %>% as.data.frame()

meta=left_join(meta,meta2)

##so meta is dataframe of 5 columns: celltype, perturbation, batch, Number of total cells of that celltype/pert/batch, and number of total cells of that pert/batch

meta=meta[meta[,"Tot"]>cutoff,]


meta["Pert"]=relevel(factor(meta[,"Pert"]),ref="GFP")

lst=list()
for(i in unique(meta[,"CellType"])){lst[[i]]=meta[meta[,"CellType"]==i,]}


print("Fit model!")
out<-lapply(lst,function(cur){
celltype=cur[1,"CellType"]

print(celltype)
cur["logTot"]=log(cur[,"Tot"])
fit<-glm(Num~offset(logTot)+Batch+Pert,data=cur,family="poisson")
tab=summary(fit)
tab=tab$coefficients
tab=data.frame(tab)
tab=tab[grep("Pert",rownames(tab)),]
tab["Gene"]=sub("Pert","",rownames(tab))
tab["CellType"]=celltype
tab=tab[,c(5,6,4,1,2,3)]
colnames(tab)[3]="pval"
return(tab)
})

tab=do.call(rbind,out)

rownames(tab)=NULL

tab=tab[order(tab[,"pval"]),]

tab["padj"]=p.adjust(tab[,"pval"],"fdr")

print("Done!")

return(tab)

}




ErrorBars_Poisson<-function(seur,celltype="CellType",perturatbations="perturbation",batch="batch",cutoff=10)
{

print("Clean Data")
meta=seur@meta.data[,c(celltype,perturatbations,batch)]
colnames(meta)=c("CellType","Pert","Batch")

meta <- meta %>% group_by(CellType,Pert,Batch) %>% summarise(Num=length(CellType)) %>% as.data.frame()

meta2<-meta %>% group_by(Pert,Batch) %>% summarise(Tot=sum(Num)) %>% as.data.frame()

meta=left_join(meta,meta2)

##so meta is dataframe of 5 columns: celltype, perturbation, batch, Number of total cells of that celltype/pert/batch, and number of total cells of that pert/batch

meta=meta[meta[,"Tot"]>cutoff,]


meta["Pert"]=relevel(factor(meta[,"Pert"]),ref="GFP")

lst=list()
for(i in unique(meta[,"CellType"])){for(j in unique(meta[,"Pert"])){lst[[paste(i,j,sep="_")]]=meta[meta[,"CellType"]==i & meta[,"Pert"]==j,]}}
#ieta["Pert"]=relevel(factor(meta[,"Pert"]),ref="GFP")


print("Fit model!")
out<-lapply(lst,function(cur){
celltype=cur[1,"CellType"]

#print(celltype)
cur["logTot"]=log(cur[,"Tot"])
fit<-glm(Num~offset(logTot)+1,data=cur,family="poisson")
tab=summary(fit)
tab=tab$coefficients
tab=data.frame(tab)
tab["Num"]=dim(cur)[1]
#tab=tab[grep("Pert",rownames(tab)),]
#tab["Gene"]=sub("Pert","",rownames(tab))
#tab["CellType"]=celltype
#tab=tab[,c(5,6,4,1,2,3)]
#colnames(tab)[3]="pval"
return(tab)
})

tab=do.call(rbind,out)

tab=data.frame(tab)

tab["Name"]=rownames(tab)
tab["CellType"]=as.character(lapply(rownames(tab),function(x){strsplit(x,"_")[[1]][1]}))
tab["Pert"]=as.character(lapply(rownames(tab),function(x){strsplit(x,"_")[[1]][2]}))

return(tab)


rownames(tab)=NULL

tab=tab[order(tab[,"pval"]),]

tab["padj"]=p.adjust(tab[,"pval"],"fdr")

print("Done!")

return(tab)

}

if(!interactive())
{
print("Load Seurat")
load("../SeurObj_R3.5/key.500genes.Robj")
tab=CellComp_Poisson(key,celltype="CellType",perturatbations="perturbation",batch="batch",cutoff=10)
save(tab,file="CellComposition.Poisson.Robj")
write.table(tab,file="CellComposition.Poisson.txt",sep="\t",quote=F,row.names=F)
p=ggplot(tab,aes(x=Gene,y=CellType,fill=Estimate,size=-log(pval,10)))+geom_point(shape=21)+scale_fill_gradient2(low="blue",mid="white",high="red")+xlab("")+ylab("")+theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("CellComposition.Poisson.pdf",p,width=14,height=5)



#tab<-ErrorBars_Poisson(key)
#save(tab,file="ForError.Bars.Robj")

}
