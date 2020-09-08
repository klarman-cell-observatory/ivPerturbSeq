library(Seurat)
library(dplyr)
library(tidyr)
print("Load data")
load("updated.with.ASD.Robj")
load("../../SeuratObjects/key.500genes.for.analysis.Robj")

print("Get Expression Info")
#comb=comb[]
#tab<-comb %>% group_by(CellType,Mouse_gene_name) %>% summarise() %>% as.data.frame()
celltypes=unique(key@ident)
out<-lapply(celltypes,function(x){print(x);return(rowMeans(key@data[,key@ident==x]>0))})
out<-lapply(out,data.frame)
names(out)=celltypes
for(x in celltypes){out[[x]]["CellType"]=x}
out<-lapply(out,function(x){x["Gene"]=rownames(x);return(x)})
lapply(out,head)
tab=do.call(rbind,out)
colnames(tab)=c("Percent","CellType","Mouse_gene_name")
#tab["Mouse_gene_name"]=rownames(tab)
print(head(tab))
print(dim(comb))
comb=inner_join(comb,tab)
print(dim(comb))

print("Bin!")
comb["Bin"]=ceiling(comb[,"Percent"]*20)

tab=comb[comb[,"Type"]=="log2FoldChange",] %>% group_by(CellType,Gene,Bin) %>% summarise(logFC_bulk=mean(logFC)) %>% as.data.frame()
comb=comb[comb[,"Type"]!="log2FoldChange",]
tab["sign"]=sign(tab[,4])
print(head(tab))
print(dim(comb))
comb=inner_join(comb,tab)
print(dim(comb))
comb["logFC_sign"]=comb[,"logFC"]*comb[,"sign"]

print("Get Sig genes")
tab<-comb[comb[,"SigLevel"]!="Not",] %>% group_by(Gene,CellType,Bin) %>% summarise() %>% as.data.frame()
colnames(tab)[1]="ComparableGene"

comb<-left_join(tab,comb)

comb<-comb %>% unite(Bin2,CellType,Bin,ComparableGene,remove=F)
tab<-comb %>% group_by(Bin2,Gene,ComparableGene) %>% summarise(Mean=median(logFC_sign)) %>% as.data.frame()
#tab<-comb %>% group_by(Bin2,Gene,ComparableGene) %>% summarise(Mean=mean(logFC_sign>log(1.2))) %>% as.data.frame()
tab2=tab[tab[,2]==tab[,3],]
tab2=tab2[,c("Bin2","Mean")]
colnames(tab2)[2]="Mean_Eff"
tab<-left_join(tab,tab2)

print(head(tab))
res<-tab %>% group_by(ComparableGene,Bin2) %>% summarise(pval=mean(Mean>=Mean_Eff)) %>% as.data.frame()
res["padj"]=p.adjust(res[,"pval"],"fdr")
res=res[order(res[,"pval"]),]
print(res)

write.table(res,"pvalues.ASD.genes.txt",sep="\t",col.names=T,row.names=F,quote=F)





