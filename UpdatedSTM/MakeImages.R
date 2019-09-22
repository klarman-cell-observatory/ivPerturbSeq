library(Seurat)
source("/psych/genetics_data/ssimmons/Perturb/BasicPipeline/RunSTM_clean.R")
source("/psych/genetics_data/ssimmons/SingleCell2/niceDraw.R")
library(cowplot)
library(grid)
library(gridExtra)
library(stm)
library(dplyr)
library(tidyr)


##Takes in the nice Seurat object and the output of RunSTM
GetTSNE_all<-function(seur,out,outdir)
{
savename=paste(outdir,"/FeaturePlot_Topic_all.pdf",sep="")

res=list()

for(Topic in 1:5)
{
print("Clean Up Seurat")
inter=intersect(names(seur@ident),names(out[[3]]@ident))

seur<-SubsetData(seur,inter)



seur@meta.data[inter,paste("Topic",toString(Topic),sep="")]=out[[3]]@meta.data[inter,paste("X",toString(Topic),sep="")]

print("Plot")

p=niceFeaturePlot(seur,paste("Topic",toString(Topic),sep=""),datainfo=T,low="grey",size=1)+ggtitle(paste("Topic",toString(Topic)))+theme(legend.position="none")

res[[Topic]]=p

}


p=plot_grid(res[[1]],res[[2]],res[[3]],res[[4]],res[[5]],ncol=1,nrow=5)

ggsave(savename,p,height=30,width=6)

print("Done")

return(p)

}


##Takes in the nice Seurat object and the output of RunSTM
GetTSNE<-function(seur,out,Topic,outdir)
{
savename=paste(outdir,"/FeaturePlot_Topic_",toString(Topic),".pdf",sep="")

print("Clean Up Seurat")
inter=intersect(names(seur@ident),names(out[[3]]@ident))

seur<-SubsetData(seur,inter)



seur@meta.data[inter,paste("Topic",toString(Topic),sep="")]=out[[3]]@meta.data[inter,paste("X",toString(Topic),sep="")]

print("Plot")

p=niceFeaturePlot(seur,paste("Topic",toString(Topic),sep=""),datainfo=T,low="grey",size=2)+ggtitle(paste("Topic",toString(Topic)))

ggsave(savename,p,height=7,width=7)

print("Done")

}


##Does the dotplot, only requires results of testSTM
GetDotPlot<-function(ret,outdir)
{
savename=paste(outdir,"/DotPlot.pdf",sep="")


p=ggplot(ret,aes(x=Topic,y=perturbation,fill=Effect_Size,size= -log(pval_pert,10)))+geom_point(shape=21)+scale_fill_gradient2(low="blue",high="red",mid="white")

ggsave(savename,p,width=7,height=10)

return(p)

}


##Prints the logPvalue
GetPvaluePlot<-function(out,outdir)
{
savename=paste(outdir,"/PvaluePlot.pdf",sep="")

 res<-testSTM_all(out[[1]],out[[2]],logistic=F,form1 = Topic ~ batch + perturbation + nUMI,form2 = Topic ~ batch + nUMI,test="normal")

res["Significant"]=res[,"padj"]<.05

res["Topic"]=sub("Topic","",res[,"Topic"])

p=ggplot(res,aes(x=Topic,y=-log(pvalue,10),color=Significant))+geom_point(size=4)+ylab("-log10 p-value")+scale_color_manual(values=rev(c("red","grey")))+theme(legend.position="none")


ggsave(savename,p,height=7,width=7)

return(p)

}



##Gets TopGenes
GetTopGenes_all_v2<-function(out,outdir,max.words=10)
{

savename=paste(outdir,"/TopWords.Topic.all.pdf",sep="")

res<-list()

print("Get Scores")

#lab=labelTopics(out[[1]],n=max.words)$score[Topic,]

wordcounts=out[[1]]$settings$dim$wcounts$x

logbeta=out[[1]]$beta$logbeta[[1]]


emp.prob <- log(wordcounts) - log(sum(wordcounts))

lift <- logbeta - rep(emp.prob, each = nrow(logbeta))

ldascore <- exp(logbeta) * (logbeta - rep(colMeans(logbeta), each = nrow(logbeta)))

lift=ldascore

#lift=log(lift+1)

lift=t(lift)

lift=data.frame(lift)

colnames(lift)=1:dim(lift)[2] 

rownames(lift)=out[[1]]$vocab

lift["Genes"]=out[[1]]$vocab

tab<-lift %>% gather(Topic,Score,-Genes)

tab=tab[order(tab[,"Score"],decreasing=T),]

tab=tab[!duplicated(tab[,"Genes"]),]


tab["Use"]=0

print(head(tab))

for(i in 1:5){
num=sum(tab[,"Topic"]==i)
val=num-max.words
print(val)
tab[tab[,"Topic"]==i,"Use"]=c(rep(1,max.words),rep(0,val))

}


print(head(tab))

tab=tab[tab[,"Use"]==1,]


print(head(tab))

#tab=data.frame(lift[lab,Topic])

#tab["Genes"]=lab


#colnames(tab)[1]="Score"


tab["Topic"]=sub("^","Topic ",tab[,"Topic"])

tab[,"Genes"]=factor(tab[,"Genes"],tab[order(tab[,"Score"]),"Genes"])


p=ggplot(tab,aes(x=Genes,y=Score))+geom_bar(stat="identity")+ylab("Score")+coord_flip()+theme(legend.position="none")+facet_wrap(~Topic,ncol=1,scale="free")+theme(strip.background=element_blank())





ggsave(savename,p,height=30,width=6)

return(p)

}



##Gets TopGenes
GetTopGenes_all<-function(out,outdir,max.words=10)
{

savename=paste(outdir,"/TopWords.Topic.all.pdf",sep="")

res<-list()

for(Topic in 1:5)
{
print("Get Scores")

lab=labelTopics(out[[1]],n=max.words)$score[Topic,]

wordcounts=out[[1]]$settings$dim$wcounts$x

logbeta=out[[1]]$beta$logbeta[[1]]


emp.prob <- log(wordcounts) - log(sum(wordcounts))

lift <- logbeta - rep(emp.prob, each = nrow(logbeta))

ldascore <- exp(logbeta) * (logbeta - rep(colMeans(logbeta), each = nrow(logbeta)))

lift=ldascore

#lift=log(lift+1)

lift=t(lift)

lift=data.frame(lift)

colnames(lift)=1:dim(lift)[2] 

rownames(lift)=out[[1]]$vocab

tab=data.frame(lift[lab,Topic])

tab["Genes"]=lab


colnames(tab)[1]="Score"


tab[,2]=factor(tab[,2],tab[order(tab[,1]),2])

p=ggplot(tab,aes(x=Genes,y=Score))+geom_bar(stat="identity")+ylab("Score")+coord_flip()+ggtitle(paste("Topic",Topic))+theme(legend.position="none")

res[[Topic]]=p
}



p=plot_grid(res[[1]],res[[2]],res[[3]],res[[4]],res[[5]],ncol=1,nrow=5)

ggsave(savename,p,height=30,width=6)

return(p)

}


##Gets TopGenes
GetTopGenes<-function(out,Topic,outdir,max.words=10)
{

savename=paste(outdir,"/TopWords.Topic.",toString(Topic),".pdf",sep="")

print("Get Scores")

lab=labelTopics(out[[1]],n=max.words)$score[Topic,]

wordcounts=out[[1]]$settings$dim$wcounts$x

logbeta=out[[1]]$beta$logbeta[[1]]


emp.prob <- log(wordcounts) - log(sum(wordcounts))

lift <- logbeta - rep(emp.prob, each = nrow(logbeta))

ldascore <- exp(logbeta) * (logbeta - rep(colMeans(logbeta), each = nrow(logbeta)))

lift=ldascore

#lift=log(lift+1)

lift=t(lift)

lift=data.frame(lift)

colnames(lift)=1:dim(lift)[2] 

rownames(lift)=out[[1]]$vocab

tab=data.frame(lift[lab,Topic])

tab["Genes"]=lab


colnames(tab)[1]="Score"


tab[,2]=factor(tab[,2],tab[order(tab[,1]),2])

p=ggplot(tab,aes(x=Genes,y=Score))+geom_bar(stat="identity")+ylab("Score")+coord_flip()+ggtitle(paste("Topic",Topic))

ggsave(savename,p)

}



##Compares to WGNCA results
CompareWGCNA<-function(out,celltype,outdir,test=F)
{
savename=paste(outdir,"/WGCNA.and.STM.pdf",sep="")


load("/psych/genetics_data/ssimmons/Perturb/Updated_Data_March_11_2019/WGCNA_DE_Results/lists.tidy.current.Robj")

print("Extract Topics")
topics=out[[3]]@meta.data[,grep("^X",colnames(out[[3]]@meta.data))]

colnames(topics)=sub("^X","Topic",colnames(topics))

print("getWGCNA")

tab=tab[tab[,4]==celltype,]

uniq=unique(tab[,1])

seur=out[[3]]

lst<-list()

for(i in 1:length(uniq))
{
lst[[uniq[i]]]=tab[tab[,1]==uniq[i],2]
}

print("WGCNA Score")
seur<-AddModuleScore(seur,lst,enrich.name="WGCNA")

wgcna=seur@meta.data[,grep("^WGCNA",colnames(seur@meta.data))]

print(head(topics))

print(head(wgcna))
print(head(seur@meta.data))


print(uniq)

comb<-cbind(topics,wgcna)

print("Get Correlation")
COR=data.frame(cor(topics,wgcna,method="spearman"))

COR["Topic"]=rownames(COR)

print("Plot")
tab<- COR %>% gather(WGCNA,Score,-Topic)


print(head(tab))

colnames(tab)[3]="Correlation"

#names(uniq)=names(lst)

lst<-names(lst)
names(lst)=colnames(seur@meta.data)[grep("^WGCNA",colnames(seur@meta.data))]

tab["WGCNA"]=lst[tab[,"WGCNA"]]

if(length(uniq)==1){tab[,"WGCNA"]=uniq[1]}

if(test){ret<-list();ret[[1]]=wgcna;ret[[2]]=topics;ret[[3]]=COR;ret[[4]]=tab;print(length(ret));return(ret)}

#p<-ggplot(tab,aes(x=WGCNA,y=Correlation,fill=Topic))+geom_tile()+scale_fill_gradient2(low="blue",mid="white",high="red")+xlab("")+ylab("")+theme(axis.ticks=element_blank())

p<-ggplot(tab,aes(x=WGCNA,fill=Correlation,y=Topic))+geom_tile()+scale_fill_gradient2(low="blue",mid="white",high="red")+xlab("")+ylab("")+theme(axis.ticks=element_blank(),axis.line=element_blank())+theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(savename,p)

return(p)

}


##Runs for all!
##Takes in:
##seur, the clean Seurat object for this celltyp
##out, the output of RunSTM for this cell type
##res, the output of testSTM for this cell type, run on runSTM
##outdir, the outdirectory to pring images to for this cell type.
##numTopics, the number of topics used
RunAll<-function(seur,out,res,outdir,numTopics=5,name="")
{

print("WGCNA")

p0<-CompareWGCNA(out,name,outdir)

print("First is the dotplot!")



p3=GetDotPlot(res,outdir)

print("Next are the TSNE FeaturePlots and Top words!")
for(Topic in 1:numTopics)
{
print(Topic)
GetTSNE(seur,out,Topic,outdir)
GetTopGenes(out,Topic,outdir,max.words=10)
}

p2=GetTSNE_all(seur,out,outdir)
#p1=GetTopGenes_all(out,outdir,max.words=10)
p1=GetTopGenes_all_v2(out,outdir,max.words=10)


print("Next are the p-value plots!")
p4=GetPvaluePlot(out,outdir)


tot=list()
tot[["TopGenes"]]=p1
tot[["TSNE"]]=p2
tot[["DotPlot"]]=p3
tot[["PvaluePlot"]]=p4
tot[["WGCNA"]]=p0


p0=plot_grid(tot[[4]],tot[[5]],ncol=2,labels=c('','E.'))
p2=plot_grid(tot[[3]],p0,ncol=1,rel_heights=c(3,1),labels = c('C.', 'D.'))
p1=plot_grid(tot[[1]],tot[[2]],ncol=2,labels = c('A.', 'B.'))

p=plot_grid(p1,p2,ncol=2,rel_widths=c(1,1.5))

tot[["Topic"]]=p1
tot[["Tests"]]=p2
tot[["All"]]=p

savename=paste(outdir,"/Combined",name,".pdf",sep="")

ggsave(savename,p,width=14,height=14)


savename=paste(outdir,"/list.images.Robj",sep="")


#save(tot,file=savename)

}





if(!ineractive())
{

print("Load Seurat")
#load("/psych/genetics_data/ssimmons/Perturb/RunSTM/Final/STM_5Topics/Seurat.objects.Robj")


load("/psych/genetics_data/ssimmons/Perturb/Updated_Data_March_11_2019/SeuratObjects/CellType/all.celltypes.Robj")

#load("/psych/genetics_data/ssimmons/Perturb/Updated_Data_March_11_2019/Results5Topics/CellType/all.celltypes.Robj")

#dir="/psych/genetics_data/ssimmons/Perturb/RunSTM/Final/STM_5Topics"

dir="/psych/genetics_data/ssimmons/Perturb/Updated_Data_March_11_2019/UpdatedSTM/Results5Topics"

for(i in names(seurs))
{
print(i)
print("load data")
seur=seurs[[i]]
load(paste(dir,"/stm.",i,".Robj",sep=""))
load(paste(dir,"/test.",i,".Robj",sep=""))
outdir=paste(dir,"/Results_new/",i,sep="")

system(paste("mkdir",outdir))

print("Run!")
RunAll(seur,out,res,outdir,name=i)

}



}

