#library(asbio)
library(stm)
library(Seurat)
library(dplyr)
library(tidyr)
library(Matrix)



source("../CodeDirectory_Current/PertLM.R",chdir=T)

##
##Runs STM for a particular CellType
##
runSTM<-function(seur,celltypes,numberPCs=0,forForm=c("perturbation","batch"),numTopics=10,genes=c(),minPert=30,getMeta=F)
{
print("Subset!")
if(length(celltypes)>0)
{
seur=SubsetData(seur,WhichCells(seur,celltypes))



seur<-FindVariableGenes(seur,x.low.cutoff=1,do.plot=F)
}
#seur<-ScaleData(seur,genes.use=seur@var.genes,vars.to.regress=c("nUMI"))
print("Get Common Perts!")
lst=table(seur@meta.data$perturbation)
#if(length(lst)>28){minPert=max(minPert,as.numeric(lst[order(as.numeric(lst),decreasing=T)[28]]-1))}


seur<-SubsetData(seur,names(seur@ident)[seur@meta.data$perturbation %in% names(lst)[lst>minPert]])
print("Number Perts before:")
print(length(lst))



print("Number Perts after:")
print(sum(lst>minPert))

#print("Clean up Perts with big overrepresentation")
#lst=table(seur@meta.data$perturbation)
#med=median(lst)
#max_size=2*med
#seur=SubsetData(seur,sample(names(seur@ident)))
#seur<-SetAllIdent(seur,"perturbation")
#seur=batchSample(seur,total=max_size)
#seur=SubsetData(seur,names)
#seur=SubsetData(seur,max.cells.per.ident=max_size)
#print("Number Cells:")
#print(length(seur@ident))



mn=rowMeans(as.matrix(seur@data)>0)

dat=seur@raw.data[rownames(seur@data)[mn>.05 & mn<.9],names(seur@ident)]



dat=dat[grep("^mt-",rownames(dat),invert=T),]
dat=dat[grep("^Rp[s,l]",rownames(dat),invert=T),]




#if(length(genes)>0)
#{
#print("Use supplied gene list!")
#genes=intersect(genes,rownames(seur@data))
#dat=seur@data[genes,]
#}


print("Remove Batch Effect Genes")

print(dim(dat))

tab=data.frame(t(as.matrix(dat)))
colnames(tab)=rownames(dat)
tab["batch"]=seur@meta.data[rownames(tab),"batch"]
print("get info")
tab <- tab %>% gather(Gene,Expression,-batch) %>% group_by(Gene,batch) %>% summarise(Mean=mean(Expression>0)) %>% group_by(Gene) %>% summarise(x=max(Mean),y=min(Mean),diff=x-y) %>% as.data.frame() 

print("Filter")

print(head(tab))


genes=tab[tab[,"y"]>0,"Gene"]


genes=intersect(genes,rownames(dat))

dat=dat[genes,]

dat=t(dat)

print(dim(dat))

print("Make formula!")

if(numberPCs>0)
{
for(i in 1:numberPCs)
{
forForm<-c(forForm,paste("s(PC",as.character(i),")",sep=""))
}
}



form=paste(forForm,collapse="+")
form=paste("~",form,sep="")
form=as.formula(form)

meta=seur@meta.data
pc=seur@dr$pca@cell.embeddings
#meta[colnames(pc)]=pc

print(form)

print("Run STM!")
print(dim(dat))
print(dim(meta))
#dat=data.frame(as.matrix(dat)
#dat=log(dat+1)
res<-stm(dat,K=numTopics,prevalence=form,data=data.frame(meta),LDAbeta=T,interactions=F)
if(getMeta)
{
lst=list()
lst[["res"]]=res
lst[["meta"]]=meta
seur@meta.data[colnames(data.frame(res$theta))]=data.frame(res$theta)
lst[["seur"]]=seur
return(lst)
}
return(res)

}


##
##Tests for differences in each coefficient of the STM output
##
testSTM<-function(res,meta,form="~batch+perturbation",reference="GFP",variable="perturbation",as.df=T,rawtheta=T,logistic=T,numRep=1000,topics=c(),getResid=F)
{
print("Get theta")
theta=res$theta

if(!rawtheta)
{
theta=thetaPosterior(res,1)

theta <- do.call(rbind, theta)
}
if(length(topics)>0){theta=theta[,topics];theta=data.frame(t(theta));for(i in colnames(theta)){theta[i]=theta[,i]/sum(theta[,i])};theta=t(theta)}

if(logistic)
{
theta=log(theta/(1-theta))
}
print("Combine with metadata")
#theta=data.frame(theta)


numTopics=dim(theta)[2]

meta=cbind(theta,meta)

meta=data.frame(meta)


print("Fix Perturbation")

print(head(meta))
meta["perturbation"]=factor(meta[,"perturbation"])

if(nchar(reference)>0)
{
meta <- within(meta, perturbation <- relevel(perturbation, ref = reference))
}

results=list()

form=as.character(form)
#print(form)
form=paste("Topic",form,sep="")
form=as.formula(form)
print(form)

print("Perform linear regression!")

meta_backup=meta
models=list()
for(i in 1:numTopics)
{
print(paste("Topic:",toString(i)))
colnames(meta)[i]="Topic"
fit=lm(form,data=meta)
models[[i]]=fit
coef=summary(fit)$coefficients
coef=coef[grep(variable,rownames(coef)),]

print(head(coef))
rownames(coef)=sub(variable,"",rownames(coef))
coef=lm.pert(meta,form,variable,useBatch=T,dep_var="Topic",ref=reference,verbose=F,numRep=numRep)
results[[i]]=coef
meta=meta_backup
}
if(as.df)
{
ret=c()
for(i in 1:length(results)){
cur=data.frame(results[[i]])
cur["Topic"]=i
#cur["Perturbation"]=rownames(cur)
#cur=cur[,c(6,1,2,3,4,5)]
rownames(cur)=NULL
ret=rbind(ret,cur)
}
#colnames(ret)=c("Perturbation","Effect_Size","SE","t","pval","Topic")
ret["padj"]=p.adjust(ret[,"pval_pert"],"fdr")
ret=ret[order(ret[,"pval_pert"]),]

lst=list()
lst[[1]]=ret
lst[[2]]=models
if(getResid)
{
return(lst)
}

return(lst[[1]])

}


return(results)
}


library(ggplot2)
library(dplyr)
library(tidyr)
##draws the Topic heatmap
drawRes<-function(res,meanNorm=F,addRef=F,reference="GFP")
{


colnames(res)[2]="Log_Odds_Ratio"

if(meanNorm)
{
norm<-res %>% group_by(Topic) %>% summarise(Mean=mean(Log_Odds_Ratio)) %>% as.data.frame()
res<-left_join(res,norm)
res["Log_Odds_Ratio"]=res[,"Log_Odds_Ratio"]-res[,"Mean"]

}

mx=max(abs(res[,2]))+.01
p=ggplot(res,aes(x=factor(Topic),y=perturbation,fill=Log_Odds_Ratio))+geom_tile()+scale_fill_gradient2(low="blue",high="red",mid="white",lim=c(-mx,mx))+xlab("")+xlab("Topic")
return(p)
}



testTotalSTM<-function(res,meta,form)
{
theta=res$theta

theta=log(theta/(1-theta))

theta=data.frame(theta)


num=dim(theta)[1]

meta=cbind(theta,meta)



}

##samples so that evenly distributed between batches
batchSample<-function(seur,pert="GFP",batchs="batch",total=100,perturbation="perturbation")
{
seur<-SubsetData(seur,names(seur@ident)[!is.na(seur@meta.data[,perturbation])])

seur<-SetAllIdent(seur,perturbation)

seur2<-seur



ret<-c()

for(pert in unique(seur2@ident))
{

seur<-SubsetData(seur,WhichCells(seur,pert))

uniq=unique(seur@meta.data[,batchs])
seur@meta.data["samp"]=0

for(batch in uniq)
{
names=names(seur@ident)[seur@meta.data[,batchs]==batch]

seur@meta.data[names,"samp"]=sample(1:length(names))

}

names=names(seur@ident)[order(seur@meta.data[,"samp"])[1:total]]

seur<-seur2

ret<-c(ret,names)
}

return(SubsetData(seur,ret))

}


if(!interactive())
{
set.seed(1)
args=commandArgs(trailingOnly=TRUE)
seur_file=args[1]
outdir=args[2]
numTop=as.numeric(args[3])
print("Read in seurat object!")
load(seur_file)
seur=key
print("Clean up Perts")
#seur<-SubsetData(seur,names(seur@ident)[!is.na(seur@meta.data$perturbation)])

#seur<-SubsetData(seur,names(seur@ident)[!(seur@meta.data$perturbation %in% c("Suv420h1","Map1a","Ankrd11"))])

#seur<-SubsetData(seur,names(seur@ident)[seur@meta.data$Good])
#print("Remove bad batch!")
#seur<-SubsetData(seur,names(seur@ident)[seur@meta.data$batch!="ctx180122"])
#seur<-SubsetData(seur,names(seur@ident)[seur@meta.data$batch!="ctx170723"])

print("Iterate through cell types with enough cells!")
lst=summary(seur@ident)
celltypes=names(lst)[lst>500]

seur@meta.data["nGene"]=scale(seur@meta.data[,"nGene"])

seur@meta.data["batch"]=seur@meta.data[,"orig.ident"]

for(cell in celltypes)
{
print(cell)


#out<-runSTM(seur,c(cell),numberPCs=0,forForm=c("batch","perturbation","nGene"),numTopics=numTop,genes=c(),minPert=10,getMeta=T)

#save(out,file=paste(outdir,"/stm.",as.character(cell),".Robj",sep=""))

load(paste(outdir,"/stm.",as.character(cell),".Robj",sep=""))

print("Test Cell Type!")
res<-testSTM(out[[1]],out[[2]],form="~batch+perturbation+nGene",reference="GFP",variable="perturbation",as.df=T,rawtheta=T,logistic=F,numRep=10000,topics=c(),getResid=F)
save(res,file=paste(outdir,"/test.",as.character(cell),".Robj",sep=""))
print("next!")
}

}


##
##Tests for differences in each coefficient of the STM output
##
testSTM_LMM<-function(res,meta,form="~batch+perturbation+(1|batch:perturbation)",reference="GFP",variable="perturbation",as.df=T,rawtheta=T,logistic=T,numRep=1000,topics=c(),getResid=F)
{
library(lme4)
print("Get theta")
theta=res$theta

if(!rawtheta)
{
theta=thetaPosterior(res,1)

theta <- do.call(rbind, theta)
}
if(length(topics)>0){theta=theta[,topics];theta=data.frame(t(theta));for(i in colnames(theta)){theta[i]=theta[,i]/sum(theta[,i])};theta=t(theta)}

if(logistic)
{
theta=log(theta/(1-theta))
}
print("Combine with metadata")
#theta=data.frame(theta)


numTopics=dim(theta)[2]

meta=cbind(theta,meta)

meta=data.frame(meta)


print("Fix Perturbation")

print(head(meta))
meta["perturbation"]=factor(meta[,"perturbation"])

if(nchar(reference)>0)
{
meta <- within(meta, perturbation <- relevel(perturbation, ref = reference))
}

results=list()

form=as.character(form)
#print(form)
form=paste("Topic",form,sep="")
form=as.formula(form)
print(form)

print("Perform linear regression!")

meta_backup=meta
models=list()
for(i in 1:numTopics)
{
print(paste("Topic:",toString(i)))
colnames(meta)[i]="Topic"
fit=lmer(form,data=meta)
models[[i]]=fit
coef=summary(fit)$coefficients
coef=coef[grep(variable,rownames(coef)),]

print(head(coef))
rownames(coef)=sub(variable,"",rownames(coef))
coef=data.frame(coef)
coef["perturbation"]=rownames(coef)
#coef=lm.pert(meta,form,variable,useBatch=T,dep_var="Topic",ref=reference,verbose=F,numRep=numRep)
results[[i]]=coef
meta=meta_backup
}
if(as.df)
{
ret=c()
for(i in 1:length(results)){
cur=data.frame(results[[i]])
cur["Topic"]=i
cur=cur[,c(dim(cur)[2],1:(dim(cur)[2]-1))]
#cur["Perturbation"]=rownames(cur)
#cur=cur[,c(6,1,2,3,4,5)]
rownames(cur)=NULL
ret=rbind(ret,cur)
}
#colnames(ret)=c("Perturbation","Effect_Size","SE","t","pval","Topic")
#ret["padj"]=p.adjust(ret[,"pval_pert"],"fdr")
#ret=ret[order(ret[,"pval_pert"]),]

lst=list()
ret["Est_P"]<-2 * (1 - pnorm(abs(ret$t.value)))
ret["padj"]=p.adjust(ret[,"Est_P"],"fdr")

ret=ret[order(ret$Est_P),]

lst[[1]]=ret
lst[[2]]=models
if(getResid)
{
return(lst)
}

return(lst[[1]])

}


return(results)
}


library(ggplot2)





##
##Tests for differences in each coefficient of the STM output
##
testSTM_all<-function(res,meta,toUse=c(),form="~batch+perturbation",reference="GFP",variable="perturbation",as.df=T,rawtheta=T,logistic=T,numRep=1000,topics=c(),getResid=F,form1=Topic~batch+perturbation+nUMI,form2=Topic~batch+nUMI,getFit=F)
{
print("Get theta")
theta=res$theta

if(!rawtheta)
{
theta=thetaPosterior(res,1)

theta <- do.call(rbind, theta)
}
if(length(topics)>0){theta=theta[,topics];theta=data.frame(t(theta));for(i in colnames(theta)){theta[i]=theta[,i]/sum(theta[,i])};theta=t(theta)}

if(logistic)
{
theta=log(theta/(1-theta))
}
print("Combine with metadata")
#theta=data.frame(theta)


numTopics=dim(theta)[2]

meta=cbind(theta,meta)

meta=data.frame(meta)


print("Fix Perturbation")

print(head(meta))
if(length(toUse)>0){meta=meta[toUse,]}

meta["perturbation"]=factor(meta[,"perturbation"])

if(nchar(reference)>0)
{
meta <- within(meta, perturbation <- relevel(perturbation, ref = reference))
}

results=list()

#form=as.character(form)
#print(form)
#form=paste("Topic",form,sep="")
#form=as.formula(form)
#print(form)

print("Perform linear regression!")

meta_backup=meta
models=list()
r2<-c()
for(i in 1:numTopics)
{
print(paste("Topic:",toString(i)))
colnames(meta)[i]="Topic"
#form1=Topic~batch+perturbation+nGene
#form2=Topic~batch+nGene
#library(lme4)
fit1=lm(form1,data=meta)
fit2=lm(form2,data=meta)
fit<-anova(fit1,fit2,test="LRT")
#r2<-c(r2,partial.R2(fit1,fit2))
results[[i]]=fit
if(getFit){results[[i]]=fit1}

meta=meta_backup
}


names=c()
pvals<-c()
for(i in 1:length(results)){print(i);pval<-results[[i]]$"Pr(>Chi)"[2];pvals<-c(pvals,pval);names<-c(names,sub("^","Topic",i))}

ret=data.frame(cbind(names,pvals))

colnames(ret)=c("Topic","pvalue")

ret["pvalue"]=as.numeric(as.character(ret[,"pvalue"]))

#ret["r2"]=as.numeric(as.character(ret[,"r2"]))

ret["padj"]=p.adjust(ret[,"pvalue"],"fdr")

return(ret)
}




##
##Tests for differences in each coefficient of the STM output
##
testSTM_pheno<-function(res,meta,toUse=c(),form="~batch+nGene+pheno",pheno="seizure",reference="GFP",variable="perturbation",as.df=T,rawtheta=T,logistic=T,numRep=1000,topics=c(),getResid=F,form1=Topic~batch+perturbation+nGene,form2=Topic~batch+nGene)
{

meta["pheno"]=meta[,pheno]


print("Get theta")
theta=res$theta

if(!rawtheta)
{
theta=thetaPosterior(res,1)

theta <- do.call(rbind, theta)
}
if(length(topics)>0){theta=theta[,topics];theta=data.frame(t(theta));for(i in colnames(theta)){theta[i]=theta[,i]/sum(theta[,i])};theta=t(theta)}

if(logistic)
{
theta=log(theta/(1-theta))
}
print("Combine with metadata")
#theta=data.frame(theta)


numTopics=dim(theta)[2]

meta=cbind(theta,meta)

meta=data.frame(meta)

meta=meta[!is.na(meta[,"pheno"]),]

print(dim(meta))
print(head(meta))

results=list()


print("Perform linear regression!")

meta_backup=meta
models=list()


form=paste("Topic",form,sep="")


for(i in 1:numTopics)
{
print(paste("Topic:",toString(i)))
colnames(meta)[i]="Topic"

print(head(meta))



form=as.formula(form)

print(form)

#form1=Topic~batch+perturbation+nGene
#form2=Topic~batch+nGene
#library(lme4)
fit=lm(form,data=meta)
results[[i]]=fit
meta=meta_backup
}


return(results)
}





GetTopGenes<-function(res,n)
{
lab=labelTopics(res,n=n)


frex=data.frame(lab$frex)
colnames(frex)=sub("^","Gene",as.character(1:dim(frex)[2]))
rownames(frex)=sub("^","Topic",as.character(1:dim(frex)[1]))


n=dim(frex)[2]

top=data.frame(lab$prob)
colnames(top)=sub("^","Gene",as.character(1:n))
rownames(top)=sub("^","Topic",as.character(1:dim(top)[1]))

lift=data.frame(lab$lift)
colnames(lift)=sub("^","Gene",as.character(1:n))
rownames(lift)=sub("^","Topic",as.character(1:dim(lift)[1]))



score=data.frame(lab$score)
colnames(score)=sub("^","Gene",as.character(1:n))
rownames(score)=sub("^","Topic",as.character(1:dim(score)[1]))

res=list()

res[["top"]]=top
res[["lift"]]=lift
res[["score"]]=score
res[["frex"]]=frex

for(i in names(res)){dat=res[[i]];res[[i]]["Topic"]=rownames(dat);n=dim(res[[i]])[2];res[[i]]=res[[i]][,c(n,1:(n-1))]}

return(res)


}


