library(dplyr)
library(tidyr)


runAllMatrices<-function(lst)
{
for(outdir in lst){print(outdir);makeMatrices(outdir)}
}

makeMatrices<-function(outdir,pert="/psych/genetics_data/ssimmons/Perturb/dialout_3_26_18/pbc_rev.csv")
{
#args = commandArgs(trailingOnly=TRUE)
#outdir=args[1]
#pert=args[2]
savename=paste(outdir,"Matrix.tsv",sep="/")


print("Load and combine data")

withcbc=read.table(paste(outdir,"withcbc.txt",sep="/"),sep="\t",header=F)
withpbc=read.table(paste(outdir,"unique.txt",sep="/"),sep="\t",header=F)

rownames(withpbc)=withpbc[,1]
rownames(withcbc)=withcbc[,1]

colnames(withpbc)=c("Name","Perturbation")
colnames(withcbc)=c("Name","CBC","UMI")

withcbc=withcbc[,2:3]


inter=intersect(rownames(withpbc),rownames(withcbc))

dat=cbind(withpbc[inter,],withcbc[inter,])


print("Clean Up")


dat=dat[grep("^CB:Z:",dat[,"CBC"]),]
dat=dat[grep("^UB:Z:",dat[,"UMI"]),]
dat=dat[grep("^Perturb",dat[,"Perturbation"]),]

dat["CBC"]=sub("^CB:Z:","",dat[,"CBC"])
dat["UMI"]=sub("^UB:Z:","",dat[,"UMI"])


print("Get count data!")
tab <- dat %>% group_by(CBC,UMI,Perturbation) %>% summarise(count=1) %>% group_by(CBC,Perturbation) %>% summarise(count=length(UMI)) %>% spread(Perturbation,count,fill=0)

tab=data.frame(tab)
rownames(tab)=tab[,1]
tab=tab[,2:dim(tab)[2]]

lst=scan(pert,what="")
cols=colnames(tab)
cols=sub("^Perturb","",cols)
cols=as.numeric(cols)
colnames(tab)=lst[cols]


print("save")
write.table(tab,savename,sep="\t",quote=F)


}


