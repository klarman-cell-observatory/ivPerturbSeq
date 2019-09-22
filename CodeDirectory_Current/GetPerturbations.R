library(dplyr)
library(tidyr)
library(Biostrings)


##
##Given a tidyr dataframe (dat), with one lane for cell name (NameColumn), one lane for PBC (PertColumn), and one lane
##for UMI Counts (CountColumn), gets the perturbation for each cell with the largest number of UMIs. 
##Only considers cells where the most common perturbation divided by the next most common perturbation is > a given 
##quantity (Ratio, 1.4 by default)
##

getPerts<-function(dat,NameColumn,PertColumn,CountColumn,minRatio=1.4,revComp=F,seur=NA,saveAs="")
{
print("Clean Data")
dat=dat[,c(NameColumn,PertColumn,CountColumn)]

#print(head(dat))

colnames(dat)=c("Name","PBC","Count")


if(revComp)
{
print("Get Revers!")
dat["PBC"]=RevComp(dat[,"PBC"])
}

tab <- dat %>% spread(PBC,Count,fill=0)


#print(head(tab))

ratio<-apply(tab[,2:dim(tab)[2]],1,function(x){x=x[order(x,decreasing=T)];if(x[2]==0){return(as.numeric("inf"))};rat=x[1]/x[2];return(rat)})

pbc<-apply(tab[,2:dim(tab)[2]],1,function(x){return(colnames(tab[,2:dim(tab)[2]])[which.max(x)])})

max_count=apply(tab[,2:dim(tab)[2]],1,max)

res<-cbind(tab[,1],pbc)

res=data.frame(res)

colnames(res)=c("Name","PBC")

res["Ratio"]=ratio

res["Max_Count"]=max_count

res=res[res$Ratio>minRatio,]

print(head(res))

if(nchar(saveAs)>0)
{
seur@meta.data[saveAs]=NA

rownames(res)=res[,1]

inter=intersect(rownames(res),names(seur@ident))

seur@meta.data[inter,saveAs]=as.character(res[inter,"PBC"])
print(class(res[inter,"Ratio"]))
seur@meta.data[inter,paste(saveAs,"_Ratio",sep="")]=res[inter,"Ratio"]

print(class(res[inter,"Max_Count"]))
seur@meta.data[inter,paste(saveAs,"_Count",sep="")]=res[inter,"Max_Count"]

return(seur)
}

return(res)

}




RevComp<-function(dna)
{

dna<-reverseComplement(DNAStringSet(dna))
dna=as.character(dna)

return(dna)
}

