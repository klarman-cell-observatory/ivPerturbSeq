


##
##Performs perurbations to get p-values for lm
##
lm.pert<-function(dat,form,perturbations="perturbation",batch="batch",numRep=1000,useBatch=F,dep_var="DPT",ref="",minPert=10,verbose=T)
{

dat=dat[!is.na(dat[,perturbations]),]

lst=table(dat[,perturbations])

dat=dat[dat[,perturbations] %in% names(lst)[lst>minPert],]

if(nchar(ref)>1)
{
dat[perturbations]=factor(dat[,perturbations])
dat <- within(dat, perturbation <- relevel(perturbation, ref = ref))
}
dat[dep_var]=scale(dat[,dep_var])


mat=dat
fit<-lm(form,data=dat)


tvals=abs(summary(fit)$coefficients[,3])


res=data.frame(cbind(tvals))

for(i in 1:numRep){
if(verbose)
{
if(i%%50 ==0)
{
print(i)
}
}
if(!useBatch)
{
dat[perturbations]=sample(dat[,perturbations])
}
else
{
batchs=unique(dat[,batch])
for(i in batchs)
{
dat[dat[,batch]==i,perturbations]=sample(dat[dat[,batch]==i,perturbations])
}
}

fit<-lm(form,data=dat)

tvals=abs(summary(fit)$coefficients[,3])

res<-cbind(res,tvals)

}


pv=apply(res,1,function(x){sum(x>=x[1])/length(x)})

fit<-lm(form,data=mat)
res=summary(fit)$coefficients

res=data.frame(res)
res["pert_pval"]=pv

colnames(res)=c("Effect_Size","SD","T","pval","pval_pert")

res=res[grep(perturbations,rownames(res)),]

res["perturbation"]=sub(perturbations,"",rownames(res))

res=res[order(res[,"pval_pert"]),]

res["padj"]=p.adjust(res[,"pval_pert"],"fdr")

res=res[,c(6,1:6,7)]

return(res)

}









##
##Performs glm
##
glm.friendly<-function(dat,form,perturbations="perturbation",batch="batch",numRep=1000,useBatch=F,dep_var="DPT",ref="",minPert=10,verbose=T)
{

dat=dat[!is.na(dat[,perturbations]),]

lst=table(dat[,perturbations])

dat=dat[dat[,perturbations] %in% names(lst)[lst>minPert],]

if(nchar(ref)>1)
{
dat[perturbations]=factor(dat[,perturbations])
dat <- within(dat, perturbation <- relevel(perturbation, ref = ref))
}
#dat[dep_var]=scale(dat[,dep_var])


mat=dat
fit<-glm(form,data=dat,family="binomial")
fit<-summary(fit)
fit<-fit$coefficients
fit=fit[grep("pert",rownames(fit)),]
fit=fit[order(fit[,4]),]
rownames(fit)=sub("perturbation","",rownames(fit))
fit=data.frame(fit)
fit["padj"]=p.adjust(fit[,4],"fdr")
return(fit)

}

