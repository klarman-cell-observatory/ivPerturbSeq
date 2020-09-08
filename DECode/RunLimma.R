library(Seurat)
library(Matrix);
library(limma)
library(edgeR)

##Runs limma with gseurat object in given celltype
##usese given formula for design and one to be tested
runLimma<-function(seur,celltype,form=~0+batch+perturbation,test="perturbation",ref="GFP")
{

print("Get count data!")
dat=seur@raw.data[,names(seur@ident)[seur@ident==celltype]]

meta=seur@meta.data[seur@ident==celltype,]

print("Filter Genes")
mn=Matrix::rowMeans(dat>0)
dat=dat[mn>.05,]
dat=as.matrix(dat)


print("Set up for testing!")
dge <- DGEList(dat)

dge<-calcNormFactors(dge)

meta[test]=factor(meta[,test])

if(!(ref %in% meta[,test])){return(c())}

meta[test]=relevel(meta[,test],ref)

##Build design matrix!
design <- model.matrix(form,data=meta)

print(head(design))


print("Fit Matrix!")
v <- voom(dge, design, plot = F)
fit <- lmFit(v, design)

print("Test")

pert=unique(meta[,test])
pert=pert[pert!=ref]
out<-lapply(pert,
function(x)
{
print(x);
toTest=paste(test,x,"-",test,ref,sep="")
print(toTest)
#toTest=as.formula(toTest)
print(toTest)

ret<-tryCatch({
contr <- makeContrasts(contrasts=toTest, levels = colnames(coef(fit)))



mrk <- contrasts.fit(fit, contr)
mrk<-eBayes(mrk)

top.table <- topTable(mrk, sort.by = "P", n = Inf)

return(top.table)
},
error = function(e){print("Yuck");return(c())}
)

return(ret)


})

names(out)=pert

lst<-c()

print("Combine")

for(i in names(out)){if(is.data.frame(out[[i]])){lst<-c(lst,i)}}

out<-out[lst]

for(i in names(out)){out[[i]]["comp"]=i;out[[i]]["Gene"]=rownames(out[[i]])}


out=do.call(rbind,out)

#print(head(out))


print("Reorder")

out=out[order(out$P.Value),]

out["Padj"]=p.adjust(out[,"P.Value"],"fdr")

return(out)

}


##Run Limma for all cell types!
runLimma_all<-function(seur,form=~0+perturbation+cngeneson+batch,test="perturbation",ref="GFP")
{

seur@meta.data["cngeneson"]=scale(seur@meta.data[,"nGene"])

out<-lapply(unique(seur@ident),function(x){print(x);if(sum(seur@ident==x)<50){return(c())};mrk=runLimma(seur,x,form,test,ref);mrk=data.frame(mrk);mrk["CellType"]=x;return(mrk)})

#mrk=do.call(rbind,out)

return(out)

}




