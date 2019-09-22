library(openxlsx)
library(stm)


lst<-system("ls -1 Results5Topics_cur2/stm.*.Robj", inter=T)


ret=list()

for(l in lst)
{
print(l)
load(l)

celltype=sub(".Robj","",sub("Results5Topics_cur2/stm.","",l,fixed=T),fixed=T)

print(celltype)

res=out[[1]]

mat= t(res$beta$logbeta[[1]])
mat=data.frame(mat)
colnames(mat)=sub("^","Topic",1:5)

mat["Gene"]=res$vocab

mat=mat[,c(6,1:5)]

ret[[paste("Gene Weights_",celltype,sep="")]]=mat

mat=res$theta

mat=data.frame(mat)

colnames(mat)=sub("^","Topic",1:5)

mat["Cell"]=rownames(t(res$mu$mu))

mat=mat[,c(6,1:5)]

ret[[paste("Cell Weights_",celltype,sep="")]]=mat

}

write.xlsx(ret,"Model.Params.xlsx")
