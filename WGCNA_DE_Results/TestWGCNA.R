source("../CodeDirectory_Current/TestWGNCA_eigen.R",chdir=T)

load("../SeuratObjects/key.500genes.for.analysis.Robj")

load("lists.tidy.current.Robj")

key@meta.data["batch"]=as.character(key@meta.data[,"orig.ident"])

res<-TestWGCNA(key,tab,numRep=10000)

save(res,file="DE.new.Robj")
