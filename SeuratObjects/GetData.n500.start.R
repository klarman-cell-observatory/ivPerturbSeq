set.seed(1)
savedir="."
codedir="../CodeDirectory_Current"

source(paste(codedir,"/load_Seurat.R",sep=""),chdir=T)
source(paste(codedir,"/SharedVariable.R",sep=""),chdir=T)
source(paste(codedir,"/niceDraw.R",sep=""),chdir=T)

load("dat.Robj")


seur<-dir10X(dat=dat)


#Basic processing
seur<-SharedVariable(seur)

seur<-ScaleData(seur,genes.use=seur@var.genes,vars.to.regress="nUMI")

seur<-RunPCA(seur,pcs.compute=40)


##
##We found different version of R led to difference irlba results, so here we include the ones we got originally
#load("pc.3.3.Robj")

#seur@dr$pca@cell.embeddings=pc_results

pdf(paste(savedir,"/ElbowPlot.pdf",sep=""))
PCElbowPlot(seur,40)
dev.off()


#TSNE and clustering
seur<-RunTSNE(seur,dims.use=1:26,verbose=T,pca=F)

seur<-FindClusters(seur,dims.use=1:26,resolution=1.2)

#save(seur,file="temp.Robj")

##Shift cluster numbering so starts at 1
seur@meta.data["Clustering"]=as.numeric(as.character(seur@ident))+1
seur<-SetAllIdent(seur,"Clustering")

##Plot Clustering
pdf(paste(savedir,"/TSNE.clustering.all.pdf",sep=""))
TSNEPlot(seur,T)
dev.off()



##markers are identified by DropViz and DE between clusters in this dataset (not shown)
pdf(paste(savedir,"/FeaturePlot.markers.pdf",sep=""))
markers=c("Slc17a7","Slc17a6","Gad2","Gad1","Plp1","Mobp","Pdgfra","Top2a","Slc1a3","Flt1","Vtn","Bgn","Slc47a1","Csf1r","Eomes","Slc4a5","Dnah12","Neurod6","Ndnf","Slc22a6","Aqp4","Pllp")
celltypes=c("Excitatory","Excitatory","Inhibitory","Inhibitory","Oligo","Mature Oligo","OPC","Cycling","Astroglia","Vascular","Vascular","Vascular","Fibroblast-like","Microglia","IPC","Choroid Plexus","Ependyma","Cortical Excitatory","Hippocampus","Fibroblast-like","Astroglia","ODC")


source(paste(codedir,"/niceDraw.R",sep=""))

for(i in 1:length(markers))
{
print(markers[i])
p=niceFeaturePlot(seur,markers[i])+ggtitle(paste(celltypes[i],markers[i],sep=": "))
plot(p)
}

dev.off()

##get Scrublet
writeMM(seur@raw.data[,names(seur@ident)],paste(savedir,"/MM.for.scrublet.txt",sep=""))

#RUN SCRUBLET IN PYTHON, see script ....

system(paste("python ",codedir,"/","RunScrublet.py ",savedir,"/MM.for.scrublet.txt ",savedir,"/scrublet.txt",sep=""))

#Actual Code Used## scrub<-scan(paste(savedir,"/scrublet.txt",sep=""))
scrub<-scan(paste(savedir,"/scrublet.init.txt",sep="")) ## For Reproducible results



seur@meta.data["scrublet"]=scrub

p=niceFeaturePlot(seur,"scrublet",datainfo=T)
ggsave(paste(savedir,"/scrublet.pdf",sep=""),p)

save(seur,file=paste(savedir,"/all.pre.celltype.Robj",sep=""))

#q(save="no")
load("mapping.Robj")
##get a list where the ith entry is the cell type label of cluster i
mapToCellType<-rep("None",38)


#mapToCellType[c(2,10,14,16,34,8)]="Excitatory"
#mapToCellType[c(11,15,22,24,29)]="Inhibitory"
#mapToCellType[c(24)]="Inhibitory"

#mapToCellType[c(36)]="IPC"
#mapToCellType[c(3,30,25)]="ODC"

#mapToCellType[c(23,12,4,5,6,9)]="Astroglia"
#mapToCellType[c(17,35,20)]="Fibroblast-like"
#mapToCellType[c(33)]="Ependyma"
#mapToCellType[c(26)]="Hippocampus"
#mapToCellType[c(32)]="Choroid Plexus"
#mapToCellType[c(28,7,18)]="Vascular"
#mapToCellType[c(38)]="Doublet"
#mapToCellType[c(1,13,19,31,21,37,27)]="Microglia"
#save(seur,file=paste(savedir,"/all.500genes.Robj",sep=""))
mapToCellType=lst


seur@meta.data["CellType"]=mapToCellType[as.numeric(as.character(seur@ident))]

##save result
save(seur,file=paste(savedir,"/all.500genes.Robj",sep=""))

##Take out mixed cluster and subcluster to get better labelling
seur<-SetAllIdent(seur,"CellType")

#mix<-SubsetData(seur,WhichCells(seur,"Mixed"))

#mix<-SharedVariable(mix)

#mix<-ScaleData(mix,genes.use=mix@var.genes,vars.to.regress=c("nUMI","orig.ident"))


#mix<-RunTSNE(mix,dims.use=1:10,verbose=T,pca=F)

#mix<-FindClusters(mix,dims.use=1:10,force.recalc=T)

#mapCellTypeMixed=c()


#mapCellTypeMixed=rep("Inhibitory",8)

#mapCellTypeMixed[c(0,1,7,4)+1]="Astroglia" 



#mix@meta.data["CellType_new"]=mapCellTypeMixed[as.numeric(as.character(mix@ident))+1]

#mix<-SetAllIdent(mix,"CellType_new")

#seur@meta.data["CellType_initial"]=seur@meta.data[,"CellType"]

#seur@meta.data[names(mix@ident),"CellType"]=as.character(mix@ident)

#save(mix,file=paste(savedir,"/mixed.500genes.Robj",sep=""))

##Figure out scrublet cutoff
p=ggplot(seur@meta.data,aes(x=scrublet))+geom_histogram(bins=100)+geom_vline(xintercept=.15)

ggsave(paste(savedir,"/histogram.scrublet.pdf",sep=""),p)

tab<-seur@meta.data %>% group_by(orig.ident) %>% summarise(Num=length(nGene),Expected=Num/100000,actual=mean(scrublet>.15)) %>% as.data.frame()

p=ggplot(tab,aes(x=Expected,y=actual))+geom_point()+geom_abline(intercept=0,slope=1)+xlab("Expected Doublet Rate")+ylab("Scrublet Doublet Rate")


ggsave(paste(savedir,"/Expected.vs.Actual.Rate.Scrublet.pdf",sep=""),p)

seur<-SubsetData(seur,names(seur@ident)[seur@meta.data$scrublet<.15])




seur<-SubsetData(seur,names(seur@ident)[seur@ident!="Doublet"])



save(seur,file=paste(savedir,"/cleaned.500genes.Robj",sep=""))

key<-SubsetData(seur,WhichCells(seur,unique(seur@ident)[c(1,3,4,5,6)]))

cutoffs<-cbind(c("Astroglia","Excitatory","Inhibitory","Microglia","ODC"),c(950,2000,500,700,600))
cutoffs=data.frame(cutoffs)
cutoffs[2]=as.numeric(as.character(cutoffs[,2]))
colnames(cutoffs)=c("CellType","Cutoff")

p=ggplot(key@meta.data,aes(x=nGene,fill=CellType))+geom_histogram(bins=100)+facet_wrap(~CellType,ncol=1)+geom_vline(data=cutoffs,aes(xintercept=Cutoff))

ggsave(paste(savedir,"/cutoffs.nGene.new.pdf",sep=""),p)

rownames(cutoffs)=cutoffs[,1]
key@meta.data[,"Cutoff"]=cutoffs[as.character(key@ident),"Cutoff"]

key@meta.data["Keep"]=key@meta.data$nGene>key@meta.data$Cutoff

key<-SubsetData(key,names(key@ident)[key@meta.data[,"Keep"]])


key<-SubsetData(key,names(key@ident)[key@meta.data$orig.ident!=17])

pdf(paste(savedir,"/TSNE.key.pdf",sep=""))
TSNEPlot(key,T)
dev.off()

##Add dialout perturbations!

dialoutdir="../Dialout_Results"
load(paste(dialoutdir,"perturbation.table.Robj",sep="/"))
source(paste(dialoutdir,"GetPerturbations.R",sep="/"),chdir=T)

key<-getPerts(dat,"Name","gene","Count",1.3,seur=key,saveAs="perturbation_dialout")


load(paste(dialoutdir,"perturbation.10X.table.Robj",sep="/"))

#key<-getPerts(dat,"Name","gene","Count",1.3,seur=key,saveAs="perturbation_10X",revComp=T)

rownames(tab)=tab[,"Name"]

inter=intersect(names(key@ident),rownames(tab))

key@meta.data[inter,"perturbation_10X"]=tab[inter,"gene"]

key@meta.data[inter,"perturbation_10X_Raio"]=tab[inter,"Ratio"]

key@meta.data[inter,"perturbation_10X_Count"]=tab[inter,"Max_Count"]

key@meta.data["perturbation"]=key@meta.data[,"perturbation_dialout"]

key@meta.data[is.na(key@meta.data["perturbation"]),"perturbation"]=key@meta.data[is.na(key@meta.data["perturbation"]),"perturbation_10X"]


key@meta.data[!is.na(key@meta.data[,"perturbation_10X"]) & !is.na(key@meta.data[,"perturbation_dialout"]) & key@meta.data[,"perturbation_10X"]!=key@meta.data[,"perturbation_dialout"] ,"perturbation"]=NA

key@meta.data["batch"]=key@meta.data[,"orig.ident"]
save(key,file=paste(savedir,"/key.500genes.Robj",sep=""))

##Iterate over all cell types to get initial subclustering, in a list of Seurat objects
seurs<-list()

 for(i in unique(key@ident)){cur<-SubsetData(key,WhichCells(key,i));cur<-SharedVariable(cur);cur<-ScaleData(cur,genes.use=cur@var.genes,vars.to.regress=c("nUMI","batch"));cur<-RunPCA(cur);cur<-RunTSNE(cur,dims.use=1:15,verbose=T,pca=F);cur<-FindClusters(cur,dims.use=1:15);seurs[[i]]=cur;print(i)}



##Get Subset For Actual Analysis
PerturbCode=codedir
source(paste(PerturbCode,"BatchSample.v2.R",sep=""),chdir=T)

key@meta.data["batch"]=key@meta.data[,"orig.ident"]


key<-SubsetData(key,names(key@ident)[!is.na(key@meta.data$perturbation)])

key<-SubsetData(key,names(key@ident)[!(key@meta.data$perturbation %in% c("Suv420h1","Map1a","Ankrd11"))])

keep<-c()

for(celltype in unique(key@ident))
{
cells<-SubsetData(key,WhichCells(key,celltype))

cells<-SetAllIdent(cells,"perturbation")

lst<-table(as.character(cells@ident))


cells<-SubsetData(cells,names(cells@ident)[cells@ident %in% names(lst)[lst>10]])


lst<-table(as.character(cells@ident))

med=2*median(lst)

cells<-batchSample(cells,total=med,maxPerc=1.0)




lst<-table(as.character(cells@ident))


cells<-SubsetData(cells,names(cells@ident)[cells@ident %in% names(lst)[lst>10]])

keep<-c(keep,names(cells@ident))

}


key<-SetAllIdent(key,keep)






save(key,file=paste(savedir,"/key.500genes.for.analysis.Robj",sep=""))
