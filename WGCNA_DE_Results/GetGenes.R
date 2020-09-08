library(WGCNA)
library(Seurat)




##
##Basic code used for WGCNA 
##
RunWGCNA<-function(seur,output,myPower)
{
	print("Filter Data!")
	dat=as.matrix(seur@data)
	 dat<-dat[!(rowMeans(dat, na.rm = TRUE) < 0.15),] #get rid of lowly expressed genes
	 dat<-dat[!(rowMeans(dat, na.rm = TRUE) > 9),] # get rid of highly expressed genes
	
	dat <- dat[!(rownames(dat) %in% grep(pattern = "^Rp", x = rownames(dat), value = TRUE)),] ##Remove ribosomal genes
	
	dat=t(dat)
	
	#save(colnames(dat),file=paste(output,".genes.Robj",sep=""))

	return(colnames(dat))

	print("Remaining data:")
	print(head(colnames(dat)))
	print(head(rownames(dat)))
	
	print(dim(dat))
	
	print("Get modules!")
	
	
	networkType = "signed"
	minModuleSize = 7 
	mergeCutHeight = 0.15
	net = blockwiseModules(
    # Input data
    dat,
    
    # Data checking options
    checkMissingData = TRUE,
    
    # Options for splitting data into blocks
    blocks = NULL,
    maxBlockSize = 5000,   # xin jin change this from 5000 to 2000
    randomSeed = 59069,
    
    # load TOM from previously saved file?
    loadTOM = FALSE,
    
    # Network construction arguments: correlation options
    corType = "bicor", # more robust for non-normal? Song 2012
    maxPOutliers = 0.1,    # XinJin changed from dedault 0.9
    quickCor = 0,
    pearsonFallback = "individual",
    cosineCorrelation = FALSE,
    
    # Adjacency function options
    # this is where power gets input
    power = myPower,
    networkType = networkType,
    
    # Topological overlap options
    TOMType = "signed",
    TOMDenom = "min",
    
    # Saving or returning TOM
    getTOMs = NULL,
    saveTOMs = TRUE,
    saveTOMFileBase = "blockwiseTOM",
    
    # Basic tree cut options
    deepSplit = 3,   # xin jin changed from 2 to 3 to break down module size
    detectCutHeight = 0.995,
    minModuleSize = minModuleSize,
    
    # Advanced tree cut options
    maxCoreScatter = NULL, minGap = NULL,
    maxAbsCoreScatter = NULL, minAbsGap = NULL,
    minSplitHeight = NULL, minAbsSplitHeight = NULL,
    useBranchEigennodeDissim = FALSE,
    minBranchEigennodeDissim = mergeCutHeight,
    pamStage = TRUE, pamRespectsDendro = TRUE,
    
    # Gene reassignment, module trimming, and module "significance" criteria
    reassignThreshold = 1e-6,
    minCoreKME = 0.5,
    minCoreKMESize = minModuleSize/3,
    minKMEtoStay = 0.3,
    
    # Module merging options
    mergeCutHeight = 0.15,
    impute = TRUE,
    trapErrors = TRUE,
    
    # Output options
    numericLabels = FALSE,
    
    # Options controlling behaviour
    nThreads = 0,
    verbose = 3, 
    indent = 0)
    
	mods=NULL
	
	print("Save modules!")
	save(net,file=output)
	
	print("Done!")
	
	return(net)
	
}


if(!interactive())
{
	print("Load data and get arguments")
	load("../SeuratObjects/key.500genes.Robj")
	args = commandArgs(trailingOnly=TRUE)
	
	celltype=args[1]
	seur<-SubsetData(key,WhichCells(key,celltype))
	
	myPower=as.numeric(args[2])
	
	savename=paste("Results/Genes.",celltype,".Robj",sep="")
	
	print("Run WGCNA!")
	
	net<-RunWGCNA(seur,savename,myPower)
	save(net,file=savename)	
	print("Done!")
}

