##Relies on:
#R
#samtools-1.8
#UGER
#python 2
#java
#hdf5


echo "Start by clustering and filtering the matrix, loading Perturbation information, etc"
cd ../SeuratObjects
Rscript GetData.n500.start.R 
Rscript GetData.n500.end.R


echo "Perform Cell Type composition comparision"
cd ../CellComposition
Rscript AdamsMethod.R ##Fits a poisson regression model to see which perturbations effect which cell type perturbations


echo "Next, perform WGCNA"
cd ../WGCNA_DE_Results
#Generate WGCNA and plot (need stuff from Xin)
Rscript GetGenes.R
Rscript TestWGCNA.R #Test WGCNA (easy)

echo "Next, perform STM analysis"
cd ../UpdatedSTM
mkdir Results5Topics
Rscript ../CodeDirectory_Current/RunSTM.R ../SeuratObjects/key.500genes.for.analysis.Robj Results5Topics	5 ##Run STM with 5 topics, using perturbation and batch info tp fit model 
cp MakeImages.R Results5Topics/MakeImages.R
cd Results5Topics
Rscript MakeImages.R ##Make images
cd ..





