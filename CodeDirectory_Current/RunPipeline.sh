#! /bin/bash

#$ -cwd
#$ -q broad
#$ -P regevlab

#$ -e err.7.err
#$ -o out.7.log

#$ -l h_vmem=90g,h_rt=60:00:00,os="RedHat7"

source /broad/software/scripts/useuse
use R-3.5
use .samtools-1.8
use UGER
use .python-2.7.9-sqlite3-rtrees
use .java-jdk-1.8.0_181-x86-64
use .hdf5-1.8.16


export R_LIBS_USER=/psych/genetics_data/ssimmons/R/x86_64-pc-linux-gnu-library/3.5_Pert


echo "Start by clustering and filtering the matrix, loading Perturbation information, etc"
cd ../SeurObj_R3.5
Rscript GetData.n500.start.R 
Rscript GetData.n500.end.R

#echo "Get QC information"
#cd ../QCPlots
#Rscript QCPlots.R ##Makes some simple QC plots

echo "Perform Cell Type composition comparision"
cd ../CellComposition
Rscript AdamsMethod.R ##Fits a poisson regression model to see which perturbations effect which cell type perturbations


#echo "Cluster individual cell types and perform basic QC analysis on them"
#cd ../CellType_Figs
#Rscript ClusterCellTypes.R 	


echo "Next, perform WGCNA"
cd ../WGCNA_DE_Results
#Generate WGCNA and plot (need stuff from Xin)
Rscript TestWGCNA.R #Test WGCNA (easy)

echo "Next, perform STM analysis"
cd ../UpdatedSTM
mkdir Results5Topics
Rscript ../CodeDirectory_Current/RunSTM.R ../SeurObj_R3.5/key.500genes.for.analysis.Robj Results5Topics	5 ##Run STM with 5 topics, using perturbation and batch info tp fit model 
cp MakeImages.R Results5Topics/MakeImages.R
cd Results5Topics
Rscript MakeImages.R ##Make images
cd ..

echo "Check WGCNA with other DS schemes"
cd ../DownSampleCheck
Rscript DS.WGCNA.Check.R

#echo "Perform gene level analysis"
#cd ../RemoveCellTypeGenes
#Rscript GetCellTypeGenes_correct_subcluster.R ##gets genes most effected by perturbations

#echo "Perform analysis on known gene sets"
#cd ../Outside_Gene_Sets
#Rscript TestAll.anova.R ##Tests to see which GO terms are most effected



