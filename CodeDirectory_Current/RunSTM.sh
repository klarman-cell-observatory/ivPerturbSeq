#! /bin/bash

#$ -cwd
#$ -q broad
#$ -P regevlab

#$ -e err.7.err
#$ -o out.7.log

#$ -l h_vmem=90g,h_rt=60:00:00,os="RedHat7"

source /broad/software/scripts/useuse
use R-3.5

#Seurat object
echo $1
#outdir
echo $2
#Number of topic
echo $3

export R_LIBS_USER=/psych/genetics_data/ssimmons/R/x86_64-pc-linux-gnu-library/3.5_scumi


Rscript /psych/genetics_data/ssimmons/Perturb/Updated_Data_March_11_2019/CodeDirectory_Current/RunSTM.R $1 $2 $3
#Rscript /psych/genetics_data/ssimmons/Perturb/BasicPipeline/RunSTM.R $1 $2 $3 
