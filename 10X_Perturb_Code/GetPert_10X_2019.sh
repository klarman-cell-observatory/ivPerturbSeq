#! /bin/bash

#$ -cwd
#$ -q broad
#$ -P regevlab
#$ -e err.err
#$ -o out.log
#$ -l h_vmem=90g,h_rt=24:00:00,os=RedHat6
#$ -t 1-2

source /broad/software/scripts/useuse
reuse UGER
reuse .samtools-1.3.1
reuse BEDTools
reuse R-3.3
reuse .star-2.5.2b

export R_LIBS_USER=/psych/genetics_data/ssimmons/R/x86_64-pc-linux-gnu-library/3.3

SEEDFILE=$1
indir=$2
outdir=$3
pertfile=$4

tenXdir=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $1}')
name=$tenXdir

tenXdir=$indir/$tenXdir

$pertfile=/psych/genetics_data/ssimmons/Perturb/dialout_3_26_18/pbc_rev.csv

outdir=$3/$name

mkdir $outdir

echo $tenXdir
echo $pertfile
echo $outdir
CurDir=$(dirname "${BASH_SOURCE[0]}")

${CurDir}/mapIt.sh $tenXdir $pertfile $outdir


