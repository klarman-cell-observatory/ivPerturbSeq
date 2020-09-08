#! /bin/bash

#$ -cwd
#$ -q broad
#$ -P regevlab
#$ -e err.err
#$ -o out.log
#$ -l h_vmem=90g,h_rt=72:00:00,os=RedHat6

source /broad/software/scripts/useuse
reuse UGER
use .python-2.7.9-sqlite3-rtrees

SEEDFILE=$1
statfile="stats.combined.txt"

fq1_list=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $1}') ##read1 fastq in comma separated list
fq2_list=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $2}') ##read2 fastq in comma separated list
barcodes=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $3}') ##Cell barcode in list, one per line
output=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $4}')	
primer=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $5}')	##which primer was used in this experiment

pert=perturbations.txt ##list of perturbations



fq1=${output}.R1.fastq
fq2=${output}.R2.fastq

echo Combine to fastq!

fq1_list=$(echo $fq1_list | sed 's/,/ /g')
fq2_list=$(echo $fq2_list | sed 's/,/ /g')

cat $fq1_list > $fq1
cat $fq2_list > $fq2

echo Extract info!

python Perturb_mat_filt.py $fq1 $fq2 $pert $barcodes $output $primer &> ${output}.err

echo $output $(wc -l $output) >> $statfile

echo Clean up!
rm $fq1
rm $fq2
