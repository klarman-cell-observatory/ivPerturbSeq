#! /bin/bash

#$ -cwd
#$ -q broad
#$ -P regevlab
#$ -e err.err
#$ -o out.log
#$ -l h_vmem=90g,h_rt=72:00:00,os=RedHat6
#$ -t 1-18

source /broad/software/scripts/useuse
reuse UGER
use .python-2.7.9-sqlite3-rtrees

SEEDFILE=/psych/genetics_data/ssimmons/Perturb/Dialout_code/AllFastqs/sample.sheet.Combined.txt

statfile="/broad/hptmp/ssimmons/Perturb/dialout/stats.combined.txt"

fq1_list=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $1}')
fq2_list=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $2}')
barcodes=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $3}')
output=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $4}')
primer=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $5}')

pert=/psych/genetics_data/ssimmons/Perturb/Dialout_code/perturbations.txt



fq1=${output}.R1.fastq
fq2=${output}.R2.fastq

echo Combine to fastq!

fq1_list=$(echo $fq1_list | sed 's/,/ /g')
fq2_list=$(echo $fq2_list | sed 's/,/ /g')

cat $fq1_list > $fq1
cat $fq2_list > $fq2

echo Extract info!

python /psych/genetics_data/ssimmons/Perturb/Dialout_code/Perturb_mat_filt.py $fq1 $fq2 $pert $barcodes $output $primer &> ${output}.err

echo $output $(wc -l $output) >> $statfile

echo Clean up!
rm $fq1
rm $fq2
