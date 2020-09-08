#! /bin/bash

##########################################
##Extracting perturbation barcodes from 10X
#########################################
##Takes in a 10X directory (the directory outout by the cellrange counts function (at least in cellranger 2.1.0) one level about the 'outs' directory)
##A list of perturbation barcodes (be careful of not getting reverse compliment
##the out directory (note each sample requires its own outdirectory)
##
##It will then extract unmapped 10X reads, map them to BFP+available perturbation barcodes+trailing sequence, extract uniquely mapped reads, and turn them into a cell by perturbation UMI based count matrix 
##
##Required to be installed/ on the path:
##1) samtools
##2) bedtools
##3) STAR aligner
##4) R, including dplyr and tidyr
##########################################


tenXdir=$1
pert=$2
outdir=$3

tenXbam=$tenXdir/outs/possorted_genome_bam.bam

echo "Setup!"

mkdir $outdir

echo "Get unmapped reads!"
samtools view -f 4 $tenXbam | grep UB | grep CB | grep 'AGAATT\|CCTAGA'  | samtools view -S -b  > $outdir/unmapped.bam

echo "Make into fastq"
bedtools bamtofastq -i $outdir/unmapped.bam -fq $outdir/unmapped.fq


echo "Make STAR Reference"
makeRef.sh $pert $outdir 

echo "Map to reference!"
STAR --genomeDir $outdir/STAR --readFilesIn $outdir/unmapped.fq --outFileNamePrefix $outdir/STAR_out --outSAMtype BAM SortedByCoordinate --alignIntronMax 1

echo "Clean up fastq file"
#$outdir/unmapped.fq

echo "Get uniquelly mapped + high quality reads"
#samtools view -F 2820 $outdir/STAR_outAligned.sortedByCoord.out.bam  | awk '$5==255' | awk '$4>655' | awk '$4<714'| sed 's/AS:i://g'| awk '$14 > 45' | awk '{print $1"\t"$3}' >  $outdir/unique.txt
samtools view -F 2820 $outdir/STAR_outAligned.sortedByCoord.out.bam  | awk '$5==255' | awk '$4>655' | awk '$4<714'| sed 's/nM:i://g'| awk '$15 < 3' | awk '{print $1"\t"$3}' >  $outdir/unique.txt
samtools view -F 2820 $outdir/STAR_outAligned.sortedByCoord.out.bam -b  > $outdir/keep.bam.tmp



##This next step is required to get the CBC and UMI from the original BAM file, since does not get retained in mapping steps
echo "Get mapped reads"
awk '{print $1}' $outdir/unique.txt > $outdir/mapped.txt 
samtools view $outdir/unmapped.bam | fgrep -w -f $outdir/mapped.txt | awk '{print $1"\t"$21"\t"$24}' > $outdir/withcbc.txt

echo "Extract count matrix!"
Rscript GetMatrix.R $outdir $pert

###QUESTION: Do we want to clean up everything? Only some?
echo Clean up!
#rm $outdir/*bam
#rm $outdir/*sam
#rm $outdir/*txt
#rm $outdir/*fq
#rm $outdir/*fa
#rm -r $outdir/STAR*
#mv $outdir/keep.bam.tmp $outdir/keep.bam








