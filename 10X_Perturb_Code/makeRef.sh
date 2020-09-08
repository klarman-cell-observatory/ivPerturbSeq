#!/bin/bash


pert=$1
outdir=$2
outfile=$outdir/comb.pert.fa

seqfive=$(cat fiveprime.fa)
seqthree=$(cat 3prime.fa)

CURLINE=0

rm $outfile

echo "Make fasta reference using perturbed barcodes!"

while read pert; do
	CURLINE=$((CURLINE+1))
	echo ">Perturb$CURLINE" >> $outfile
	echo "$seqfive$pert$seqthree" >> $outfile
done < $pert

echo "Make STAR reference from fasta file!"

mkdir $outdir/STAR

STAR --runMode genomeGenerate --genomeDir $outdir/STAR --genomeSAindexNbases 6 --genomeFastaFiles $outfile

echo "Done making reference!"

