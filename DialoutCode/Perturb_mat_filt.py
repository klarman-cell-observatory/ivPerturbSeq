import sys
import re
import pandas as pd
import numpy as np
import Levenshtein

##Reverse compliment!
def RevComp(dna):
	dna_dict={"A":"T","T":"A","C":"G","G":"C","N":"N"}
	new_dna=""
	for base in dna:
		new_dna=dna_dict[base]+new_dna
	return(new_dna)


##Takes in pointer to read file
##produces pandas data frame for it
def flattenFastq(read):
	mat=pd.read_csv(read,sep=' ',header=None)
	mat1=mat[1::4]
	mat2=mat[0::4]
	mat1.index=mat2.index
	mat1=pd.concat([mat2,mat1],axis=1)
	mat1.columns=['name','i7','seq','empty']
	
	mat1=mat1[['name','seq']]

	mat1=mat1.set_index('name')


	return(mat1)



##Filters based on TPT
def filterTPT(Reads,cutoff=.2):
	index=Reads.index
	Reads_input=Reads
	Reads["counts"]=[1]*len(Reads)
	Counts=Reads.groupby(["umi","cbc","seq"]).sum()
	#Counts=Counts[Counts["counts"]>1]
	BCUMI=Counts.groupby(["umi","cbc"]).sum()
	BCUMI.columns=["counts_tot"]
	Counts.reset_index(inplace=True)
	BCUMI.reset_index(inplace=True)
	Counts=pd.merge(Counts,BCUMI,on=['cbc','umi'])
	Counts["TPT"]=((1.0*Counts['counts'])/(Counts['counts_tot']))
	#Counts["TPT"].to_csv("TPT.txt")
	Counts=Counts[Counts["TPT"]>.2]
	Reads=Reads_input
	Reads.reset_index(inplace=True)
	Reads=pd.merge(Reads,Counts,on=["umi","cbc","seq"])
	Reads=Reads.set_index("name")
	Reads=Reads[["umi","cbc","seq"]]
	return Reads

def error_correct(barcode,poss_bar,maxDist=2):
	#dist=[Levenshtein.hamming(barcode,bar_test) for bar_test in poss_bar]
	dist=[Levenshtein.distance(barcode,bar_test) for bar_test in poss_bar]
	if min(dist)>maxDist:
		return "";
	correct_ind=dist.index(min(dist))

	correct_bar=poss_bar[correct_ind]

	return(correct_bar)


##
##Takes in two read files, produces perturb barcode by cell barcode matrix
##realbarcodes: file containing barcode in data
def getMatrix(read1,read2,realbarcodes,realCBC,barcode="GGCACAAGCTTAATTAAGAATT",pos=33,getMid=False,getEarly=False):
	poss_bar=realbarcodes

	##loads reads
	print("Load read 1")
	R1=flattenFastq(read1)
	print("Load read 2")
	R2=flattenFastq(read2)


	pertlen=len(realbarcodes[1])
        print("Perturbation Length: "+str(pertlen))


	print("Get CBC and UMI")
	##For R1, get UMI and CBC 
	cbc=[seq[0:16] for seq in R1['seq']]
	umi=[seq[16:26] for seq in R1['seq']]

	R1["umi"]=umi
	
	R1["cbc"]=cbc

	R1=R1[["umi","cbc"]]

	
	print("Combine!")
	##combine into one data frame
	Reads=pd.concat([R1,R2],axis=1,join="inner")
	#print(Reads.shape)
	#print("Filter!")
	#Reads=filterTPT(Reads,cutoff=.2)
	if getEarly:
		return Reads
	
	print("check for primer")
	#checksubstr=[val in x for x in Reads['seq']]
	checksubstr=[barcode==x[0:len(barcode)] for x in Reads['seq']]
	print(sum(checksubstr))
	Reads=Reads[checksubstr]
	
	Reads=Reads[["AGAATT" in x for x in Reads['seq']]]
	print(Reads.shape)
	Reads=Reads[["CCTAG" in x for x in Reads['seq']]]
	
	print(Reads.shape)
	print("Filter!")
	Reads=filterTPT(Reads,cutoff=.2)

	print(Reads.shape)

	print("Get perturb barcode")

	#loc=[x.find('GGCACAAGCTTAATTAAGAATT')+len('GGCACAAGCTTAATTAAGAATT') for x in Reads['seq']]

	loc=[re.search(barcode,x).end() for x in Reads['seq']]

	#print(loc[1:10])

	Reads["location_pert"]=loc

	#Reads=Reads[[l==pos for l in loc]]
	#pos=len("TAGCAAACTGGGGCACAAGCTTAATtaagaatt")

	pert=[RevComp(Reads['seq'][i][pos:(pos+pertlen)]) for i in range(0,len(Reads['seq']))]

	Reads["pert"]=pert


	if getMid:
		return Reads;
	#print(Reads.shape)
	##error correct!
	print("Error Correct Perturbation barcodes!")

	#print(Reads.head())
	Reads["corrected_pert"]=[error_correct(barcode,poss_bar,maxDist=2) for barcode in Reads["pert"]]
	Reads=Reads[Reads["corrected_pert"]!=""]

	print(Reads.shape)
	print("Error Correct cbc")
	correctCBC=True
	poss_cbc=realCBC
	#print(Reads.head())
	if correctCBC:
		#cbc=list(set(Reads["cbc"].tolist()))
		#cbc_new=[error_correct(barcode,poss_cbc,maxDist=2) for barcode in cbc]
		

		Reads["corrected_cbc"]=[error_correct(barcode,poss_cbc,maxDist=1) for barcode in Reads["cbc"]]
	
		Reads=Reads[Reads["corrected_cbc"]!=""]
		Reads["cbc"]=Reads["corrected_cbc"]
	
	
	print(Reads.shape)
	##Take Care of duplicates!!
	print("Flatten data")
	##produces umi+cellbarcode vs perturb barcode matrix
	Reads["comb_bar"]=[x+'-'+y for x,y in zip(Reads['umi'],Reads['cbc'])]
	Reads["one"]=1
	Reads=pd.pivot_table(Reads,values="one",index="comb_bar",columns=["corrected_pert"],aggfunc=np.sum,fill_value=0)
	Reads["comb_bar"]=Reads.index
	Reads=Reads.melt(id_vars="comb_bar")
	#print(Reads.head())
	print(Reads.shape)
	#print(Reads.head())
	Reads=Reads[Reads["value"]>0]

	print("Count duplicates")
	lst=Reads["comb_bar"].tolist()
	Reads["Num"]=[lst.count(bar) for bar in Reads["comb_bar"]]

	print(Reads.shape)
	print("Remove duplicates")
	Reads=Reads[Reads["Num"]==1]

	print(Reads.shape)
	

	print("Flatten")
	##makes cellbarcode by perturb barcode count matrix
	Reads["umi"]=[a.split("-")[0] for a in Reads["comb_bar"]]
	Reads["cbc"]=[a.split("-")[1] for a in Reads["comb_bar"]]
	Reads=pd.pivot_table(Reads,values="Num",index="cbc",columns=["corrected_pert"],aggfunc=np.sum,fill_value=0)
	
	print(Reads.shape)
	print("done!")
	return(Reads)







##
##Passed arguments: location read1, location read2, location perturbation barcodes, location cell barcodes, save file name 
if __name__ == '__main__':
	args=sys.argv
	if(len(args)<7):
		print("Not enough arguments!")
		exit()
	print(args[1])
	print(args[2])
	print(args[3])
	print(args[4])
	print(args[5])
	print(args[6])
	fil=open(args[3]);bar=fil.readlines();fil.close();bar=[b.strip() for b in bar]
	fil=open(args[4]);cbc=fil.readlines();fil.close();cbc=[b.strip() for b in cbc]
	cbc=[cb[0:16] for cb in cbc]
	primer=args[6]
	
	pos=33
	if primer=="GCAAACTG":
		pos=31
	print(pos)
	mat=getMatrix(args[1],args[2],bar,cbc,pos=pos,barcode=primer)
	mat.to_csv(args[5])
