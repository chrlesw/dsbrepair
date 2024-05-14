#!/usr/bin/env   python

import sys, os, string, datetime

#############################
###varScorer.py is a tidied-up version of _variantScorer_v6bis.py . C .White 13/03/2022  

## call as varScorer.py Chr cs_1 [cs_2] halfWindow outNameRoot
## note : halfWindow is half of the window around the cut (ie if = 5 will check a 10 base window centered on the cut)

## if run standalone, it needs to be run from the _vars folder produced by analDSB.py

## note that the Chromosome name comes from the input parameter and might be incorrect if any sequences map elsewhere
#############################

if  (len(sys.argv) == 5)  :
	theChr = sys.argv[1]
	t1_cs = sys.argv[2]
	offset = sys.argv[3]
	fileNameRoot = sys.argv[4]
	numSites = 1

elif  (len(sys.argv) == 6)  :
	theChr = sys.argv[1]
	t1_cs = sys.argv[2]
	t2_cs = sys.argv[3]
	offset = sys.argv[4]
	fileNameRoot = sys.argv[5]
	numSites = 2

else :
	print()
	print("varScorer.py " )
	print(" run from _vars folder produced by  analDSB.py")
	print("inputs :  Chromosome(eg Chr1) t1_cs [t2_cs] offset fileNameRoot   " )
	print("for FnCas12a : cs = end of protospacer - 3nt ")
	print()
	print("varScorer.py, C.White 28/11/2022 ")
	print()
	sys.exit()

##############################################

os.system("echo \"type\tchr\tstart\tend\tlen\tr1p\tr1m\tins_seq\" > " + fileNameRoot + "_indelsList.txt")

os.system("awk  ' {FS=\"\t\"}  $4==\"INS\" {print $4, \"" + theChr +  "\" , $2 , $3, length($5), $6, $7, $5 }' OFS='\t' " +  fileNameRoot + "_vars.txt   | sort --key=2 -h >> " + fileNameRoot + "_indelsList.txt")
os.system("awk  ' {FS=\"\t\"} $4==\"SUB\" {print $4, \"" + theChr +  "\", $2 , $3, length($5), $6, $7, $5  }' OFS='\t' " +  fileNameRoot + "_vars.txt   | sort --key=2 -h >> " + fileNameRoot + "_indelsList.txt")
os.system("awk  ' {FS=\"\t\"} $4==\"DEL\" {print $4, \"" + theChr +  "\", $2 , $3, $3 - $2, $6, $7  }' OFS='\t' " +  fileNameRoot + "_vars.txt   | sort --key=2 -h >> " + fileNameRoot + "_indelsList.txt")

os.system("awk  ' {FS=\"\t\"} {print \"DEL\", $6, $7 , $8, $8 - $7, $9, $10, $5  }' OFS='\t' ../_junctions/" + fileNameRoot + "_delsJunctions.txt   | sort --key=2 -h >> " + fileNameRoot + "_mHomolDelsList.txt")

##############################################

if (numSites == 1) :

	os.system("echo \"type\tchr\tstart\tend\tt1_del\tr1p\tr1m\tmhomol\"  > " + fileNameRoot + "_all_dels.txt")
	os.system("echo \"type\tchr\tstart\tend\tt1_ins\tr1p\tr1m\tt1_ins_seq\"  > " + fileNameRoot + "_all_ins.txt")
	os.system("echo \"type\tchr\tstart\tend\tt1_sub\tr1p\tr1m\tt1_sub_seq\"  > " + fileNameRoot + "_all_sub.txt")

	#####pull out the t1 modifs #####
#	os.system("cat " + fileNameRoot + "_indelsList.txt | awk '( (NR>1) && ($3 <= " + t1_cs + " + " + offset + " )  && ($4 >= " + t1_cs + " - " + offset + " ) ) {print}' | awk '($1==\"DEL\") {print}' OFS='\t' >> " + fileNameRoot + "_t1_dels.txt ")
	os.system("cat " + fileNameRoot + "_mHomolDelsList.txt | awk '( (NR>1) && ($3 <= " + t1_cs + " + " + offset + " )  && ($4 >= " + t1_cs + " - " + offset + " ) ) {print}' | awk '($1==\"DEL\") {print}' OFS='\t' >> " + fileNameRoot + "_all_dels.txt ")
	os.system("cat " + fileNameRoot + "_indelsList.txt | awk '( (NR>1) && ($3 <= " + t1_cs + " + " + offset + " )  && ($4 >= " + t1_cs + " - " + offset + " ) ) {print}' | awk '($1==\"INS\") {print}' OFS='\t' >> " + fileNameRoot + "_all_ins.txt " )
	os.system("cat " + fileNameRoot + "_indelsList.txt | awk '( (NR>1) && ($3 <= " + t1_cs + " + " + offset + " )  && ($4 >= " + t1_cs + " - " + offset + " ) ) {print}' | awk '($1==\"SUB\") {print}' OFS='\t' >> " + fileNameRoot + "_all_sub.txt")

elif (numSites == 2) :
	os.system("echo \"type\tchr\tstart\tend\tall_dels\tr1p\tr1m\tmhomol\"  > " + fileNameRoot + "_all_dels.txt")
	os.system("echo \"type\tchr\tstart\tend\tall_ins\tr1p\tr1m\tall_ins_seq\"  > " + fileNameRoot + "_all_ins.txt")
	os.system("echo \"type\tchr\tstart\tend\tall_sub\tr1p\tr1m\tall_sub_seq\"  > " + fileNameRoot + "_all_sub.txt")

	os.system("echo \"type\tchr\tstart\tend\tt1_del\tr1p\tr1m\tmhomol\"  > " + fileNameRoot + "_t1_dels.txt")
	os.system("echo \"type\tchr\tstart\tend\tt1_ins\tr1p\tr1m\tt1_ins_seq\"  > " + fileNameRoot + "_t1_ins.txt")
	os.system("echo \"type\tchr\tstart\tend\tt1_sub\tr1p\tr1m\tt1_sub_seq\"  > " + fileNameRoot + "_t1_sub.txt")

	os.system("echo \"type\tchr\tstart\tend\tt2_del\tr1p\tr1m\tmhomol\"  > " + fileNameRoot + "_t2_dels.txt")
	os.system("echo \"type\tchr\tstart\tend\tt2_ins\tr1p\tr1m\tt2_ins_seq\"  > " + fileNameRoot + "_t2_ins.txt")
	os.system("echo \"type\tchr\tstart\tend\tt2_sub\tr1p\tr1m\tt2_sub_seq\"  > " + fileNameRoot + "_t2_sub.txt")

##### pull out all deletions (t1 and t2) #####
#	os.system("cat " + fileNameRoot + "_indelsList.txt | awk '( (NR>1) && ( ($3 <= " + t1_cs + " + " + offset + " )  && ($4 >= " + t1_cs + " - " + offset + " ) ) || (( $3 <= " + t2_cs + " + " + offset + " )  && ($4 >= " + t2_cs + " - " + offset + " ) )) {print}' | awk '($1==\"DEL\") {print}' OFS='\t' >> " + fileNameRoot + "_all_dels.txt ")
	os.system("cat " + fileNameRoot + "_mHomolDelsList.txt | awk '( (NR>1) && ( ($3 <= " + t1_cs + " + " + offset + " )  && ($4 >= " + t1_cs + " - " + offset + " ) ) || (( $3 <= " + t2_cs + " + " + offset + " )  && ($4 >= " + t2_cs + " - " + offset + " ) )) {print}' | awk '($1==\"DEL\") {print}' OFS='\t' >> " + fileNameRoot + "_all_dels.txt ")
	os.system("cat " + fileNameRoot + "_indelsList.txt | awk '( (NR>1) && ( ($3 <= " + t1_cs + " + " + offset + " )  && ($4 >= " + t1_cs + " - " + offset + " ) ) || (( $3 <= " + t2_cs + " + " + offset + " )  && ($4 >= " + t2_cs + " - " + offset + " ) )) {print}' | awk '($1==\"INS\") {print}' OFS='\t' >> " + fileNameRoot + "_all_ins.txt ")
	os.system("cat " + fileNameRoot + "_indelsList.txt | awk '( (NR>1) && ( ($3 <= " + t1_cs + " + " + offset + " )  && ($4 >= " + t1_cs + " - " + offset + " ) ) || (( $3 <= " + t2_cs + " + " + offset + " )  && ($4 >= " + t2_cs + " - " + offset + " ) )) {print}' | awk '($1==\"SUB\") {print}' OFS='\t' >> " + fileNameRoot + "_all_sub.txt ")
	#####pull out the t1 modifs#####
#	os.system("cat " + fileNameRoot + "_indelsList.txt | awk '( (NR>1) && ($3 <= " + t1_cs + " + " + offset + " )  && ($4 >= " + t1_cs + " - " + offset + " ) ) {print}' | awk '($1==\"DEL\") {print}' OFS='\t' >> " + fileNameRoot + "_t1_dels.txt ")
	os.system("cat " + fileNameRoot + "_mHomolDelsList.txt | awk '( (NR>1) && ($3 <= " + t1_cs + " + " + offset + " )  && ($4 >= " + t1_cs + " - " + offset + " ) ) {print}' | awk '($1==\"DEL\") {print}' OFS='\t' >> " + fileNameRoot + "_t1_dels.txt ")
	os.system("cat " + fileNameRoot + "_indelsList.txt | awk '( (NR>1) && ($3 <= " + t1_cs + " + " + offset + " )  && ($4 >= " + t1_cs + " - " + offset + " ) ) {print}' | awk '($1==\"INS\") {print}' OFS='\t' >> " + fileNameRoot + "_t1_ins.txt " )
	os.system("cat " + fileNameRoot + "_indelsList.txt | awk '( (NR>1) && ($3 <= " + t1_cs + " + " + offset + " )  && ($4 >= " + t1_cs + " - " + offset + " ) ) {print}' | awk '($1==\"SUB\") {print}' OFS='\t' >> " + fileNameRoot + "_t1_sub.txt")

	#####pull out the t2 modifs#####
#	os.system("cat " + fileNameRoot + "_indelsList.txt | awk '( (NR>1) && ($3 <= " + t2_cs + " + " + offset + " )  && ($4 >= " + t2_cs + " - " + offset + " ) ) {print}' | awk '($1==\"DEL\") {print}' OFS='\t' >> " + fileNameRoot + "_t2_dels.txt ") 
	os.system("cat " + fileNameRoot + "_mHomolDelsList.txt | awk '( (NR>1) && ($3 <= " + t2_cs + " + " + offset + " )  && ($4 >= " + t2_cs + " - " + offset + " ) ) {print}' | awk '($1==\"DEL\") {print}' OFS='\t' >> " + fileNameRoot + "_t2_dels.txt ") 
	os.system("cat " + fileNameRoot + "_indelsList.txt | awk '( (NR>1) && ($3 <= " + t2_cs + " + " + offset + " )  && ($4 >= " + t2_cs + " - " + offset + " ) ) {print}' | awk '($1==\"INS\") {print}' OFS='\t' >> " + fileNameRoot + "_t2_ins.txt " )
	os.system("cat " + fileNameRoot + "_indelsList.txt | awk '( (NR>1) && ($3 <= " + t2_cs + " + " + offset + " )  && ($4 >= " + t2_cs + " - " + offset + " ) ) {print}' | awk '($1==\"SUB\") {print}' OFS='\t' >> " + fileNameRoot + "_t2_sub.txt")




