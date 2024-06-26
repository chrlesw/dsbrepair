#!/usr/bin/env   python


######################## _dsb_varAnal.py #########################################
###### C. White 2020-2024
###### _dsb_varAnal.py is a tidied-up version of _alignVarArc_v3.py
# 
# It should be run from the directory containing the .fastq or fastq.gz files generated by generated by the nanopore sequencer or from merging paired-end Illumina-type sequencing output.
# It steps through the files in the current directory and for each fastq file which includes the fileNameRoot string in its name : 
#      - for MiSeq data, the prefix and suffix seqs are appended to each sequence (if no prefix and suffix seqs are given, empty strings (nothing) will be appended). 
#            This anchors the alignments of sequences with deletions and is recommended for short (Illumina) sequences.
#      - the seqs are aligned to the reference genome with minimap2, using Nanopore optimisations for nanopore seqs
#      -  secondary and supplementary alignments are removed
#      -  variants are called and lists of selected variants which affect a window centered on the cut site(s) built
#      -  microhomologies potentially implicated in the deletions are quantified
#      -  the deletion lists are expanded to have one deletion per line (ie X occurrences of a given deletion will result in X lines)
#            this expansion produces two text files - one with all deletions and one with the number of deletions 
#            defined by the "numLines" parameter and deletion lengths >= to the "minDel" parameter (selected from full list after randomisation of the lines)
#      - these deletion files are then used to draw arc diagram plots with the delsArcs.R  script.

#       note that inversions will need to be aligned to a different reference seq
# ####################################################################################################
# the script calls (and needs) the following scripts which need to be in the $PATH :
#      - varScorer.py  
#      - varBuildFreqs.py  
#      - varexpandList.py 
#      - junctMhomol.py
#      - delsArcs.R  (the full path to this is passed as an argument)
#      - mhomolMeansVdels.py 
####################################################################################################
#### note - the script will deal properly with calls with one or two cut-sites #####
####      - the cutsiteCenter_2  parameter should be left out if only one cut-site is used #####
#### note - -preSeq and -postSeq are appended at the beginning and end (respectively) of sequences to anchor mapping to the reference
####           (this is v)ery useful for small amplicons (eg. Illumina) with large deletions) 
####      - for nanopore sequence analyses, the prefix and suffix seqs are ignored
####      - they can be left out if not needed (need both or neither)

#### note - the -targ parameter is passed to the delsArcs.R script for selection of the Arc map formatting setup (see the delsArcs.R script for details) 


#  -npore or --isNanoporeSeq : required=True (are these nanopore sequences (yes or no))
#  -nameSel or --fileNameSelector :	required=True	(text in filename for choice of files to be analysed for alignments and variant analyses)
#  -ref or --refFasta :	required=True	(reference sequence or genome (fasta file))
#  -chr or --theChr :	required=True	(chromosome name of target site (eg Chr1) )
#  -cs1 or --cutsiteCenter_1 :	required=True	(coordinate of the center of cut-site 1)
#  -cs2 or --cutsiteCenter_2 :	required=False	(coordinate of the center of cut-site 2)
#  -win or --halfWindow :	required=True	(the window centered on the cut-site is +/- halfWindow)
#  -nl or --numLines, type=int :	required=True	(number of deletions to extract for arc diagrams)
#  -md or --minDel, type=int :	required=True	(minimum deletion length for arc diagrams)
#  -inclFocus or --includeArcFocus : required=True, choices=["yes", "no"], (also generate arc diagrams focussed locally on the cut-site(s))
#  -plotMhomol or --includeMicrohomolArcs : required=True, default="yes", choices=["yes", "2colours", "no"], help="colour arc diagrams with microhomology data, yes/2colours/no (2colours gives µhomol ≤ 2 = blue, ≥3 = red;  no = b&w arcs)")
#  -ampLen or --ampliconLength : required=False, default="0", help="integer length of the amplicon")
#  -targ or --targetLocus : required=True, choices=["site2_915" , "chr3_92" , "ercc1_2ab" , "ercc1_2cd" , "xpf_2ab" , "send1_2cd" , "send1_2ab" , "iugus"], 
#             help="the target locus (see list)"). 
#  -rp or --rscriptPath :	required=True	(path to the delsArcs.R script)
#  -ampLen or --ampliconLength : required=False
#  -preSeq or --prefix_seq:	required=False	(upstream seq for appending to sequences  ())
#  -postSeq or --postfix_seq :	required=False	(downstream seq for appending to sequences)

####################################################################################################
############## for example (nanopore target site2_915 ), call as:   #############################
# _dsb_varAnal.py -npore yes -targ site2_915 -nameSel aString -chr Chr1 -cs1 5187987 -cs2 5188243 -win 5 -nl 1000 -md 1 -inclFocus no \
#  -plotMhomol no  -ampLen 2800 -ref /Users/self/_informatics/seqs_ref/tair10/TAIR10_all.fa \
#  -rp /Users/self/.bioinf/ngsscripts/cas_anal/delsArcs.R 

############## for example (miseq target site2_915), call as: #############################
# _dsb_varAnal.py -npore no -targ site2_915 -nameSel aString -chr Chr1 -cs1 5187987 -cs2 5188243 -win 5 -nl 1000 -md 1 -inclFocus no \
# -plotMhomol no  -ampLen 400 -ref /Users/self/_informatics/seqs_ref/tair10/TAIR10_all.fa \
# -rp /Users/self/.bioinf/ngsscripts/cas_anal/delsArcs.R \
#  -preSeq GAGATTTTCCGTTTCCCGCTTTAATACAGTGCCCCAATTCGCGCGACACATAGAGTGTAGAGACGCTTTCACGAGCGTTTCCGACGTCGGACTTTCAGCTCATCATCTCCACATCTTTAACGGTAAAGGTACACTCTCATCGTCTTCTCTCCTATTGATCACTGATGAGTATTGAGACACAGGTTATAGTGGAATTTAGT \
# -postSeq GTAATACATTTTGGATTTTCCTGATTTCTTTTATGCAACATGTAAATCTTGGACGTGGTCTTTACAAGTGTTTTTTTTTTTATGACAGTTGATAGCTGTTGTGCTACCTTTTGCTGTCATTTGTGTTTACTACTTCATTAGAAATGATGTTTATGACCTGCATCATGCAATACTAGGTACCTCTTTCTTCTTCCAATTGT

### preseq and postseq are 200nt upstream and downstream sequences flanking the amplicon that help to anchor the alignments (can be of user-determined length)
####################################################################################################
####################################################################################################

import sys, os, string, datetime, argparse, gzip
from Bio import SeqIO 
from Bio.SeqRecord import SeqRecord

print()
print("######################## _dsb_varAnal.py #########################################")
print("                                                                C. White 2020-2024")
print()
print("##### This script should be run from the directory containing the fastq or fastq.gz files")
print("#####     generated by the Nanopore sequencer or by merging paired-end Illumina-type sequencing ouput")
print()
###########################################################################################################
totalSeqs = totalFiles = 0
theChr = numsites = ""
###########################################################################################################
parser = argparse.ArgumentParser()

parser.add_argument("-npore", "--isNanoporeSeq", required=True, choices=["yes", "no"], help="is this nanopore sequences (yes or no)")
parser.add_argument("-targ", "--targetLocus", required=True, choices=["site2_915" , "chr3_92" , "ercc1_2ab" , "ercc1_2cd" , "xpf_2ab" , "send1_2cd" , "send1_2ab" , "iugus" , "chr2_184"], help="the target locus (see list)")
parser.add_argument("-nameSel", "--fileNameSelector", required=True, help="text in filename for choice of files for alignments and variant analyses")
parser.add_argument("-ref", "--refFasta", required=True, help="reference sequence or genome (fasta file)")
parser.add_argument("-chr", "--theChr", required=True, help="chromosome name of target site (eg Chr1) ")
parser.add_argument("-cs1", "--cutsiteCenter_1", required=True, help="coordinate of the center of cut-site 1")
parser.add_argument("-cs2", "--cutsiteCenter_2", required=False, help="coordinate of the center of cut-site 2")
parser.add_argument("-win", "--halfWindow", required=True, help="the window centered on the cut-site is +/- halfWindow")
parser.add_argument("-nl", "--numLines", type=int, required=True, help="number of deletions to extract for arc diagrams")
parser.add_argument("-md", "--minDel", type=int, required=True, help="minimum deletion length for inclusion in arc diagrams")
parser.add_argument("-inclFocus", "--includeArcFocus", required=True, choices=["yes", "no"], help="also generate arc diagrams focussed locally on the cut-site(s)")
parser.add_argument("-plotMhomol", "--includeMicrohomolArcs", required=True, default="yes", choices=["yes", "2colours", "no"], help="colour arc diagrams with microhomology data, yes/2colours/no (2colours gives µhomol ≤ 2 = blue, ≥3 = red;  no = b&w arcs)")
parser.add_argument("-ampLen", "--ampliconLength", required=False, default="0", help="length of the amplicon")

parser.add_argument("-rp", "--rscriptPath", required=True, help="path to the delsArcs.R script")

parser.add_argument("-preSeq", "--prefix_seq", required=False, default='None', help="upstream seq for appending to sequences")
parser.add_argument("-postSeq", "--postfix_seq", required=False, default='None', help="downstream seq for appending to sequences")

args = parser.parse_args()

nanoporeSeq = args.isNanoporeSeq
targ = args.targetLocus
theFileRoot = args.fileNameSelector
theRefSeq = args.refFasta
theChr = args.theChr
cutsiteCenter_1 = args.cutsiteCenter_1

halfWindow = args.halfWindow
numLines = args.numLines
minDel = args.minDel

rscriptPath = args.rscriptPath
includeFocussedArcDiags = args.includeArcFocus
plotMhomol = args.includeMicrohomolArcs
ampLen = args.ampliconLength



if (str(args.cutsiteCenter_2) == "None") :
	cutsiteCenter_2 = 0
	numsites = "1"
	print("One cut-site centered at ", theChr + ":" +  str(cutsiteCenter_1) )
else:
	cutsiteCenter_2 = args.cutsiteCenter_2
	numsites = "2"
	print("Two cut-sites centered at ", theChr + ":" +  str(cutsiteCenter_1) + " , " + theChr + ":" + str(cutsiteCenter_2))

if ((str(args.prefix_seq) == "None") or (str(args.postfix_seq) == "None") ) :
	prefixSeq = suffixSeq = ""
	print("no prefix or postfix seqs appended")
else:
	prefixSeq = args.prefix_seq
	suffixSeq = args.postfix_seq
	print("prefix and postfix seqs appended")

################### set things up ###################################################

seqs_read = 0
top = os.getcwd()

today = datetime.date.today()
theDate = str(f'{today:%Y_%m_%d}')

if (nanoporeSeq == "no"):
	theAlignsFolder = "_alignments_variants_"  + theDate + "/"
elif (nanoporeSeq == "yes"):
	theAlignsFolder = "_alignments_variants_ont_"  + theDate + "/"

try:
	os.makedirs(theAlignsFolder)
except OSError:
	print("the alignments folder already exists")
	sys.exit()

try:
	os.makedirs(theAlignsFolder + "_bam_files")
except OSError:
	print("the bam folder already exists")
	sys.exit()

try:
	os.makedirs(theAlignsFolder + "_vcf_files")
except OSError:
	print("the vcf folder already exists")
	sys.exit()

try:
	os.makedirs(theAlignsFolder + "_vars")
except OSError:
	print("the vars folder already exists")
	sys.exit()

try:
	os.makedirs(theAlignsFolder + "_aln_files")
except OSError:
	print("the _aln_files folder already exists")
	sys.exit()

try:
	os.makedirs(theAlignsFolder + "_junctions")
except OSError:
	print("the _junctions folder already exists")
	sys.exit()

try:
	os.makedirs(theAlignsFolder + "_arcs")
except OSError:
	print("the _arcs folder already exists")
	sys.exit()

######################################################################################
######################################################################################

for root, dirs, files in os.walk(top) :
	for fname in files:
		if (theFileRoot in fname) :
			totalFiles +=1
			outFileName = ""
			seqs_read = 0
			fnameRoot = fname.split(".")[0]

			outVarsFileRoot = theAlignsFolder + fname.split(".")[0]

			print()
			print("######################################################################")
			print("mapping sequences and calling variants  :  " , "file = ", fname.split(".")[0]) 
			print()

######################################################################################
###########  analyses of merged paired-end Illumina sequences ########################

			if (nanoporeSeq == "no"):
				################# for short amplicons, append seq blocks before and after the seqs, write to fasta file #####

				outFileName = theAlignsFolder + fname.split(".")[0] + "_appended.fasta"
				f_out = open(outFileName, 'a')

				######################################################################
				###########  first append prefix and suffix sequences ################

				if ("fastq.gz" in fname):
					for seq_record in SeqIO.parse(gzip.open(fname, mode='rt'), 'fastq'):
						seqs_read +=1 
						totalSeqs +=1
						seq_record = prefixSeq + seq_record
						seq_record += suffixSeq
						r=SeqIO.write(seq_record, f_out, 'fasta')
						if r!=1: print('Error while writing sequence:  ' + seq_record.id)
				else :
					for seq_record in SeqIO.parse(open(fname, mode='r'), 'fastq'):
						seqs_read +=1 
						totalSeqs +=1
						seq_record = prefixSeq + seq_record
						seq_record += suffixSeq
						r=SeqIO.write(seq_record, f_out, 'fasta')
						if r!=1: print('Error while writing sequence:  ' + seq_record.id)
				
				f_out.close() 

				######################################################################
				###########  map sequences to reference and call  ################

				os.system("minimap2 --cs=long -a " + theRefSeq + " " + outFileName + " > " + outFileName + "temp.sam")
			      ### use samtools view to remove secondary and supplementary alignments
				os.system("samtools view -h -F 0x900  " + outFileName + "temp.sam > " + outFileName + ".sam")
				os.system("paftools.js sam2paf " + outFileName + ".sam > " + outFileName + ".paf")
				os.system("samtools sort  " + outFileName + ".sam -o " + outFileName + ".sorted.bam && samtools index  " + outFileName + ".sorted.bam")			
				os.system("callvariants.sh in=" + outFileName + ".sorted.bam  out=" + outVarsFileRoot + "_vars.txt 32bit=t vcf=" + outVarsFileRoot + ".vcf.gz  ref=" + theRefSeq + " ploidy=2 clearfilters bgzip=t extended=f ")

######################################################################################
###########  analyses of Nanopore sequences ##########################################

			elif (nanoporeSeq == "yes"):
				outFileName = theAlignsFolder + fname
				outFileRoot = theAlignsFolder + fnameRoot
				os.system("minimap2 --cs=long -ax map-ont " + theRefSeq + " " + fname + " > " + outFileRoot + "temp.sam")
					### use samtools view to remove secondary and supplementary alignments
				os.system("samtools view -h -F 0x900  " + outFileRoot + "temp.sam > " + outFileRoot + ".sam")
				os.system("paftools.js sam2paf " + outFileRoot + ".sam > " + outFileRoot + ".paf")
				os.system("samtools sort  " + outFileRoot + ".sam -o " + outFileRoot + ".sorted.bam && samtools index  " + outFileRoot + ".sorted.bam")
				os.system("callvariants.sh in=" + outFileRoot + ".sorted.bam  out=" + outFileRoot + "_vars.txt 32bit=t vcf=" + outFileRoot + ".vcf.gz  ref=" + theRefSeq + " ploidy=2 clearfilters bgzip=t extended=f ")

##############################################################################
#### generate file with a formatted list of deletions, sorted on deletion length ##########	

			with open(outVarsFileRoot + "_vars.txt") as input:
				fileVarsListFile = open(outVarsFileRoot + "_delsList.txt", 'a')
				for line in input:
					if(line[0] != "#"):
						fields = line.split('\t')
						if (fields[3] == "DEL"):
							print( theChr, fields[1], fields[2], int(fields[1]) - int(fields[2]), fields[5], fields[6], sep='\t', file=fileVarsListFile)
				fileVarsListFile.close()
				
				# sort the output file lines by the deletion length (4th column)
				os.system("sort --key=4 -h -o " + outVarsFileRoot + "_delsList.txt " + outVarsFileRoot + "_delsList.txt") 
				# add column titles line to output file
				os.system("echo \"#type\tchr\tstart\tend\tlength_del\tr1p\tr1m\"  >> " + outVarsFileRoot + "_delsList.txt")
			
			### do microhomology analyses ######		
			os.system("junctMhomol.py " + theRefSeq + " " + outVarsFileRoot + "_delsList.txt " + outVarsFileRoot + "_delsJunctions.txt " + theChr + " 10" )
			os.system("tabix -p vcf "+ outVarsFileRoot + ".vcf.gz")

	break  ##to avoid recursive file walk (will not analyse files in sub-directories)

#########################################################################
############## build microhomology stats file ###########################

mhomolFile = open( theAlignsFolder + "_stats_mhomol_Junctions.txt", 'a')
print("##use of microhomology at deletion junctions. mh1=1nt, mh2=2nt...", file=mhomolFile )
print("#name", "total", "mh0", "mh1", "mh2", "mh3", "mh4", "mh5", "mh6", "mh7", "mh8", "mh9", \
	"fr_mh0", "fr_mh1", "fr_mh2", "fr_mh3", "fr_mh4", "fr_mh5", "fr_mh6", "fr_mh7", "fr_mh8", "fr_mh9", sep="\t", file=mhomolFile )
mhomolFile.close()

os.system("cat " + theAlignsFolder + "*Junctions.txt | awk  '{if(($1==\"#@\") && ($3>0)) {print $2 ,$3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $4/$3, $5/$3, $6/$3, $7/$3, $8/$3, $9/$3, $10/$3, $11/$3, $12/$3, $13/$3} else if ($1==\"#@\") {print $2 ,$3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13} } ' OFS='\t' >> " + theAlignsFolder + "_stats_mhomol_Junctions.txt" )

##############################################################################
###########################     clean things up    ###########################
##############################################################################
print()
print("cleaning things up")
print()

os.system("mv " + theAlignsFolder + "*.vcf.gz  " + theAlignsFolder + "_vcf_files/. ")
os.system("mv " + theAlignsFolder + "*.vcf.gz.tbi " + theAlignsFolder + "_vcf_files/. ")
os.system("mv " + theAlignsFolder + "*Junctions.txt " + theAlignsFolder + "_junctions/. ")
os.system("mv " + theAlignsFolder + "*delsList.txt " + theAlignsFolder + "_junctions/. ")
os.system("mv " + theAlignsFolder + "*.paf " + theAlignsFolder + "_aln_files/. ")
os.system("mv " + theAlignsFolder + "*.bam*  " + theAlignsFolder + "_bam_files/.")
os.system("mv " + theAlignsFolder + "*.txt " + theAlignsFolder + "_vars/. ")

os.system("rm " + theAlignsFolder + "*.sam ")
if (nanoporeSeq == "no") :
	os.system("rm " + theAlignsFolder + "*.fasta ")

############################################################################
os.chdir( theAlignsFolder + "_vars")

for filename in os.listdir( '.' ) :
	thelist = filename.split("_")
	if thelist[len(thelist)-1] == "vars.txt":
		theRoot = filename.split("_vars.txt")[0]
		if cutsiteCenter_2 == 0 :
			os.system("varScorer.py " + theChr + " " + str(cutsiteCenter_1) + " " + str(halfWindow) + " " + theRoot)
		else : 
			os.system("varScorer.py " + theChr + " " + str(cutsiteCenter_1) + " " + str(cutsiteCenter_2) + " " + str(halfWindow) + " " + theRoot)

os.system("varBuildFreqs.py _dels")
os.system("varBuildFreqs.py _ins")
os.system("varBuildFreqs.py _sub")
os.system("varBuildFreqs.py _all")

os.system("varexpandList.py all_dels " + str(numLines) + " " + str(minDel) )

os.chdir("../..")
arcsPath = os.getcwd() + "/" + theAlignsFolder + "_arcs"
os.system("mv  " + theAlignsFolder + "_vars/_exp*.txt " + arcsPath )

if(ampLen != 0) :
	os.chdir(theAlignsFolder + "/_arcs")
	os.system("mhomolMeansVdels.py  _exp " + ampLen )
	os.chdir("../..")

# call the R script to draw the arc diagram plots
os.system("Rscript " + rscriptPath + " " + arcsPath + " " + nanoporeSeq + " " + targ + " " + includeFocussedArcDiags + " " + numsites + " " + plotMhomol)

##############################################################################
print("*******************************************************")
print(totalFiles, " files processed")
print("*******************************************************")
##############################################################################
