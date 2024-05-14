#!/usr/bin/env   python


#  junctMhomol.py is a tidied-up version of _junct_mhomol_v5.py C. White 2020-2023

# starting with a list of deletions reformatted from the variant caller output, this script extracts and formats junctions of deletions, 
##    scores presence of microhomology potentially implicated in the deletions and writes output to files.
##    The microhomology scoring counts down from a chosen length (usually 9nt) and scores the longest microhomol seq it finds potentially involved at the junction


### theChromosome must be in the format of the refseq fasta file - eg "Chr1"
### theLength = the length of the sequence blocks to extract and to consider for microhomologies 
### note. if the deletion is <= theLength, the whole deletion seq is used and the sequence blocks are limited to the deletion length.

## output file has "_delsJunctions.txt" suffix
## the last two columns are the numbers of occurrences of this deletion found by the variant caller (r1p and r1m) 
##    an example output line : AAAGCTGATG ( ACAGA ) CACACATCTT  del =   -5  microhomol =    1   Chr1    5189174 5189179 0   1

# call as junctMhomol.py refseq.fa delsList.txt outfile.txt theChromosome theLength
 

import sys, os, re
from Bio import SeqIO

if  (len(sys.argv) != 6) :
    print("starting with a list of deletions reformatted from the variant caller output, this script extracts and formats junctions of deletions,") 
    print("scores presence of microhomology at junctions and writes output to files.")
    print(" call as junctMhomol.py genomerefSeqFile.fa listFile.txt outFile.txt theChromosome windowLength")
    print("     note. if the deletion is <= theLength, the whole deletion seq is used and the sequence blocks analysed are limited to the deletion length.")
    print()
    print("cw 28/11/2022")

    sys.exit()

# print("starting with a list of deletions reformatted from the variant caller output, this script extracts and formats junctions of deletions,") 
# print("scores presence of microhomology at junctions and writes output to files.")
# print("     note. if the deletion is <= theLength, the whole deletion seq is used and the sequence blocks analysed are limited to the deletion length.")
# print("C. White 2020-2023")
# print()


refSeqFile = open(sys.argv[1],'r')
listFile = open(sys.argv[2],'r')
outFile = open(sys.argv[3],'a')
theChr = sys.argv[4]
theLength = int(sys.argv[5])

count = 0
mh = [0]*10 
mhLabel = []

for i in range(0, theLength):
    mhLabel+= ["mh"+str(i)]

del5_start = del5_end = upstrStart = upstrEnd = del3_end = del3_start = 0

theLine = []
lines = listFile.readlines()

theSample = re.split("[/.]", sys.argv[2])
theSampleName = theSample[len(theSample) - 2]

print(" junction", "txt1", "del", "txt2", "mhomol", "Chr", "start", "end" , "r1p", "r1m", sep="\t", file=outFile )

for record in SeqIO.parse( refSeqFile.name , "fasta"):
        if (record.id == theChr) :
            for theLine in lines :
                info = theLine.split()
                if (info and (info[0][0] != "#")) :
                    theChr = info[0]
                    theDelLen = int(info[3])
                    delSeq = del5_seq = del3_seq = ""
                    counter = 0
                    mhomolFound = ""

                    r1p = info[4]
                    r1m = info[5]

                    del5_start = int(info[1]) 
                    del5_end = del5_start + theLength
                    upstrStart = del5_start - theLength
                    del3_end = int(info[2]) 
                    del3_start = del3_end - theLength

                    upstrSeq = record.seq[upstrStart:del5_start]
                    downstrSeq = record.seq[del3_end:del3_end + theLength]

                    if (abs(theDelLen) <= theLength) :
                        delSeq = record.seq[del5_start:del3_end]
                        print( upstrSeq + " ( " + delSeq  + " ) " +  downstrSeq, "del = ", del5_start - del3_end , sep="\t", end="\t", file=outFile )
                        count += 1
                        counter = len(delSeq)
                        while ((counter > 0) and  (mhomolFound == "")):
                            if ((upstrSeq[0:counter] == delSeq[0:counter]) or (delSeq[0:counter] == downstrSeq[0:counter]) ):
                                print("microhomol = ", counter, theChr, del5_start, del3_end, r1p, r1m, sep="\t", file=outFile )
                                mh[counter] +=1 
                                mhomolFound = "yes"
                            counter -=1

                        if mhomolFound == "" :
                                mh[0] +=1
                                print("microhomol = ", "0" , theChr, del5_start, del3_end, r1p, r1m, sep="\t", file=outFile )
                    else :
                        del5_seq = record.seq[del5_start:del5_end]
                        del3_seq = record.seq[del3_start:del3_end]
                        if (abs(theDelLen) <= 2*theLength):
                            delSeq = record.seq[del5_start:del3_end]
                            print( upstrSeq + " ( " + delSeq  + " ) " +  downstrSeq, "del = ", del5_start - del3_end ,  sep="\t", end="\t", file=outFile )
                        else :
                            print( upstrSeq + " ( " + del5_seq + ".." + del3_seq + " ) " +  downstrSeq, "del = ", del5_start - del3_end ,  sep="\t", end="\t", file=outFile )
                        count += 1
                        counter = theLength - 1
                        while ((counter > 0) and  (mhomolFound == "")):
                            if ((upstrSeq[0:counter] == del3_seq[0:counter]) or (del5_seq[0:counter] == downstrSeq[0:counter])) :
                                print("microhomol = ", counter, theChr, del5_start, del3_end, r1p, r1m, sep="\t", file=outFile )
                                mh[counter] +=1
                                mhomolFound = "yes"
                            counter -=1

                        if mhomolFound == "" :
                                mh[0] +=1
                                print("microhomol = ", counter, theChr, del5_start, del3_end, r1p, r1m, sep="\t", file=outFile )

print("# ", " name", "total", *mhLabel,  sep="\t", file=outFile ) 
print("#@", theSampleName, count, *mh, sep="\t", file=outFile ) 

print()
print ( str(count) + " seqs found" )

theFile = outFile.name
outFile.close()
os.system("sort " + theFile + " -o " + theFile)


