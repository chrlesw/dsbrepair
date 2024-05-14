#!/usr/bin/env   python

import sys, os, string, datetime
import pandas as pd

#C. White 2022-2023

# varBuildFreqs.py is a cleaned-up version of _varsbuild_freqs_v6.py . C. White 20/07/2022

# _varsbuild_freqs.py generates a table of numbers of events from callvariants.sh output
# builds a dictionary of the lengths from column 4 (0 based) and sums numbers of identical events from columns 5 and 6
# it scans all files in the current folder and processes only those with names containing the fileNameRoot string

# output files are named  "_freq" + fileNameRoot

# note: the keyList range must be longer than the sequence length
# note: the "r1p" and "r1m" fields of the callvariants.sh output are numbers of copies of the variant found on plus and minus strand sequences
# 		so the true numbers of events from each line are : r1p + r1m
 
# call as _varsbuild_freqs.py fileNameRoot [true/false]


theLine = []
writeIndivLists = "false"
astrList = []
if(len(sys.argv) == 2):
	fileNameRoot = sys.argv[1]

elif(len(sys.argv) == 3):
	fileNameRoot = sys.argv[1]
	writeIndivLists = sys.argv[2]

else:
	print("varBuildFreqs.py cw 20/07/2022")
	print("call as varBuildFreqs.py fileNameSelector [true/false]")
	print("    (to also write indiv Lists add facultative parameter = true")

print("write lists = ",writeIndivLists )
if(writeIndivLists == "true") :
	os.system("mkdir __results" + fileNameRoot)

newcol = ""
thedict =  {}
theNumber = 0
theVals = []

# note the key range must be longer than the seq
keyList = [*range(1, 10001, 1)]
for i in keyList :
	thedict[int(i)] = 0

#create datafram with the keyList as first column
resTable = pd.DataFrame(keyList)
normTable = pd.DataFrame(keyList)

top = os.getcwd()  #  use current folder
for root, dirs, files in os.walk(top) :
	for fname in files:
		if  ( fileNameRoot in fname)  :
			for i in keyList :
				thedict[int(i)] = 0

			infile = open(fname, 'r')
			if(writeIndivLists == "true") :
				outfile = open("_freqs_"+fname, 'a')

			for line in infile.readlines():
				theLine = []
				theLine = line.split()
				if(theLine[0] != "type"):
					theNumber =  int(theLine[5]) + int(theLine[6] )
					thedict[int(theLine[4])] += theNumber 
					theNumber = 0
			
			newcol = fname.split(".")[0]
			theVals = list(thedict.values())
			resTable[newcol] = theVals
			normTable[newcol] = resTable[newcol]/sum(theVals)

			if(writeIndivLists == "true") :
				for key, value in thedict.items() : 
					print (key, value, sep='\t', file=outfile )

			infile.close()
			if(writeIndivLists == 1) :
				outfile.close()
	break #avoids recursive os.walk

astrList = sorted(resTable.columns ,key=str)
resTable = resTable.reindex(columns=astrList)

astrList = sorted(normTable.columns ,key=str)
normTable = normTable.reindex(columns=astrList)

print(resTable.head())

resTable.to_csv('__freqs' + fileNameRoot + '.tsv', index= False, sep='\t')
normTable.to_csv('__normfreqs' + fileNameRoot + '.tsv', index= False, sep='\t')

if(writeIndivLists == "true") :
	os.system( "mv *freqs*  __results" + fileNameRoot + "/.")

