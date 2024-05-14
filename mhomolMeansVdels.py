#!/usr/bin/env   python

import sys, os, string, datetime
import pandas as pd

# mhomolMeansVdels.py is a tidied-up version of _mhomol_meansVSdels_v2.py 
# C. White 01/02/2024

# This script parses the _all_dels.txt list file generated for the arc-plots, and outputs a table of mean nt of microhomology used for each deletion length
# output file is "__mean_Mhomol_" + fileNameRoot
# note the keyList range (deletion lengths) must be longer than the seq and is set to 1-400nt by default.
# It can be changed by including a maximum del length in the call

#call as mhomolMeansVdels.py fileNameSelector [maxDelLen] 


theLine = []
astrList = []
maxDel = 0

if(len(sys.argv) == 3):
	fileNameRoot = sys.argv[1]
	maxDel = int(sys.argv[2])

else:
	print()
	print()
	print("###################################################################")
	print("###########  mhomolMeansVdels.py  cw 01/02/2024  #############")
	print("call as mhomolMeansVdels.py fileNameSelector maxDelLen ")
	print()
	print("    note  maxDelLen  must be longer than the longest deletion. It is 400nt by default,")
	print("      but can be changed by adding the facultative second parameter")
	print("###################################################################")
	print()
	print("This script parses the _all_dels.txt list file generated for the arc-plots ")
	print("to generate a table of mean nt of microhomology used for each deletion length.")
	print()
	exit()


newcol = ""
thedict =  {}
thecounts = {}
theNumber = 0
theVals = []

# note the key range must be longer than the seq
keyList = [*range(1, maxDel+1, 1)]
for i in keyList :
	thedict[int(i)] = 0
	thecounts[int(i)] = 0

#create dataframe with the keyList as first column
resTable = pd.DataFrame(keyList)
normTable = pd.DataFrame(keyList)

top = os.getcwd()  #  use current folder
for root, dirs, files in os.walk(top) :
	for fname in files:
		if  ( fileNameRoot in fname)  :
			for i in keyList :
				thedict[int(i)] = 0
				thecounts[int(i)] = 0

			infile = open(fname, 'r')

			for line in infile.readlines():
				theLine = []
				theLine = line.split('\t')
				if( (theLine[0][0] != " ") and (theLine[0][0] != "#")  ):
					thedict[abs(int(theLine[4]))] += int(theLine[5]) 
					#print(theLine[4], theLine[5])
					thecounts[abs(int(theLine[4]))] += 1 
			newcol = fname.split(".")[0]
			theVals = list(thedict.values())
			numDels = list(thecounts.values())
			resTable[newcol] = theVals
			normTable[newcol] = [0 if j==0 else i/j for i, j in zip(theVals, numDels)]

			infile.close()
	break #avoids recursive os.walk

astrList = sorted(resTable.columns ,key=str)
resTable = resTable.reindex(columns=astrList)

astrList = sorted(normTable.columns ,key=str)
normTable = normTable.reindex(columns=astrList)

print(resTable.head())
normTable.to_csv('_delLen_vs_mean_Mhomol_' + fileNameRoot + '.tsv', index= False, sep='\t')


