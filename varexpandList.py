#!/usr/bin/env   python

import sys, os, string, datetime

# varExpandList.py is a tidied-up version of _vars_expand_List_v4.py  C.White 31/03/2023

# this script expands lists of deletions produced by CallVariants.sh script to give one deletion per line
# each line is formatted as: Type (DEL)   Chromosome   Start   End   Length   nt_microhomology
# call as varExpandList.py fileNameSelector [numLines minDel]
#      if numLines and minDel are specified
#             - generates 2 files per input file - one with all deletions >= minDelLen 
#                  and one with the specified number of dels >= minDelLen (from a randomised list of all dels) 
#      if only the filenameselector is specified
#             - generates a list with all dels (in this case minDelLen = 1)
#      includes the  number of lines written in the filename (only changes if there are less than numLines deletions)

theLine = []
numOccur = numLines = limitNumLines = minDel = linesWritten = 0
outfilename = limitOutFileName = ""

if(len(sys.argv) == 2):
	fileNameRoot = sys.argv[1]
elif(len(sys.argv) == 4):
	fileNameRoot = sys.argv[1]
	numLines = int(sys.argv[2])
	minDel = int(sys.argv[3])

	limitNumLines = 1
	print()
	print("lines = " , numLines, ";  min deletion length = " , minDel)
	print()

else:
	print("_vars_expand_List_v4.py cw 28/11/2022")
	print("call as _vars_expand_List_v4.py fileNameSelector [numLines minDel]")
	print("")
	sys.exit()

top = os.getcwd()  #  use current folder
for root, dirs, files in os.walk(top) :
	for fname in files:
		if (( fileNameRoot in fname) ) :
			infile = open(fname, 'r')
			outfilename = "_exp_" + fname

			outfile = open(outfilename, 'a')
			
			if ("_ins"  in fname) :
				print ("#type" , "Chr" , "start", "end", "length", "ins_seq", sep='\t', file=outfile ) 
			elif ("_sub" in fname) :
				print ("#type" , "Chr" , "start", "end", "length", "sub_nt", sep='\t', file=outfile ) 
			elif ("_del"  in fname) :
				print ("#type" , "Chr" , "start", "end", "length", "mhomol", sep='\t', file=outfile ) 

			linesWritten = 0

			for line in infile.readlines():
				theLine = []
				theLine = line.split("\t")
				numOccur = 0
				if( (theLine[0] != "type") and (len(theLine) == 8) and (int(theLine[4]) >= minDel) ):
					numOccur =  int(theLine[5]) + int(theLine[6] )
					linesWritten = linesWritten + numOccur
					while numOccur :
						print (theLine[0], theLine[1], theLine[2], theLine[3], theLine[4], theLine[7], sep='\t', end='', file=outfile ) 
						numOccur -= 1

			infile.close()
			outfile.close()

			if limitNumLines :
				if (numLines <= linesWritten) :
					linesWritten = numLines
				else :
					print ("only " + str(linesWritten) + "  lines written for file : " + outfilename)

				limitOutFileName = "_exp_" + str(linesWritten) + "_" + fname
				os.system("shuf -n " + str(linesWritten) + " " + outfilename  + " > " + limitOutFileName )
				print("wrote ", limitOutFileName)


	break #avoids recursive os.walk

