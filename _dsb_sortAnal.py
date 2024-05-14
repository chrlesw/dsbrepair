#!/usr/bin/env   python

#####################  _dsb_sortAnal.py    ######################################
###### C. White 
###### 31/10/2023.

# _dsb_sortAnal.py is a tidied-up version of _cas_1_or_2_sites_sortNstats.py
# It uses amplicon sequence length and python text matching to count and sort sequences with different classes of modifications of their cut-sites into seperate files ...
# If used with F1 hybrids, the script will analyse/sort the deletions/inversions in the "genome1" sequences", and also report on those in "genome2"
# this is done by using the genome-specific"U1" sequences (to be chosen to include a sequence polymorphism between the two genomes)"


# the script should be run from the folder containing the fastq files output from sequencing a CRISPR-Cas experiment.
# note that fastq.gz files need to be uncompressed beforehand

#############################################################################################################
# The parameters are:

# -numSites or --number_of_cut_sites : type=int, required=True, choices=[1,2], help="number of target sites (1 or 2)")

# -maxSeqLen or --maximum_len_seq : required=True, type=int, help="maximum_length of seq to analyse")
# -minSeqLen or --minimum_len_seq : required=True, type=int, help="minimum_length of seq to analyse")

# -targGenome1_1 or --site_1_seq_genome1 : required=True,  help="sequence to be tested at cut site 1 in Genome1")
# -targGenome1_2 or --site_2_seq_genome1 : required=False,  help="sequence to be tested at cut site 2 in Genome1")
# -targGenome2_1 or --site_1_seq_genome2 : required=False,  help="sequence to be tested at cut site 1 in Genome2")
# -targGenome2_2 or --site_2_seq_genome2 : required=False,  help="sequence to be tested at cut site 2 in Genome2")

# -fNameRoot or --filenameroot : required=True, help="string common to names of files to be analysed")
# -doubleCutDelMax or --double_cut_del_len : required=False, type=int, help="maximum_length of seq with 2-site deletion")
# -invSeq or --middleInvSeq : required=False, default='None', help="sequence to check for inversion")
# -u1Genome1 or --seq_u1_genome1 : required=True, default='None', help="Sequence to identify genome 1")
# -u1Genome2 or --seq_u1_genome2 : required=False, default='None', help="Sequence to identify genome 2")

#############################################################################################################
# for example for a 400 bp amplicon from site2_915 :
  # _dsb_sortAnal.py  -numSites 2 -maxSeqLen 398 -minSeqLen 50 -targGenome1_1 TTTGAGCTCGACCAAGTTCCCTTTGTTTTC  -targGenome1_2 TTTTGTTGGAGAAGATATGCTCACTGATCTC \
  # -targGenome2_1 aseq  -targGenome2_2 aseq  -doubleCutDelMax 200 -invSeq GACTGTTACTCCATGGGATC  -fNameRoot _pooled -u1Genome1  GCATCTGGAATTG -u1Genome2 TTGGTTGTTG
#############################################################################################################

import sys, os, string, datetime, argparse, gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

print()
print("################  _dsb_sortAnal.py    ############################")
print()
print("######## C. White") 
print("######## 31/10/2023")
print("#################")
print("This script analyses fastq files (fastq.gz files need to be uncompressed beforehand) for seqs modified at the cut-site(s) using python text matching.")
print("   It counts and sorts seqs into separate files based on the modifications of their cut-sites ...") 
print("     and will also take into account the genome in F1 hybrids")
print("  it should be run from the folder containing the fastq (not fastq.gz) files output from sequencing a CRISPR-Cas experiment. ")
print()
print()

numOfSites = maxLength = minLength = fullDelLength =  0
totalSeqs = totalFiles = 0
num_genome1_chr = num_genome2_chr = 0

#################################################
parser = argparse.ArgumentParser()

parser.add_argument("-numSites", "--number_of_cut_sites", type=int, required=True, choices=[1,2], help="number of target sites (1 or 2)")

parser.add_argument("-maxSeqLen", "--maximum_len_seq", required=True, type=int, help="maximum_length of sequences to analyse")
parser.add_argument("-minSeqLen", "--minimum_len_seq", required=True, type=int, help="minimum_length of sequences to analyse")

parser.add_argument("-targGenome1_1", "--site_1_seq_genome1", required=True, default='aseq', help="sequence to be tested at cut site 1 in Genome1")
parser.add_argument("-targGenome1_2", "--site_2_seq_genome1", required=False, default='aseq',  help="sequence to be tested at cut site 2 in Genome1")
parser.add_argument("-targGenome2_1", "--site_1_seq_genome2", required=False, default='aseq',  help="sequence to be tested at cut site 1 in Genome2")
parser.add_argument("-targGenome2_2", "--site_2_seq_genome2", required=False, default='aseq', help="sequence to be tested at cut site 2 in Genome2")

parser.add_argument("-fNameRoot", "--filenameroot", required=True, help="string common to names of files to be analysed")
parser.add_argument("-doubleCutDelMax", "--double_cut_del_len", required=False, type=int, help="maximum_length of sequences with 2-site deletion")
parser.add_argument("-invSeq", "--middleInvSeq", required=False, default='None', help="sequence to check for inversions")
parser.add_argument("-u1Genome1", "--seq_u1_genome1", required=True, default='None', help="Genome1 U1 seq to distinguish between Genome1 and Genome2 sequences")
parser.add_argument("-u1Genome2", "--seq_u1_genome2", required=False, default='None', help="Genome2 U1 seq to distinguish between Genome1 and Genome2 sequences")
##################
args = parser.parse_args()

numOfSites = args.number_of_cut_sites
maxLength = args.maximum_len_seq
minLength = args.minimum_len_seq
target_1_genome1_seq = args.site_1_seq_genome1
target_2_genome1_seq =args.site_2_seq_genome1
target_1_genome2_seq = args.site_1_seq_genome2
target_2_genome2_seq = args.site_2_seq_genome2

fileNameRoot = args.filenameroot
fullDelLength = args.double_cut_del_len
middleInvertedSeq = args.middleInvSeq
seq_u1_genome1 = args.seq_u1_genome1
seq_u1_genome2 = args.seq_u1_genome2
#################################################

depth = 0

today = datetime.date.today()
theDate = str(f'{today:%Y_%m_%d}')
theOutputFolder = "_dels_output_" + theDate + "/"

theMergedOutputFolder = os.getcwd() + "/"
top = os.getcwd()

try:
    os.makedirs(theOutputFolder)
except OSError:
    print("the output folder already exists")
    sys.exit()

f_out_stats = open(theOutputFolder + "_theStats_" + theDate + ".txt" , 'a') #note - opened with "a" for appending

if numOfSites == 1 :
    print("#  " + str(today) + "  Single cut site del analyses: seq maxLength = " + str(maxLength) + "  ; minLength = " + str(minLength), file=f_out_stats )
    print("#  Target_genome1 = ", target_1_genome1_seq, "Target_genome2 = ", target_1_genome2_seq,   file=f_out_stats )
    print("#name" , "seqs_read" , "% dels_genome1_all"  , "% dels_hors_site" , "% dels_genome2" , "% dels_nonGenome1_nonGenome2" , "% dels_target_1" , \
     "dels_genome1_all", "no_targ_del", "dels_nonGenome1_nonGenome2" , "dels_target_1"  , sep="\t" , file=f_out_stats)
elif numOfSites == 2 :
    print("#  " + str(today) + "  Double cut site del analyses: seq maxLength = " + str(maxLength) + "  ; minLength = " + str(minLength) + "  ; max_full_deletion = " + str(fullDelLength) , file=f_out_stats)
    print("#  Target_1_genome1 = ", target_1_genome1_seq, "Target_1_genome2 = ", target_1_genome2_seq, "  Target_2_genome1 = ", target_2_genome1_seq, "Target_2_genome2 = ", target_2_genome2_seq, "middleInvertedSeq = ", middleInvertedSeq, file=f_out_stats )
    print("#name" , "seqs_read" , "% dels_genome1_all"  , "% dels_hors_site" , "% dels_genome2" , "%dels_nonGenome1_nonGenome2", "% dels_target_1_only" , "% dels_target_2_only" ,\
     "% dels_targets_1_and_2" , "% dels_targets_1_and_2_indiv" , "% dels_targets_1_to_2" , "% dels_targets_1_to_2_inversion"  , \
     "% dels_1_and_2_indiv/1and2" , "% dels_targets_1_to_2/1and2" , "% dels_targets_1_to_2_inversion/1and2"  , \
     "dels_genome1_all", "no_targ_del",  "dels_genome2" , "dels_nonGenome1_nonGenome2", "dels_target_1_only" , "dels_target_2_only" , "dels_targets_1_and_2" , \
     "dels_targets_1_and_2_indiv" , "dels_targets_1_to_2" , "dels_targets_1_to_2_inversion"  , sep="\t" , file=f_out_stats)
######################################################
totalSeqs = 0
totalFiles = 0

for root, dirs, files in os.walk(top) :
    for fname in files:
        fileName = fname.split(".")
        baseName = fileName[0].split(fileNameRoot)[0]
        num_genome1_chr = num_genome2_chr = 0

        ############### for a single cut-site  ###################

        if (( fileNameRoot in fname) and (numOfSites == 1) and (".fastq" in fname)) :
            no_targ_del = dels_target_1 = seqs_read = dels_genome1_all = dels_genome2 = dels_nonGenome1_nonGenome2 = 0
            totalFiles +=1
            f_out_targ_1_dels = open(theOutputFolder + baseName + "_dels_targ_1.fastq", 'w')
            f_out_genome2_dels = open(theOutputFolder + baseName + "_dels_genome2.fastq", 'w')
            f_out_dels_nonGenome1_nonGenome2 = open(theOutputFolder + baseName + "_dels_nonGenome1_nonGenome2.fastq", 'w')

            with open(theOutputFolder + baseName + "_dels_none.fastq", 'w') as f_out_no_targ_del :
                for seq_record in SeqIO.parse(open(theMergedOutputFolder + fname, mode='r'), 'fastq'):
                    seqs_read +=1
                    totalSeqs +=1
                    is_genome1_chr = 0
                    if   (seq_u1_genome1 in seq_record.seq):
                        num_genome1_chr +=1
                    if   (seq_u1_genome2 in seq_record.seq):
                        num_genome2_chr +=1
                        
                    if  ( (seq_u1_genome1 in seq_record.seq) and (len(seq_record) <= maxLength) and (len(seq_record) >= minLength) ) :
                        num_genome1_chr +=1
                        if not ((target_1_genome1_seq in seq_record.seq) or (target_1_genome2_seq in seq_record.seq) )  :
                            dels_genome1_all+=1
                            dels_target_1 +=1
                            r=SeqIO.write(seq_record, f_out_targ_1_dels, 'fastq')
                            if r!=1: print('Error while writing sequence:  ' + seq_record.id)
                        else :
                            no_targ_del+=1
                            r=SeqIO.write(seq_record, f_out_no_targ_del, 'fastq')
                            if r!=1: print('Error while writing sequence:  ' + seq_record.id)
                    elif ( (seq_u1_genome2 in seq_record.seq) and (len(seq_record) <= maxLength) and (len(seq_record) >= minLength) ) :
                        dels_genome2 +=1
                        r=SeqIO.write(seq_record, f_out_genome2_dels, 'fastq')
                        if r!=1: print('Error while writing sequence:  ' + seq_record.id)
                    elif ( (len(seq_record) <= maxLength) and (len(seq_record) >= minLength) ) :
                        dels_nonGenome1_nonGenome2 +=1
                        r=SeqIO.write(seq_record, f_out_dels_nonGenome1_nonGenome2, 'fastq')
                        if r!=1: print('Error while writing sequence:  ' + seq_record.id)        

            if (seqs_read != 0) :
                print(baseName , seqs_read , 100*dels_genome1_all/seqs_read , 100*no_targ_del/seqs_read , 100*dels_genome2/seqs_read  , 100*dels_nonGenome1_nonGenome2/seqs_read , 100*dels_target_1/seqs_read  , \
                 dels_genome1_all, no_targ_del , dels_nonGenome1_nonGenome2, dels_target_1  ,  sep="\t" , file=f_out_stats )
            else :
                print(baseName , seqs_read , "" , "" , ""  , ""  , "" , \
                 dels_genome1_all, no_targ_del , dels_nonGenome1_nonGenome2, dels_target_1  ,  sep="\t" , file=f_out_stats )

            print(fname)
            print(seqs_read, " sequences read")
            if (seqs_read > 0):
                print(dels_genome1_all, " genome1 sequences with deletion (", 100*(dels_genome1_all/seqs_read), " %")
                print(dels_genome2, " genome2 sequences with deletion (", 100*(dels_genome2/seqs_read), " %")
            else :
                print(dels_genome1_all, " genome1 sequences with deletion")
                print(dels_genome2, " genome2 sequences with deletion")

            print(dels_nonGenome1_nonGenome2, " other sequences with deletion")
            print()

        ############### for double cut-sites  ###################

        elif (( fileNameRoot in fname) and (numOfSites == 2) and (".fastq" in fname)) :
            no_targ_del = dels_target_1_only = dels_target_2_only = dels_targets_1_and_2 = dels_target_3 = 0
            seqs_read = dels_genome1_all = dels_targets_1_to_2 = dels_targets_1_to_2_inversion = dels_targets_1_and_2_indiv = 0
            dels_genome1 = dels_genome2 = dels_genome1_all = dels_rec = dels_nonGenome1_nonGenome2 = 0 

            totalFiles +=1

            f_out_genome2_dels = open(theOutputFolder + baseName + "_dels_genome2.fastq", 'w')
            f_out_dels_nonGenome1_nonGenome2 = open(theOutputFolder + baseName + "_dels_nonGenome1_nonGenome2.fastq", 'w')

            f_out_targ_1_dels = open(theOutputFolder + baseName + "_dels_targ_1.fastq", 'w')
            f_out_targ_2_dels = open(theOutputFolder + baseName + "_dels_targ_2.fastq", 'w')
            f_out_1and2_dels = open(theOutputFolder + baseName + "_dels_targets_1_and_2.fastq", 'w')
            f_out_1to2_dels = open(theOutputFolder + baseName + "_dels_targets_1_to_2.fastq", 'w')
            f_out_1and2_indiv_dels = open(theOutputFolder + baseName + "_dels_targets_1_and_2_indiv.fastq", 'w')
            f_out_1to2_inversions = open(theOutputFolder + baseName + "_dels_targets_1_to_2inversions.fastq", 'w')
            f_out_no_targ_del = open(theOutputFolder + baseName + "_dels_neither.fastq", 'w')

            with  f_out_no_targ_del :
                
                for seq_record in SeqIO.parse(open(theMergedOutputFolder + fname, mode='r'), 'fastq'):
                    seqs_read +=1
                    totalSeqs +=1
                    is_genome1_chr = 0
                    if   (seq_u1_genome1 in seq_record.seq):
                        num_genome1_chr +=1
                    if   (seq_u1_genome2 in seq_record.seq):
                        num_genome2_chr +=1

                    if  ( (seq_u1_genome1 in seq_record.seq) and (len(seq_record) <= maxLength) and (len(seq_record) >= minLength) ) :
                        if not  ((target_1_genome1_seq in seq_record.seq) or (target_1_genome2_seq in seq_record.seq) ) :
                            dels_genome1_all+=1
                            if ( (target_2_genome1_seq in seq_record.seq) or (target_2_genome2_seq in seq_record.seq ) ) :
                                dels_target_1_only +=1
                                r=SeqIO.write(seq_record, f_out_targ_1_dels, 'fastq')
                                if r!=1: print('Error while writing sequence:  ' + seq_record.id)
                            else :
                                dels_targets_1_and_2 +=1
                                r=SeqIO.write(seq_record, f_out_1and2_dels, 'fastq')
                                if r!=1: print('Error while writing sequence:  ' + seq_record.id)
                                if ( (len(seq_record) <= fullDelLength) ) :
                                    dels_targets_1_to_2 +=1
                                    r=SeqIO.write(seq_record, f_out_1to2_dels, 'fastq')
                                    if r!=1: print('Error while writing sequence:  ' + seq_record.id)
                                elif ( (middleInvertedSeq in seq_record.seq) ) :
                                    dels_targets_1_to_2_inversion +=1
                                    r=SeqIO.write(seq_record, f_out_1to2_inversions, 'fastq')
                                    if r!=1: print('Error while writing sequence:  ' + seq_record.id)
                                else :
                                    dels_targets_1_and_2_indiv +=1
                                    r=SeqIO.write(seq_record, f_out_1and2_indiv_dels, 'fastq')
                                    if r!=1: print('Error while writing sequence:  ' + seq_record.id)

                        elif not (target_2_genome1_seq in seq_record.seq) :
                            dels_genome1_all+=1
                            dels_target_2_only +=1
                            r=SeqIO.write(seq_record, f_out_targ_2_dels, 'fastq')
                            if r!=1: print('Error while writing sequence:  ' + seq_record.id)

                        else :
                            no_targ_del+=1
                            r=SeqIO.write(seq_record, f_out_no_targ_del, 'fastq')
                            if r!=1: print('Error while writing sequence:  ' + seq_record.id)

                    elif ( (seq_u1_genome2 in seq_record.seq) and (len(seq_record) <= maxLength) and (len(seq_record) >= minLength) ) :
                        dels_genome2 +=1
                        r=SeqIO.write(seq_record, f_out_genome2_dels, 'fastq')
                        if r!=1: print('Error while writing sequence:  ' + seq_record.id)
                    elif ( (len(seq_record) <= maxLength) and (len(seq_record) >= minLength) ) :
                        dels_nonGenome1_nonGenome2 +=1
                        r=SeqIO.write(seq_record, f_out_dels_nonGenome1_nonGenome2, 'fastq')
                        if r!=1: print('Error while writing sequence:  ' + seq_record.id)                    	

            if ( (dels_targets_1_and_2 != 0) and (seqs_read != 0) ) :
                print(baseName , seqs_read , 100*dels_genome1_all/seqs_read , 100*no_targ_del/seqs_read , 100*dels_genome2/seqs_read ,100*dels_nonGenome1_nonGenome2/seqs_read , 100*dels_target_1_only/seqs_read  , 100*dels_target_2_only/seqs_read ,\
                    100*dels_targets_1_and_2/seqs_read , 100*dels_targets_1_and_2_indiv/seqs_read , 100*dels_targets_1_to_2/seqs_read , 100*dels_targets_1_to_2_inversion/seqs_read ,  \
                    100*dels_targets_1_and_2_indiv/dels_targets_1_and_2 , 100*dels_targets_1_to_2/dels_targets_1_and_2 , 100*dels_targets_1_to_2_inversion/dels_targets_1_and_2 ,  \
                     dels_genome1_all, no_targ_del , dels_genome2 , dels_nonGenome1_nonGenome2, dels_target_1_only , dels_target_2_only , dels_targets_1_and_2 , dels_targets_1_and_2_indiv  ,  dels_targets_1_to_2 , \
                    dels_targets_1_to_2_inversion ,  sep="\t" , file=f_out_stats )
            elif (seqs_read != 0) :
                print(baseName , seqs_read , 100*dels_genome1_all/seqs_read , 100*no_targ_del/seqs_read , 100*dels_genome2/seqs_read , 100*dels_nonGenome1_nonGenome2/seqs_read ,100*dels_target_1_only/seqs_read  , 100*dels_target_2_only/seqs_read ,\
                    100*dels_targets_1_and_2/seqs_read , 100*dels_targets_1_and_2_indiv/seqs_read , 100*dels_targets_1_to_2/seqs_read , 100*dels_targets_1_to_2_inversion/seqs_read ,  \
                     "", "" , "" , \
                     dels_genome1_all, no_targ_del , dels_genome2 , dels_nonGenome1_nonGenome2, dels_target_1_only , dels_target_2_only , dels_targets_1_and_2 , dels_targets_1_and_2_indiv  ,  dels_targets_1_to_2 , \
                    dels_targets_1_to_2_inversion ,  sep="\t" , file=f_out_stats )
            else :
                print(baseName , seqs_read , "" , "" , "" , ""  , "" ,\
                     "", "" , "" , "" , "" \
                     , "" , "" , "", \
                     dels_genome1_all, no_targ_del , dels_genome2 , dels_nonGenome1_nonGenome2, dels_target_1_only , dels_target_2_only , dels_targets_1_and_2 , dels_targets_1_and_2_indiv  ,  dels_targets_1_to_2 , \
                    dels_targets_1_to_2_inversion ,  sep="\t" , file=f_out_stats )

            print(fname)
            print(seqs_read, " sequences read")
            if (num_genome1_chr >0):
                print(num_genome1_chr, " genome1 sequences found, ", dels_genome1_all, " with deletion (", 100*(dels_genome1_all/num_genome1_chr), " %")
            else : 
                print(num_genome1_chr, " genome1 sequences found, ", dels_genome1_all, " with deletion")

            if (num_genome2_chr >0):
                print(num_genome2_chr, " genome2 sequences found, ", dels_genome2, " with deletion (", 100*(dels_genome2/num_genome2_chr), " %")
            else:
                print(num_genome2_chr, " genome2 sequences found, ", dels_genome2, " with deletion")

            print(dels_nonGenome1_nonGenome2, " other sequences with deletion")
            print()


    break #avoids recursive os.walk
    ######################################################

theStatsFile = f_out_stats.name
f_out_stats.close()

os.system("sort " + theStatsFile + " -o " + theStatsFile)

print()

print("tidying up")

os.chdir(top)
os.chdir(theOutputFolder)

try:
    os.makedirs("_fastq_sorted")
except OSError:
    print("the output folder already exists")
    sys.exit()

os.system("mv *.fastq  _fastq_sorted/. ")

print("*******************************************************")
print(totalFiles, " files analysed")
print(totalSeqs, " sequences analysed")
print("*******************************************************")
print()

##################################################################################################################################
