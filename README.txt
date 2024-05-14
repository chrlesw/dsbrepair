############################################# 
 Charles I. White. 06 march 2024
 CNRS UMR 6293, Inserm U1103, 
 Faculté de Médecine,Université Clermont Auvergne,
 Clermont-Ferrand, France
#############################################

A set of Python and R scripts for analyses of DSB repair from paired-end (Illumina type) or long amplicon (Oxford Nanopore) DNA sequences.

In-house Python and R scripts using publicly available tools were written to analyse the DNA sequencing data produced in this study. Briefly, the scripts analyse NGS DNA sequence outputs from both Illumina paired-end sequencing or Oxford Nanopore sequencing. The Nanopore sequences are require no pre-processing, but prior to running the analyses, the paired-end Illumina sequences are merged into individual sequences using the BBMerge tool of the BBTools package (BBMap – Bushnell B. – sourceforge.net/projects/bbmap/ Brian Bushnell (Bushnell et al, 2017) ). These merged sequence files are the input for two main Python scripts, which carries out the analysis following two approaches:

The two principal scripts take different approaches to these analyses:

1) _dsb_sortAnal.py 
 Based on text searches and precise sequence lengths, this approach is best applied only to (merged) paired-end Illumina sequences. The merged Illumina sequences are individually scanned for the presence/absence of the Cas12a target sequence(s) and the presence of inversions using Python text-search tools. This, combined with their lengths, is used to sort and enumerate the merged sequences of different deletion classes: single site deletions, individual two-site deletions, deletions spanning the two cut-sites and inversions of the sequences between the two cut-sites. 

2) _dsb_varAnal.py
Minimap2 (Li H (2018) Bioinformatics 34: 3094–3100) is used to align the Nanopore or merged paired-end Illumina sequences to the reference genome, secondary or supplementary alignments are removed and the BBMap callvariants.sh script of the BBTools package  (BBMap – Bushnell B. – sourceforge.net/projects/bbmap/ ) used to produce lists of deletions, insertions and substitutions.  Shell commands then select variants affecting user-defined windows centered on the cut-sites and sort them into lists of deletions, insertions and substitutions. The presence (and possible implication) of DNA sequence microhomologies is determined for each deletion junction and plots showing each deletion as an arc linking its start and end points are generated.
 

 More detailed information is included in the script files. 

 The other scripts: _dsb_varAnal.py calls (and needs) the following scripts which need to be in the $PATH
      - varScorer.py  
      - varBuildFreqs.py  
      - varexpandList.py 
      - junctMhomol.py
      - delsArcs.R  (the full path to this is passed as an argument)
      - mhomolMeansVdels.py 

These scripts have been writted and tested to run on Python 3.11 and R 4.2.2. 

Other dependencies are:

BBMap scripts from sourceforge https://sourceforge.net/projects/bbmap/ need to be installed and added to $PATH
minimap2
samtools
biopython
pandas
coreutils
openjdk

#############################################
##########        outputs       #############
#############################################

The following folders will be found in the _alignments_variants or _alignments_variants_ont output folder
 note that r1p = number of times found on plus strand 
           r1m = number of times found on minus strand

_aln_files
      - contains ".paf" alignments of the sequence files to the reference
_arcs
     - ".png" images of the Arc plots of the deletions affecting the cut-site windons
     - tables with one deletion per line
     - _delLen_vs_mean_Mhomol__exp.tsv gives, for each deletion length, the mean number 
          of nt of microhomology potentially involved in the deletions 
          (note that potential microhomologies are scored up to 9 nt, so the "9 nt" class will include any deletions with more than 9)

_bam_files
     - sorted, indexed BAM files with the mapping of the sequence files to the reference

_junctions
     - _delsList.txt : 
          Description of each deletion (one per line, sorted by length), showing:
               Chromosome ; Start nt ; End nt ; deletion length ; r1p ; r1m
    - _delsJunctions.txt : 
          For each analysed sequence file, shows the numbers of deletions with a given number of nt of microhomology (these are also grouped together in the _stats_mhomol_Junctions.txt file
          This is followed by a tab-separated, detailed description of each deletion, showing:
               "10_nt_upstream ( first and last 10 nt of deletion) 10_nt_downstream" ; deletion length ; nt of microhomology ; Chromosome ; Start nt ; End nt ; r1p ; r1m
          (note that potential microhomologies are scored up to 9 nt, so the "9 nt" class will include any deletions with more than 9)
     - _stats_mhomol_Junctions.txt summarises the number and the proportion of deletions afecting the cut-site windows  
          (each line extracted from the corresponding _delsJunctions.txt file) 

_vars
     - _vars.txt : output of variant caller script
     - _indelsList.txt: table of deletions, insertions and substitutions affecting the cut-site windows.
     - files with tables of deletions, insertions and substitutions affecting the cut-site windows. 3 tables are given for 2-cut-site analyses:  all dels (ins, sub) and also individual tables for dels (ins, sub) affecting target_t1 or target_t2.
     - _mHomolDelsList.txt information for each deletion affecting the cut-site windows:
            DEL     Chr  start     end  length    r1p  r1m  length_mhomol
     - __freqs_all.tsv (...): 
          summary tables of numbers of deletions, insertions and substitutions affecting the cut-site windows, for each deletion or insertion length. These are also fgiven as fractions of total events in each class (__normfreqs_all.tsv). For 2-cut-site analyses include tables with all dels (ins, sub) and also individual tables for dels (ins, sub) affecting target_t1 or target_t2.

_vcf_files : indexed vcf files 

#############################################
#############################################
