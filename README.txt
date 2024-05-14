############################################# 
 Charles I. White. 06 march 2024
 CNRS UMR 6293, Faculté de Médecine, Clermont-Ferrand, France
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
