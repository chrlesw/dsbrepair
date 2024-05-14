#!/usr/bin/env Rscript

# delsArcs.R is a tidied-up version of arcdiagram_del_mhomol_v5.R
#  C. White 19/10/2023. (based on an original idea of Miguel Hernandez Sanchez-Rebato)

## this script plots arcs joining the start and end of deletion breakpoints, 
## the lines can be coloured according to the number of nt of microhomology potentially involved in the deletion 
##   but this is less informative than originally thought and in practice is normally not done.
##
# for input the script requires deletion list text files with one event per line in format:
# "Type", "Chromosome", "Start", "End", "Length", "nt mhomol"

 # arguments in the call : theDir isNpore targ inclFocus numsites plotMhomol


# Different target sites from our current work are left in as examples, but users will need to add or edit to include their own targets. 

#######################################################################################
###########     load packages  ########################
library(conflicted)  #needed because of conflicts between tidyverse and stats packages

library("tidyverse")
conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")

library("ggraph")
library("tidygraph")
library("ggpubr")
########################################################

############  set theDir to the working directory #########################
#############  set fileroot to the regular expression for selecting files in the working directory ############  
# # # pattern examples: "" all files
#     "wt" all filenames containing "wt"
#     "*_dels\\.txt$"  all filenames ending with "_dels.txt"

fileroot <- paste("*_exp" , ".*\\.txt$", sep = "")  
#fileroot <- paste("*_exp_" , ".*nt\\.txt$", sep = "")

theDir <- "/Users/self/Desktop/temp"
########################################################
# These parameters can be set here if the script is run directly (eg in RStudio). 
# They will be overridden by the values passed with the call if the script is called via the shell (or from python...)

isNpore <- "yes"
#isNpore <- "no"


#targ <- "site2_915"
targ <- "chr3_92"
# targ <- "ercc1_2ab"
# targ <- "ercc1_2cd"
# targ <- "xpf_2ab"
# targ <- "send1_2ab"
# targ <- "iugus"
# targ <- "chr2_184"

inclFocus <- "no"
#inclFocus <- "yes"

#numsites <- "1"
numsites <- "2"

plotMhomol = "no"
#plotMhomol = "yes"

########################################################
########################################################

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  # theDir <- "_alignments_variants_ont/_arcs"
} else if (length(args)==1) { 
  theDir <- args[1]
} else if (length(args)==6) { 
  theDir <- args[1]
  isNpore <- args[2]
  targ <- args[3]
  inclFocus <- args[4]
  numsites <- args[5]
  plotMhomol <- args[6]
}

########################################################
########################################################

setwd(theDir)
all_files <- list.files(path = theDir, 
                        pattern = fileroot , 
                        full.names = TRUE)

deletions_theme <-  theme(axis.line.x=element_blank(),
                          axis.line.y=element_blank(),
                          axis.text.y=element_blank(),
                          axis.ticks.x=element_line(),
                          axis.ticks.y=element_blank(),
                          axis.title.x=element_text(family="Helvetica Neue Light", size=9, face="plain"),
                          axis.title.y=element_blank(),
                          panel.background=element_blank(),
                          panel.spacing.x = unit(0.01, "lines"),
                          plot.margin = unit(c(1, 1, -0.5, 1), "lines"),
                          panel.border=element_blank(),
                          panel.grid.major=element_blank(),
                          panel.grid.minor=element_blank(),
                          aspect.ratio = 1 / 3, 
                          
                          plot.background=element_blank(),
                          plot.title = element_text(family="Helvetica Neue Light", size=20, face="plain"),
                          plot.subtitle =  element_text(family="Helvetica Neue Light", size=12, face="plain"),
                          legend.text =  element_text(family="Helvetica Neue Light", size=7, face="plain"),
                          legend.title = element_text(family="Helvetica Neue Light", size=8, face="plain"),
                          legend.direction = "horizontal",
                          legend.position = "top", 
                          #axis.text = element_text(family="Helvetica Neue Light", face="plain", size=7),
                          axis.text.x = element_text(family="Helvetica Neue Light", face="plain", size=5, angle=45,  vjust=0.5),
                          legend.key = element_blank())

########################################################

# process each file:
lapply(X = all_files,
       FUN = function(path) {
         # read your data
         df <- read.table(path)

         # do your transformation
         names(df)=c( "Type", "Chromosome", "Start", "End", "delLength", "Microhomology")
         
         dels_edge <- df %>%
           select(Start, End, Microhomology) %>%
           group_by(Start, End, Microhomology) %>%
           na.omit()
         
         dels_tidy <- tbl_graph(edges = dels_edge, directed = FALSE)
         
#############################  NANOPORE  ###############################
if (isNpore == "yes") {
    if (targ == "site2_915"){
        if (plotMhomol == "no") {
            df_out <- ggraph(dels_tidy, layout = 'linear')  + geom_vline(xintercept = c(5187987, 5188243), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(5186670, 5189470)) + scale_x_continuous(breaks=seq(5186670, 5189470, 100))   
        } else if (plotMhomol == "2colours") { 
          df_out <- ggraph(dels_tidy, layout = 'linear')  + geom_vline(xintercept = c(5187987, 5188243), colour="darkgoldenrod1") + geom_edge_arc(strength = 1,  aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(5186670, 5189470)) + scale_x_continuous(breaks=seq(5186670, 5189470, 100)) + scale_edge_colour_gradient(limits = c(0,2), low = "blue", high = "blue" , na.value = "red")
        } else if (plotMhomol == "testing") { 
          df_out <- ggraph(dels_tidy, layout = 'linear')  + geom_vline(xintercept = c(5187987, 5188243), colour="darkgoldenrod1") + geom_edge_arc(strength = 1,  aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(5186670, 5189470)) + scale_x_continuous(breaks=seq(5186670, 5189470, 100)) + scale_edge_colour_gradient(limits = c(0,9), low = "red", high = "blue" , na.value = "black")
        } else { 
          df_out <- ggraph(dels_tidy, layout = 'linear')  + geom_vline(xintercept = c(5187987, 5188243), colour="darkgoldenrod1") + geom_edge_arc(strength = 1,  aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(5186670, 5189470)) + scale_x_continuous(breaks=seq(5186670, 5189470, 100)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red")
        }
      
      if (inclFocus == "yes") {
        df_out_focus_site1 <- ggraph(dels_tidy, layout = 'linear')  + geom_vline(xintercept = c(5187987, 5188243), colour="darkgoldenrod1") + geom_edge_arc(strength = 1,  aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(5187920, 5188050)) + scale_x_continuous(breaks=seq(5187920, 5188050, 5)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red")
        df_out_focus_site2 <- ggraph(dels_tidy, layout = 'linear')  + geom_vline(xintercept = c(5187987, 5188243), colour="darkgoldenrod1") + geom_edge_arc(strength = 1,  aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(5188180, 5188350)) + scale_x_continuous(breaks=seq(5188170, 5188300, 5)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red")
      }
    }

    else if (targ == "chr3_92"){
        if (plotMhomol == "no") {
          df_out <- ggraph(dels_tidy, layout = 'linear')  + geom_vline(xintercept = c(12113733, 12113355), colour="darkgoldenrod1") + geom_edge_arc(strength = 1,  alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(12111994, 12114973)) + scale_x_continuous(breaks=seq(12111994, 12114994, 100))  
          
        } else if (plotMhomol == "2colours") { 
          df_out <- ggraph(dels_tidy, layout = 'linear')  + geom_vline(xintercept = c(12113733, 12113355), colour="darkgoldenrod1") + geom_edge_arc(strength = 1,  aes(colour=Microhomology), alpha=.5) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(12111994, 12114973)) + scale_x_continuous(breaks=seq(12111994, 12114973, 100)) + scale_edge_colour_gradient(limits = c(0,2), low = "blue", high = "blue" , na.value = "red")
        } else if (plotMhomol == "testing") { 
          df_out <- ggraph(dels_tidy, layout = 'linear')  + geom_vline(xintercept = c(12113733, 12113355), colour="darkgoldenrod1") + geom_edge_arc(strength = 1,  aes(colour=Microhomology), alpha=.5) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(12111994, 12114973)) + scale_x_continuous(breaks=seq(12111994, 12114973, 100)) + scale_edge_colour_gradient(limits = c(0,8), low = "yellow", high = "red" , na.value = "blue")
        } else { 
          df_out <- ggraph(dels_tidy, layout = 'linear')  + geom_vline(xintercept = c(12113733, 12113355), colour="darkgoldenrod1") + geom_edge_arc(strength = 1,  aes(colour=Microhomology), alpha=.5) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(12111994, 12114973)) + scale_x_continuous(breaks=seq(12111994, 12114973, 100)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red")
        }

        if (inclFocus == "yes") {
            df_out_focus_site1 <- ggraph(dels_tidy, layout = 'linear')  + geom_vline(xintercept = c(12113733, 12113355), colour="darkgoldenrod1") + geom_edge_arc(strength = 1,  aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(12113700, 12113760)) + scale_x_continuous(breaks=seq(12113700, 12113760, 5)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red")
            df_out_focus_site2 <- ggraph(dels_tidy, layout = 'linear')  + geom_vline(xintercept = c(12113733, 12113355), colour="darkgoldenrod1") + geom_edge_arc(strength = 1,  aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(12113325, 12113385)) + scale_x_continuous(breaks=seq(12113325, 12113385, 5)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red")
        }
    }
    else if (targ == "ercc1_2ab"){
        if (plotMhomol == "no") {
            df_out <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(1479622, 1479796), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(1478104, 1481076) ) + scale_x_continuous(breaks=seq(1478104, 1481076, 50))    
        } else if (plotMhomol == "2colours") { 
            df_out <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(1479622, 1479796), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(1478104, 1481076) ) + scale_x_continuous(breaks=seq(1478104, 1481076, 50)) + scale_edge_colour_gradient(limits = c(0,2), low = "blue", high = "blue" , na.value = "red")
        }  else { 
          df_out <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(1479622, 1479796), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(1478104, 1481076) ) + scale_x_continuous(breaks=seq(1478104, 1481076, 50)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red")
        }

        if (inclFocus == "yes") {
            df_out_focus_site1 <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(1479622, 1479796), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(1478104, 1481076) ) + scale_x_continuous(breaks=seq(1478104, 1481076, 5)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red") 
            df_out_focus_site2 <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(1479622, 1479796), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(1478104, 1481076) ) + scale_x_continuous(breaks=seq(1478104, 1481076, 5)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red") 
        }    
    }
    else if (targ == "ercc1_2cd"){
    }
    else if (targ == "xpf_2ab"){
        if (plotMhomol == "no") {
            df_out <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(16476847, 16476711), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(16475300, 16478300) ) + scale_x_continuous(breaks=seq(16475300, 16478300, 50)) 
        } else if (plotMhomol == "2colours") { 
            df_out <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(16476847, 16476711), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(16475300, 16478300) ) + scale_x_continuous(breaks=seq(16475300, 16478300, 50)) + scale_edge_colour_gradient(limits = c(0,2), low = "blue", high = "blue" , na.value = "red")
        } else { 
          df_out <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(16476847, 16476711), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(16475300, 16478300) ) + scale_x_continuous(breaks=seq(16475300, 16478300, 50)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red")
        }

        if (inclFocus == "yes") {
            df_out_focus_site1 <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(16476847, 16476711), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(16476817, 16476877) ) + scale_x_continuous(breaks=seq(16476817, 16476877, 5)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red") 
            df_out_focus_site2 <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(16476847, 16476711), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(16476681, 16476741) ) + scale_x_continuous(breaks=seq(16476681, 16476741, 5)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red") 
        }    
    }
    else if (targ == "send1_2ab"){
        if (plotMhomol == "no") {
            df_out <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(18132565, 18132121), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(18130750, 18133728) ) + scale_x_continuous(breaks=seq(18130750, 18133728, 50))
        } else  if (plotMhomol == "2colours") { 
            df_out <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(18132565, 18132121), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(18130750, 18133728) ) + scale_x_continuous(breaks=seq(18130750, 18133728, 50)) + scale_edge_colour_gradient(limits = c(0,2), low = "blue", high = "blue" , na.value = "red")
        } else { 
          df_out <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(18132565, 18132121), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(18130750, 18133728) ) + scale_x_continuous(breaks=seq(18130750, 18133728, 50)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red")
        }

        if (inclFocus == "yes") {
            df_out_focus_site1 <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(18132565, 18132121), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(18130750, 18133728) ) + scale_x_continuous(breaks=seq(18130750, 18133728, 5)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red") 
            df_out_focus_site2 <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(18132565, 18132121), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(18130750, 18133728) ) + scale_x_continuous(breaks=seq(18130750, 18133728, 5)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red") 
        }    
    }
    else if (targ == "iugus"){
        #(for amplicon 2940bp)
        if (plotMhomol == "no") {
            df_out <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(713), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(1, 2940) ) + scale_x_continuous(breaks=seq(1, 2940, 100)) 
        } else  if (plotMhomol == "2colours") { 
            df_out <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(713), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(1, 2940) ) + scale_x_continuous(breaks=seq(1, 2940, 100)) + scale_edge_colour_gradient(limits = c(0,2), low = "blue", high = "blue" , na.value = "red")#, breaks=c(2,4,6,8)) 
        } else { 
          df_out <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(713), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(1, 2940) ) + scale_x_continuous(breaks=seq(1, 2940, 100)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red")#, breaks=c(2,4,6,8)) 
        }

        if (inclFocus == "yes"){
            df_out_focus_site1 <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(713), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(690, 740), ylim = c(0,10) ) + scale_x_continuous(breaks=seq(690, 740, 5)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red") #, breaks=c(2,4,6,8)) 
        }
    }
    
    else if (targ == "chr2_184"){
        if (plotMhomol == "no") {
            df_out <- ggraph(dels_tidy, layout = 'linear')  + geom_vline(xintercept = c(13655647, 13655981), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, alpha=.5) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(13654306, 13657310)) + scale_x_continuous(breaks=seq(13654306, 13657310, 100)) 
        } else  if (plotMhomol == "2colours") { 
            df_out <- ggraph(dels_tidy, layout = 'linear')  + geom_vline(xintercept = c(13655647, 13655981), colour="darkgoldenrod1") + geom_edge_arc(strength = 1,  aes(colour=Microhomology), alpha=.5) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(13654306, 13657310)) + scale_x_continuous(breaks=seq(13654306, 13657310, 100)) + scale_edge_colour_gradient(limits = c(0,2), low = "blue", high = "blue" , na.value = "red")
        } else { 
          df_out <- ggraph(dels_tidy, layout = 'linear')  + geom_vline(xintercept = c(13655647, 13655981), colour="darkgoldenrod1") + geom_edge_arc(strength = 1,  aes(colour=Microhomology), alpha=.5) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(13654306, 13657310)) + scale_x_continuous(breaks=seq(13654306, 13657310, 100)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red")
        }

        if (inclFocus == "yes") {
            df_out_focus_site1 <- ggraph(dels_tidy, layout = 'linear')  + geom_vline(xintercept = c(13655647, 13655981), colour="darkgoldenrod1") + geom_edge_arc(strength = 1,  aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(13655597, 13655697)) + scale_x_continuous(breaks=seq(13655597, 13655697, 5)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red")
            df_out_focus_site2 <- ggraph(dels_tidy, layout = 'linear')  + geom_vline(xintercept = c(13655647, 13655981), colour="darkgoldenrod1") + geom_edge_arc(strength = 1,  aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(13655931, 13656031)) + scale_x_continuous(breaks=seq(13655931, 13656031, 5)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red")
        }
    }
}
#############################  MISEQ  ############################
else if (isNpore == "no") {
    if (targ == "site2_915"){
        if (plotMhomol == "no") {
            df_out <- ggraph(dels_tidy, layout = 'linear')  + geom_vline(xintercept = c(5187987, 5188243), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(5187887, 5188343)) + scale_x_continuous(breaks=seq(5187887, 5188343, 50))    
        } else  if (plotMhomol == "2colours") {
          df_out <- ggraph(dels_tidy, layout = 'linear')  + geom_vline(xintercept = c(5187987, 5188243), colour="darkgoldenrod1") + geom_edge_arc(strength = 1,  aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(5187887, 5188343)) + scale_x_continuous(breaks=seq(5187887, 5188343, 50)) + scale_edge_colour_gradient(limits = c(0,2), low = "blue", high = "blue" , na.value = "red")
        } else {
          df_out <- ggraph(dels_tidy, layout = 'linear')  + geom_vline(xintercept = c(5187987, 5188243), colour="darkgoldenrod1") + geom_edge_arc(strength = 1,  aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(5187887, 5188343)) + scale_x_continuous(breaks=seq(5187887, 5188343, 50)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red")
        }

      #this added to get single nt mapping of very short dels  #cw 08022024
      if (inclFocus == "yes") {
        df_out_focus_site1 <- ggraph(dels_tidy, layout = 'linear')  + geom_vline(xintercept = c(5187987, 5188243), colour="darkgoldenrod1") + geom_edge_arc(strength = 1,   alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(5187974, 5188000)) + scale_x_continuous(breaks=seq(5187974, 5188000, 1)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red")
        df_out_focus_site2 <- ggraph(dels_tidy, layout = 'linear')  + geom_vline(xintercept = c(5187987, 5188243), colour="darkgoldenrod1") + geom_edge_arc(strength = 1,  alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(5188231, 5188255)) + scale_x_continuous(breaks=seq(5188231, 5188255, 1)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red")
      }
    }

    else if (targ == "ercc1_2ab"){
        if (plotMhomol == "no") {
            df_out <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(1479622, 1479796), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(1479521, 1479999) ) + scale_x_continuous(breaks=seq(1479521, 1479999, 50)) 
        } else if (plotMhomol == "2colours") { 
            df_out <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(1479622, 1479796), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(1479521, 1479999) ) + scale_x_continuous(breaks=seq(1479521, 1479999, 50)) + scale_edge_colour_gradient(limits = c(0,2), low = "blue", high = "blue" , na.value = "red")
        } else { 
          df_out <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(1479622, 1479796), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(1479521, 1479999) ) + scale_x_continuous(breaks=seq(1479521, 1479999, 50)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red")
        }

        if (inclFocus == "yes") {
            df_out_focus_site1 <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(1479622, 1479796), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(1479592, 1479652) ) + scale_x_continuous(breaks=seq(1479592, 1479652, 5)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red") 
            df_out_focus_site2 <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(1479622, 1479796), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(1479766, 1479826) ) + scale_x_continuous(breaks=seq(1479766, 1479826, 5)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red") 
        }
    }

    else if (targ == "ercc1_2cd"){
        if (plotMhomol == "no") {
            df_out <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(1479623, 1479778), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(1479521, 1479999) ) + scale_x_continuous(breaks=seq(1479521, 1479999, 50)) 
        } else  if (plotMhomol == "2colours") { 
            df_out <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(1479623, 1479778), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(1479521, 1479999) ) + scale_x_continuous(breaks=seq(1479521, 1479999, 50)) + scale_edge_colour_gradient(limits = c(0,2), low = "blue", high = "blue" , na.value = "red")
        } else { 
          df_out <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(1479623, 1479778), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(1479521, 1479999) ) + scale_x_continuous(breaks=seq(1479521, 1479999, 50)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red")
        }

        if (inclFocus == "yes") {
            df_out_focus_site1 <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(1479623, 1479778), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(1479593, 1479653) ) + scale_x_continuous(breaks=seq(1479593, 1479653, 5)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red") 
            df_out_focus_site2 <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(1479623, 1479778), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(1479748, 1479808) ) + scale_x_continuous(breaks=seq(1479748, 1479808, 5)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red") 
        }
    }

    else if (targ == "xpf_2ab"){
        if (plotMhomol == "no") {
            df_out <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(16476847, 16476711), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(16475300, 16478300) ) + scale_x_continuous(breaks=seq(16475300, 16478300, 50))
        } else  if (plotMhomol == "2colours") { 
            df_out <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(16476847, 16476711), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(16475300, 16478300) ) + scale_x_continuous(breaks=seq(16475300, 16478300, 50)) + scale_edge_colour_gradient(limits = c(0,2), low = "blue", high = "blue" , na.value = "red")
        } else { 
          df_out <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(16476847, 16476711), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(16475300, 16478300) ) + scale_x_continuous(breaks=seq(16475300, 16478300, 50)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red")
        }

        if (inclFocus == "yes") {
            df_out_focus_site1 <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(16476847, 16476711), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(16476817, 16476877) ) + scale_x_continuous(breaks=seq(16476817, 16476877, 5)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red") 
            df_out_focus_site2 <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(16476847, 16476711), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(16476681, 16476741) ) + scale_x_continuous(breaks=seq(16476681, 16476741, 5)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red") 
        }
    }

    else if (targ == "send1_2cd"){
        if (plotMhomol == "no") {
            df_out <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(18131893, 18132057), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(18131752, 18132152) ) + scale_x_continuous(breaks=seq(18131752, 18132152, 50)) 
        } else if (plotMhomol == "2colours") { 
            df_out <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(18131893, 18132057), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(18131752, 18132152) ) + scale_x_continuous(breaks=seq(18131752, 18132152, 50)) + scale_edge_colour_gradient(limits = c(0,2), low = "blue", high = "blue" , na.value = "red")
        } else { 
          df_out <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(18131893, 18132057), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(18131752, 18132152) ) + scale_x_continuous(breaks=seq(18131752, 18132152, 50)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red")
        }

        if (inclFocus == "yes") {
            df_out_focus_site1 <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(18131893, 18132057), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(18131863, 18131923) ) + scale_x_continuous(breaks=seq(18131863, 18131923, 5)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red") 
            df_out_focus_site2 <- ggraph(dels_tidy, layout = 'linear') + geom_vline(xintercept = c(18131893, 18132057), colour="darkgoldenrod1") + geom_edge_arc(strength = 1, aes(colour=Microhomology), alpha=.2) + deletions_theme + labs(x = "Coordinates (bp)", subtitle = "") + coord_cartesian(xlim = c(18132027, 18132087) ) + scale_x_continuous(breaks=seq(18132027, 18132087, 5)) + scale_edge_colour_gradient(limits = c(0,5), low = "blue", high = "orange" , na.value = "red") 
        }
    }

    else if (targ == "send1_2ab"){
    }

}
########################################################################
   
 path_out <- sub(pattern = ".txt", replacement = "_full.png", x = path)
 ggsave(path_out, df_out, bg = "transparent")

if (inclFocus == "yes") {
    path_out_focus_site1 <- sub(pattern = ".txt", replacement = "_cs1.png", x = path)
    ggsave(path_out_focus_site1, df_out_focus_site1, bg = "transparent")
    if (numsites == 2) {
        path_out_focus_site2 <- sub(pattern = ".txt", replacement = "_cs2.png", x = path)
        ggsave(path_out_focus_site2, df_out_focus_site2, bg = "transparent")
    } 
}
         
       })


#clean workspace
rm(list = ls())

