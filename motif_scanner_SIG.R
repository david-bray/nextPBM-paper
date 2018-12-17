#!/usr/bin/env Rscript

# fully customizable motif scanner
# USAGE FOR THIS CUSTOM VERSION ONLY: Rscript --vanilla motif_scanner_CUSTOM_scores_only.R <output prefix> <motif file> <bed file> <seq dir>

## PREVENT SCIENTIFIC NOTATION DURING SESSION
options(scipen=999)

## HELPER METHODS

get.motif.scores <- function(DNA_seq, motif, motif_rev){
  # scores a motif against a sequence window of the same length
  # simultaneous scores the (+) strand and (-) strand and reports the max
  pf_motif_score <- 0
  mr_motif_score <- 0
  # scores each character and adds to the running totals
  for (i in 1:nchar(DNA_seq)){
    # initialize score for this instance
    pf_score <- 0
    mr_score <- 0
    # grab the nucleotide of interest
    nucleotide <- substring(DNA_seq, i, i)
    # grab the corresponding score for the nucleotide for both strands (plus and minus), and both direction (forward and reverse)
    if (nucleotide=='A'){
      pf_score <- motif[1,i]
      mr_score <- motif_rev[4,i]
    }
    if (nucleotide=='C'){
      pf_score <- motif[2,i]
      mr_score <- motif_rev[3,i]
    }
    if (nucleotide=='G'){
      pf_score <- motif[3,i]
      mr_score <- motif_rev[2,i]
    }
    if (nucleotide=='T'){
      pf_score <- motif[4,i]
      mr_score <- motif_rev[1,i]
    }
    # update running scores for each strand and direction
    pf_motif_score <- pf_motif_score + pf_score
    mr_motif_score <- mr_motif_score + mr_score
  }
  # return both scores
  return(c(pf_motif_score, mr_motif_score))
}
# DONE and TESTED

get.rev.seq <- function(DNA_seq){
  # returns the input DNA sequence in reverse
  temp_seq <- unlist(strsplit(DNA_seq, split=''))
  temp_seq <- rev(temp_seq)
  temp_seq <- paste(temp_seq, collapse='')
  return(temp_seq)
}
# DONE and TESTED

get.motif.seq <- function(DNA_seq, DNA_dir){
  # transforms the input sequence if necessary according to the following code:
  # DNA_dir : (+, -)
  if (DNA_dir==2){
    # - strand and reverse direction
    DNA_seq <- get.rev.seq(DNA_seq)
    DNA_seq <- chartr('ACGTN', 'TGCAN', DNA_seq)
  }
  return(DNA_seq)
}
# DONE and TESTED

get.sig.motif.list <- function(DNA_seq, motif, bed_start, bed_center){
  
  ## INITIALIZE OBJECTS NEEDED FOR SCANNING AND KEEPING TRACK OF SCORES
  # save the motif length for convenience
  motif_length <- length(motif[1,])
  # reverse the motif to test the reverse sequences simultaneously
  motif_rev <- motif[,ncol(motif):1]
  # initialize number of significant motifs
  n_sig_motifs <- 0
  # initialize an empty info string to paste significant motif info to
  motif_info_string <- ""
  
  ## SCANS A FULL-SIZE PEAK AND REPORTS THE MAX
  
  # for each substring of DNA, report significant motifs and info as they arise
  for (i in 1:(nchar(DNA_seq)-motif_length+1)){
    # grab the current DNA substring to score against the motif
    sub_DNA <- substr(DNA_seq, i, i+motif_length-1)
    # NOTE: for motif scores (+, -)
    motif_scores <- get.motif.scores(sub_DNA, motif, motif_rev)
    # test both sequences (+ and -) against the significance threshold
    for (j in 1:2){
      if (motif_scores[j] > motif_thresh){ # treshold specified by .txt file that accompanies a .pwm file
        n_sig_motifs <- n_sig_motifs + 1 # update the number of significant motifs for a given peak
        sig_sub_DNA <- sub_DNA
        sig_seq_dir <- j
        sig_offset <- i-1
        
        # to report back:
        sig_motif_start <- bed_start + sig_offset # should work with UCSC                             MOTIF START
        sig_motif_end <- bed_start + sig_offset + motif_length - 1 # should work with UCSC            MOTIF END
        sig_motif_offset <- -bed_center + bed_start + sig_offset #                                    MOTIF OFFSET
        sig_motif_strand <- switch(sig_seq_dir, "+", "-") #                                           MOTIF STRAND
        sig_motif_score <- motif_scores[j] #                                                          MOTIF SCORE
        # transform the top-scoring DNA subsequence (if necessary)
        # recall that sig_seq_dir codes for (+, -)
        sig_sub_DNA <- get.motif.seq(sig_sub_DNA, sig_seq_dir) #                                      MOTIF SEQUENCE
        
        # append all of the info that should be included in the motif_info_string
        motif_info_string <- paste(motif_info_string, " ", as.character(sig_motif_start), ",", as.character(sig_motif_end), ",", as.character(sig_motif_offset), ",", as.character(sig_motif_strand), ",", as.character(sig_motif_score), ",", as.character(sig_sub_DNA), sep="")
      }
    }
  }
  
  # return a list with relevant information about the significant motifs
  sig_motif_list <- list("number" = n_sig_motifs, "info" = motif_info_string)
  return(sig_motif_list)
}
# DONE and TESTED

score.sequence <- function(bed_coord, ChIP_data=T){
  
  ## IF SCORING CHIP DATA
  if(isTRUE(ChIP_data)){
    # get the center of the bed coordinates
    bed_center <- round((as.numeric(bed_coord[2]) + as.numeric(bed_coord[3]))/2)
    
    # use the center and the flank to build new sequence coordinates
    bed_start <- bed_center - flank
    bed_end <- bed_center + flank
  
  ## IF SCORING BACKGROUND DATA
  } else {
    # bed center is not used in the motif finding but still reported (to be consistent and provide offsets)
    bed_center <- round((as.numeric(bed_coord[2]) + as.numeric(bed_coord[3]))/2)
    bed_start <- as.numeric(bed_coord[2])
    bed_end <- as.numeric(bed_coord[3])
  }
  
  # retrieve the given DNA sequence
  DNA_seq <- as.character(bed_coord[4])
  
  # handles the rare case where the DNA sequence being tested is smaller than the motif being used
  if (nchar(DNA_seq) < ncol(motif)) {
    report_row <- vector()
    report_row[1] <- bed_coord[1]
    report_row[2] <- bed_coord[2]
    report_row[3] <- bed_coord[3]
    report_row[4] <- bed_center
    report_row[5] <- 0
    report_row[6] <- ""
    
    # return the row to the parent method without executing rest
    return(report_row)
  }
  
  # EDIT THIS - statement itself is probably fine but the method called needs to be changed
  sig_motif_list <- get.sig.motif.list(DNA_seq, motif, bed_start, bed_center)
  
  # fields to report
  report_row <- vector()
  report_row[1] <- bed_coord[1]
  report_row[2] <- bed_coord[2]
  report_row[3] <- bed_coord[3]
  report_row[4] <- bed_center
  report_row[5] <- sig_motif_list$number # could be reported by the motif list above
  report_row[6] <- sig_motif_list$info # comma-separated fields with space separated entries
  
  # return the row to the parent method (interpreted as column but transposed by the parent method)
  return(report_row)
}
# DONE and TESTED

log.transform.motif <- function(motif){
  # transforms a motif for use in log-likelihood scoring
  motif <- motif/0.25
  motif <- log(motif)
  return(motif)
}
# DONE and TESTED



## ACTUAL SCRIPT

# parse the command line arguments
args <- commandArgs(trailingOnly=TRUE)

# ARGUMENTS:
#     args[1]     output file prefix
#     args[2]     motif file name
#     args[3]     path to peak file (bed format) - WITH SEQUENCES

# # TEST ARGUMENTS
# setwd("C:/Users/dbray/Google Drive/Boston University/Spring 2016/Siggers Rotation/THP1_TF_project/analysis/PBM_design/motif_scans")
# args <- c("PU1_GSE21512", "PU1_GSE21512.pwm", "test.bed")

# parse the sequence file
seq_bed <- read.table(args[3], header=F, quote="", sep='\t')
seq_bed <- seq_bed[, 1:4]

# load the motif and the scoring threshold
motif <- as.matrix(read.table(paste(getwd(), "motifs/PWM", args[2], sep='/')))
motif_thresh <- as.numeric(read.table(paste(getwd(), "/motifs/PWM/", args[1], "_thresh.txt", sep='')))

# if using log likelihood scoring re-construct motif
log_likelihood <- TRUE
if (isTRUE(log_likelihood)){
  motif <- log.transform.motif(motif)
}

# score the sequences using the given motif
seq_scores <- t(apply(seq_bed, 1, score.sequence, ChIP_data=F))
write.table(seq_scores, file=paste(getwd(), "/results/", args[1], "_sig_scores.txt", sep=''), quote=F, row.names=F, col.names=F, sep='\t')
