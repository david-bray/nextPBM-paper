#!/usr/bin/env Rscript

# fully customizable motif scanner - USED FOR PU.1 PBM DESIGN

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
  # return max of the scores (either plus strand or minus strand)
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

get.top.motif.list <- function(DNA_seq, motif){
  
  ## INITIALIZE OBJECTS NEEDED FOR SCANNING AND KEEPING TRACK OF SCORES
  # if log_likehood is True, get the background seq
  # save the motif length for convenience
  motif_length <- length(motif[1,])
  # reverse the motif to test the reverse sequences simultaneously
  motif_rev <- motif[,ncol(motif):1]
  # initialize max motif score
  max_motif_score <- -9999
  # initialize a corresponding best sub-sequence
  max_sub_DNA <- ""
  # initialize an integer to represent the direction and strand of the top motif
  max_seq_dir <- 0
  # initialize integer to hold the position of the top-scoring motif (offset relative to the reference genome location)
  max_offset <- 0
  
  ## SCANS A FULL-SIZE PEAK AND REPORTS THE MAX
  for (i in 1:(nchar(DNA_seq)-motif_length+1)){
    # grab the current DNA substring to score against the motif
    sub_DNA <- substr(DNA_seq, i, i+motif_length-1)
    # NOTE: for motif scores (+, -)
    motif_scores <- get.motif.scores(sub_DNA, motif, motif_rev)
    # test both scores returned against the current maximum
    # update the top-scoring DNA sub-sequence and other parameters if necessary
    for (j in 1:2){
      if (motif_scores[j] > max_motif_score){
        max_motif_score <- motif_scores[j]
        max_sub_DNA <- sub_DNA
        max_seq_dir <- j
        max_offset <- i-1
      }
    }
  }
  
  # transform the top-scoring DNA subsequence (if necessary)
  # recall that max_seq_dir codes for (+, -)
  max_sub_DNA <- get.motif.seq(max_sub_DNA, max_seq_dir)
  
  # return a list with relevant information about the maximum-scoring motif
  top_motif_list <- list("DNA" = max_sub_DNA, "score" = max_motif_score, "seq_dir" = max_seq_dir, "offset" = max_offset)
  return(top_motif_list)
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
    report_row[5:10] <- rep(".", 6)
    
    # return the row to the parent method without executing rest
    return(report_row)
  }
  
  # get the top motif list
  top_motif_list <- get.top.motif.list(DNA_seq, motif)
  
  # construct a row of data (characters) to return to the parent method
  report_row <- vector()
  report_row[1] <- bed_coord[1]
  report_row[2] <- bed_coord[2]
  report_row[3] <- bed_coord[3]
  report_row[4] <- bed_center
  report_row[5] <- bed_start + top_motif_list$offset # should work with UCSC (AFTER FIXING SCRIPT)
  report_row[6] <- bed_start + top_motif_list$offset + length(motif[1,]) - 1 # should work with UCSC
  report_row[7] <- -bed_center + bed_start + top_motif_list$offset
  report_row[8] <- switch(top_motif_list$seq_dir, "+", "-")
  report_row[9] <- top_motif_list$score
  report_row[10] <- top_motif_list$DNA
  
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
# setwd("C:/Users/dbray/Google Drive/Boston University/Spring 2016/Siggers Rotation/THP1_TF_project/analysis/PBM_design/motif_scans/")
# args <- c("PU1_test", "PU1_GSE21512.pwm", "test.bed")

# parse the sequence file (may need to change if there's a header)
seq_bed <- read.table(args[3], sep='\t', header=F, quote="")
seq_bed <- seq_bed[, 1:4]

# load the motif
motif <- as.matrix(read.table(paste(getwd(), "motifs/PWM", args[2], sep='/')))

# if using log likelihood scoring re-construct motif
log_likelihood <- TRUE
if (isTRUE(log_likelihood)){
  motif <- log.transform.motif(motif)
}

# score the sequences using the given motif
seq_scores <- t(apply(seq_bed, 1, score.sequence, ChIP_data=F))
write.table(seq_scores, file=paste(getwd(), "/results/", args[1], "_max_scores.txt", sep=''), quote=F, row.names=F, col.names=F, sep='\t')
