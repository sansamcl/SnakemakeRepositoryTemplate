#!/usr/bin/env Rscript

# Chris Sansam
# version 02
# April 16, 2022

# load libraries
library(bamsignals)
library(GenomicRanges)
library(Rsamtools)

# set variables passed through snakemake 'input', 'output', and 'configs' lists
bedFile <- snakemake@config[["rfd_Windows"]]
blacklistBedFile <- snakemake@config[["blacklistBedFile"]]
posBamFile <- snakemake@input[[1]]
negBamFile <- snakemake@input[[3]]
rfdBedgraph <- snakemake@output[[1]]
posBedgraph <- snakemake@output[[2]]
negBedgraph <- snakemake@output[[3]]

# set terms for removing chromosomes absent from bam
#negativeGrepTerms <- "chrX|alt|chrY|random|chrM"

# read analysis windows into dataframe and convert to genomic ranges
## this could be done in one step with rtracklayer, but we've had trouble with this package on conda
df <- read.table(bedFile)
gr <- makeGRangesFromDataFrame(df,
                         ignore.strand=T,
                         seqnames.field="V1",
                         start.field="V2",
                         end.field="V3")

# rfd cannot be calculated for chromosomes without reads, 
# so the windows in the df with seqnames not present in both .bam file are removed
sqnames <- unique(c(idxstatsBam(posBamFile)$seqnames,
  idxstatsBam(negBamFile)$seqnames))
gr <- gr[which(as.vector(seqnames(gr) %in% sqnames))]

# count reads on watson strand
watson <- bamCount(negBamFile,
                   gr,
                   ss=FALSE,
                   paired.end="midpoint")

# count reads on crick strand
crick <- bamCount(posBamFile,
                   gr,
                   ss=FALSE,
                   paired.end="midpoint")


# read blacklist windows into dataframe and convert to genomic ranges
## this could be done in one step with rtracklayer, but we've had trouble with this package on conda
bl_df <- read.table(blacklistBedFile)
bl_gr <- makeGRangesFromDataFrame(bl_df,
                         ignore.strand=T,
                         seqnames.field="V1",
                         start.field="V2",
                         end.field="V3")

# count reads on watson strand in blacklisted windows
bl_watson <- bamCount(negBamFile,
                   bl_gr,
                   ss=FALSE,
                   paired.end="midpoint")

# count reads on crick strand in blacklisted windows
bl_crick <- bamCount(posBamFile,
                   bl_gr,
                   ss=FALSE,
                   paired.end="midpoint")

# find windows overlapping blacklist regions
overlaps_df <- as.data.frame(findOverlaps(bl_gr,gr))

# add blacklisted region read counts to overlaps_df
overlaps_df$bl_watson <- bl_watson[overlaps_df$queryHits]
overlaps_df$bl_crick <- bl_crick[overlaps_df$queryHits]

# sum blacklisted reads in each window
lst <- split(overlaps_df,f=df$subjectHits)
lst2 <- lapply(lst,function(df){data.frame(
  "subjectHits" = df$subjectHits[1],
  "bl_watson" = sum(df$bl_watson),
  "bl_crick" = sum(df$bl_crick))})
overlaps_df2 <- do.call(rbind,lst2)
bl_watson2 <- rep(0,length(gr))
bl_watson2[overlaps_df2$subjectHits] <- overlaps_df2$bl_watson
bl_crick2 <- rep(0,length(gr))
bl_crick2[overlaps_df2$subjectHits] <- overlaps_df2$bl_crick

# subtract blacklisted read counts from each window
crick2 <- crick - bl_crick2
watson2 <- watson - bl_watson2

# calculate rfd
rfd <- (crick2-watson2)/(crick2+watson2)

# write bedgraphs
writeBedgraph <- function(vctr,filename=rfdBedgraph){
  score(gr) <- vctr
  gr2 <- gr[which(!is.na(vctr))]
  write.table(data.frame(as.data.frame(gr2)[,c(1:3)],score(gr2)),
              file=filename,
              quote= F,
              row.names = F,
              col.names = F,
              sep="\t")
}

writeBedgraph(rfd,rfdBedgraph)
writeBedgraph(crick2,posBedgraph)
writeBedgraph(watson2,negBedgraph)
