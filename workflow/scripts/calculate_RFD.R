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


rfd <- (crick-watson)/(crick+watson)


writeBedgraph <- function(vctr,filename=rfdBedgraph){
  score(gr) <- vctr
  gr2 <- gr[which(is.numeric(vctr))]
  write.table(data.frame(as.data.frame(gr2)[,c(1:3)],score(gr2)),
              file=filename,
              quote= F,
              row.names = F,
              col.names = F,
              sep="\t")
}

writeBedgraph(rfd,rfdBedgraph)
writeBedgraph(crick,posBedgraph)
writeBedgraph(watson,negBedgraph)
