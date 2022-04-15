#!/usr/bin/env Rscript

# Chris Sansam
# version 01
# April 14, 2022

bedFile <- snakemake@config[["rfd_Windows"]]
negativeGrepTerms <- "chrX|alt|chrY|random|chrM"
bamFile <- snakemake@input[[1]]
rfdBedgraph <- snakemake@output[[1]]
posBedgraph <- snakemake@output[[2]]
negBedgraph <- snakemake@output[[3]]
library(bamsignals)
library(GenomicRanges)
df <- read.table(bedFile)
df2 <- df[-grep(negativeGrepTerms,df$V1),]
gr <- makeGRangesFromDataFrame(df2,
                         ignore.strand=T,
                         seqnames.field="V1",
                         start.field="V2",
                         end.field="V3")
test <- bamCount(bamFile,
                   gr,
                   mapqual=30,
                   ss=TRUE,
                   paired.end="midpoint",
                   tlenFilter=c(10,1000))
rfd <- (test[1,]-test[2,])/(test[1,]-test[2,])

writeBedgraph <- function(vctr,filename=rfdBedgraph){
  score(gr) <- vctr
  write.table(data.frame(as.data.frame(gr)[,c(1:3)],score(gr)),
              file=filename,
              quote= F,
              row.names = F,
              col.names = F,
              sep="\t")
}

writeBedgraph(rfd,rfdBedgraph)
writeBedgraph(test[1,],posBedgraph)
writeBedgraph(test[2,],negBedgraph)
