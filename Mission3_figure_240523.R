## setting directory
setwd("/Users/ywkim/Desktop/bioinfo1_data/")

## load library
library(DescTools)
library(stringr)
library(dplyr)

## data import
data <- read.table("CLIP-let7g-gene.pileup")
colnames(data) <- c("chrom", "pos", "_ref", "count", "basereads", "quals")

## filtering match and substitution
matches <- gsub("[<>$*#^]", "", data$basereads)
data <- cbind(data, matches)

## filtering data
data <- data %>% filter(matches != "")

## calculate shannon entropy
entropy <- c()

for (i in data$matches){
  tmp <- substring(i, 1:nchar(i), 1:nchar(i))
  
  entropy <- c(entropy, Entropy(table(tmp)/nchar(i)))
}

## make bedGraph dataframe
bedGraph <- data.frame(chrom      = data$chrom,
                       chromStart = data$pos,
                       chromEnd   = data$pos+1,
                       dataValue  = entropy)

write.table(bedGraph, "Mission3.bedGraph", col.names = F, row.names = F, quote = F)
