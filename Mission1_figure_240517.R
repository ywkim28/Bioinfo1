## setting directory
setwd("Desktop/bioinfo1_data/")

## setting package
library(dplyr)
library(ggplot2)
library(ggpubr)

## data import
counts <- read.table("read-counts.txt", header = T)

## data qc
# remove columns of 'CLIP.let7g.bam', 'filtered.RPF.siLuc.bam'
counts <- counts[,-c(8,9)]

# remove low level transcript or low level ribosome footprints
counts <- counts %>% filter(RNA.control.bam >= 30 & RNA.siLin28a.bam >= 30 & RNA.siLuc.bam >= 30)
counts <- counts %>% filter(RPF.siLuc.bam >= 80)

# calculation enrichment and change
# log transformation
clip_enrichment.log2 <- log2(counts$CLIP.35L33G.bam / counts$RNA.control.bam)
rden_change.log2     <- log2((counts$RPF.siLin28a.bam / counts$RNA.siLin28a.bam) / (counts$RPF.siLuc.bam / counts$RNA.siLuc.bam))

counts <- cbind(counts, clip_enrichment.log2, rden_change.log2)

# remove infinite value
counts <- counts %>% filter(!is.infinite(clip_enrichment.log2))
counts <- counts %>% filter(!is.infinite(rden_change.log2))

# standardization
counts$clip_enrichment.log2 <- scale(counts$clip_enrichment.log2)
counts$rden_change.log2 <- scale(counts$rden_change.log2)

## correlation coefficient
cor(counts$clip_enrichment.log2, counts$rden_change.log2)

## visualization of Fig.4D
ggplot(counts) +
  geom_point(aes(x = clip_enrichment.log2, y = rden_change.log2),
             size = 0.3,
             colour = "black") +
  xlab("LIN28A CLIP enrichment (log2)") +
  ylab('Ribosome density change upon Lin28a KD (log2)') +
  annotate("text", x = 3, y = -3, label = "r = 0.4789216") +
  ggtitle("(Fig.4D) CLIP and ribosome footprinting upon Lin28a KD")

## data integration
local <- read.csv("mouselocalization-20210507.txt", header = T, sep = "\t")

gene_id <- gsub("[.].*$", "", counts$Geneid)
counts <- cbind(counts, gene_id)

counts <- merge(counts, local, by = "gene_id")

## visualization of Fig.5B
ggplot(counts, aes(x = clip_enrichment.log2, y = rden_change.log2,
                   color = type)) +
  geom_point() +
  xlab("LIN28A CLIP enrichment (log2)") +
  ylab('Ribosome density change upon Lin28a KD (log2)') +
  ggtitle("(Fig.5B) Protein localization")

