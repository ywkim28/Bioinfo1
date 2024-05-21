## setting directory
setwd("Desktop/bioinfo1_data/")

## setting package
library(dplyr)
library(ggplot2)

## data handling
data.raw <- read.table("fivepcounts-filtered-RPF-siLuc.txt")

relative_position <- data.raw$V2 - data.raw$V9
counts            <- data.raw$V4
data.filtered <- data.frame(relative_position, counts)

data.count        <- data.filtered %>% group_by(relative_position) %>% summarise(counts = sum(counts))
data.count$counts <- round(data.count$counts/1000,0)

## visualization Figure S5A
ggplot(data.count) +
  geom_bar(aes(relative_position, counts), stat = "identity") +
  xlim(-50, 50) +
  xlab("relative position to start codon of 5'-end of reads") +
  ylab("siLuc raw read count (x1000)") +
  geom_vline(xintercept = 0, col = "red") +
  annotate("text", x = 0, y = 60, label = "start codon")
