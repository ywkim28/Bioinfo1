################################################################################
##### bioinfo1 final presentation
################################################################################

##### setting R environment


##### setting repository
setwd("/disk2/bikyw/8.bioinfo1/")

#### fastqc (before cutadapt)
system(paste("fastqc", "--thread", "8", "./*.fastq.gz"))

#### trim reads
adapter_seq <- "AGATCGGAAGAG"

forward <- paste("-a", adapter_seq, "-o", "./DDX3X_1_trimmed.fq.gz")
reverse <- paste("-A", adapter_seq, "-p", "./DDX3X_2_trimmed.fq.gz")

system(paste("cutadapt", "--cores", "8", forward, reverse, "./DDX3X_1.fastq.gz", "./DDX3X_2.fastq.gz"))

#### fastqc (after cutadapt)
system(paste("fastqc", "--thread", "8", "./*.fq.gz"))

#### indexing ref genome (using hisat2)
indexing <- "/program/HISAT2/hisat2-build"

# system(paste("gzip", "-d", "./*.fa.gz"))

system(paste(indexing, "-p", "32", "./*.fa", "./indexed"))

#### mapping, converting and sorting (using hisat2 and samtools)
mapping   <- "/program/HISAT2/hisat2"

ref_idx   <- paste("-x", "./indexed")
forward   <- paste("-1", "./DDX3X_1_trimmed.fq.gz")
reverse   <- paste("-2", "./DDX3X_2_trimmed.fq.gz")

mapping_cmd <- paste(mapping, "-p", "16", ref_idx, forward, reverse)
convert_cmd <- paste("samtools view -Sb --threads 16")
sorting_cmd <- paste("samtools sort --threads 16")

out         <- "./DDX3X_clip.bam"

system(paste(mapping_cmd, "|", convert_cmd, "|", sorting_cmd, "-o", out))


#### mapping rate
system(paste("samtools", "flagstat", "./DDX3X_clip.bam", ">", "./flagstat.txt"))

#### find DDX3X information from gene annotation
system(paste("gzip", "-d", "./*.gtf.gz"))

system(paste("grep", "-i", "DDX3X", "./Homo_sapiens.GRCh38.112.gtf"))

## DDX3X genome annotation information
# chr : X
# range : 41333348	41364472

#### bam indexing
system(paste("samtools", "index", "./DDX3X_clip.bam", "./DDX3X_clip.bai"))

#### filter location of DDX3X
system(paste("samtools", "view", "-b", "-o", "./filtered.bam", "./DDX3X_clip.bam", "X:41333348-41364472"))
system(paste("samtools", "index", "./filtered.bam", "./filtered.bai"))

#### pileup
system(paste("samtools", "mpileup", "./filtered.bam", ">", "./filtered.pileup"))

system(paste("awk", "'$2 >= 41333348 && $2 <= 41364472 { print $0; }'", "./filtered.pileup",">", "./filtered_gene.pileup"))

#### data cleansing (filter with match and substitution)
library(data.table)
library(dplyr)

data <- fread("./filtered_gene.pileup")
colnames(data) <- c("chrom", "pos", "_ref", "count", "basereads", "quals")
data$relative_position <- seq(1:nrow(data))
data <- data %>% filter(count >= 3)



data$match <- gsub("[<>$*#^\"]", "",data$basereads)
data$match <- toupper(data$match)
data$match <- gsub("]", "", data$match)

data_qc <- data %>% filter(match != "")



## calculate shannon entropy
library(DescTools)

entropy <- c()

for (i in data_qc$match){
  tmp <- substring(i, 1:nchar(i), 1:nchar(i))

  entropy <- c(entropy, Entropy(table(tmp)/nchar(i)))
}


data_qc$entropy <- entropy

barplot(entropy)

## make bedGraph dataframe
bedGraph <- data.frame(chrom      = data_qc$chrom,
                       chromStart = data_qc$pos-1,
                       chromEnd   = data_qc$pos,
                       dataValue  = entropy)

write.table(bedGraph, "DDX3X.bedGraph", col.names = F, row.names = F, quote = F)


