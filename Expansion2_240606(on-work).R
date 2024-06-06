###### pre-processing
#### setting repository
setwd("/disk2/bikyw/7.bioinfo/")

#### fastqc (before cutadapt)
system(paste("fastqc", "--thread", "8", "./00.raw/*.gz", "--outdir", "./99.fastqc/"))

#### trim reads (using cutadapt)
adapter_seq <- "AGATCGGAAGAG"

forward <- paste("-a", adapter_seq, "-o", "./01.trimmed/trimmed_read1.fq.gz")
reverse <- paste("-A", adapter_seq, "-p", "./01.trimmed/trimmed_read2.fq.gz")

system(paste("cutadapt", "--cores", "8", forward, reverse, "./00.raw/DDX52_clip_read1.fastq.gz", "./00.raw/DDX52_clip_read2.fastq.gz"))

#### fastqc (after cutadapt)
system(paste("fastqc", "--thread", "8", "./01.trimmed/*.gz", "--outdir", "./99.fastqc/"))

#### indexing ref genome (using hisat2)
indexing <- "/program/HISAT2/hisat2-build"

system(paste("gzip", "-d", "./98.ref/*.gz"))

system(paste(indexing, "-p", "8", "./98.ref/*.fa", paste0("./98.ref/", "indexed")))

#### mapping, converting and sorting (using hisat2 and samtools)
mapping   <- "/program/HISAT2/hisat2"

ref_idx   <- paste("-x", "./98.ref/indexed")
forward   <- paste("-1", "./01.trimmed/trimmed_read1.fq.gz")
reverse   <- paste("-2", "./01.trimmed/trimmed_read2.fq.gz")

mapping_cmd <- paste(mapping, "-p", "8", ref_idx, forward, reverse)
convert_cmd <- paste("samtools view -Sb --threads 8")
sorting_cmd <- paste("samtools sort --threads 8")

out         <- "./02.bam/DDX52_clip.bam"

system(paste(mapping_cmd, "|", convert_cmd, "|", sorting_cmd, "-o", out))


#### mapping rate
system(paste("samtools", "flagstat", "./02.bam/DDX52_clip.bam", ">", "./02.bam/flagstat.txt"))

#### find DDX52 information from gene annotation
system(paste("gzip", "-d", "./98.ref/*.gz"))

system(paste("grep", "-i", "DDX52", "./98.ref/Homo_sapiens.GRCh38.112.gtf"))

## DDX52 genome annotation information
# chr : 17
# range : (37609739..37643446)

#### bam indexing
system(paste("samtools", "index", "./02.bam/DDX52_clip.bam", "./02.bam/DDX52_clip.bai"))

#### filter location of DDX52 
system(paste("samtools", "view", "-b", "-o", "./02.bam/clip_DDX52.bam", "./02.bam/DDX52_clip.bam", "17:37609000-37644000"))


