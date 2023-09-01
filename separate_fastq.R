fastq <- readLines("/home/juan/MitoR/Miami/Mouse_Amplicons/00_fastq/Ch9-TALGFP-12W/Ch9-TALGFP-12W_R2_001.fastq")
indices <- seq(2, length(fastq), by = 4)

sequences <- fastq[indices]

fastq_BC <- c()
# Termino con el indice de las secuencias que tienen el bar-code indicado
for (i in (1:length(sequences))) {
  if (grepl("^GTACAGGA", sequences[i])) {
    fastq_BC <- c(fastq_BC, i)
  }
}

fastq_BC_4 <- fastq_BC * 4
fastq_BC_3 <- fastq_BC_4 - 3
fastq_BC_2 <- fastq_BC_4 - 2
fastq_BC_1 <- fastq_BC_4 - 1

fastq_BC_total <- c(fastq_BC_4, fastq_BC_3, fastq_BC_2, fastq_BC_1)
fastq_BC_sorted <- sort(fastq_BC_total)

fastq_BC <- fastq[fastq_BC_sorted]

writeLines(fastq_BC, "/home/juan/MitoR/Miami/Mouse_Amplicons/00_fastq/Ch9-TALGFP-12W/Ch9-TALGFP-12W_R2_001_NC9-1-B.fastq")
