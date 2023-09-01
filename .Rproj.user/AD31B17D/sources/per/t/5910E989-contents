library(magrittr)
path_fastq <- "/home/juan/MitoR/Miami/Mouse_Amplicons/00_fastq/Ch9-TALGFP-12W/Ch9-TALGFP-12W_R2_001_NC9-2-B.fastq" %>%
  strsplit("\\.")

# Cargar fastq
fastq <- readLines(sprintf("%s.%s", path_fastq[[1]][1], path_fastq[[1]][2]))

# Dejo solamente las secuencias
indices <- seq(2, length(fastq), by = 4)
sequences <- fastq[indices]

# Secuencia a eliminar: CCGTATAT
fastq_trimeado <- stringr::str_remove(sequences, "TGCTCCGA")

fastq[indices] <- fastq_trimeado

writeLines(fastq, sprintf("%s_trimmed.%s", path_fastq[[1]][1], path_fastq[[1]][2]))

