library(openxlsx)
library(bioseq)
library(magrittr)

# Carga variantes del paciente para cambiarselas al genoma de referencia
path <- "/home/juan/Github/MitochondriaAnalysis/Pacientes/47286/MitoR/47286_filtered_MitoR.xlsx"
mutations <- openxlsx::read.xlsx(as.character(path), sheet = 1)[, c(2, 3, 4)]

# Primero hay que cargar el string de https://www.ncbi.nlm.nih.gov/nuccore/J01415.2?report=fasta
mito_reference <- strsplit(mito_reference, "") %>% unlist()
mito_reference <- mito_reference[mito_reference != "\n"]

# Cambio en el paciente los nucleotidos segun el analisis de variantes
paciente_47286 <- mito_reference
for (i in 1:nrow(mutations)) {
  paciente_47286[as.integer(mutations[i, "POS"])] <- mutations[i, "ALT"]
}

mito_reference <- paste(mito_reference, sep = "", collapse = "")
paciente_47286 <- paste(paciente_47286, sep = "", collapse = "")

# Pasa a ser objeto DNA
mito_reference_DNA <- bioseq::as_dna(mito_reference)
paciente_47286_DNA <- bioseq::as_dna(paciente_47286)

# Transcripcion, traduccion y vuelta a objeto string
mito_reference_AA <- bioseq::seq_transcribe(mito_reference_DNA) %>%
  bioseq::seq_translate(, code = 2) %>% as.character()
paciente_47286_AA <- bioseq::seq_transcribe(mito_reference_DNA) %>%
  bioseq::seq_translate(, code = 2) %>% as.character()


# Guarda como .txt
save(mito_reference_AA, file = "/home/juan/Juan/Mitocondria_Ref_AA.txt")
save(paciente_47286_AA, file = "/home/juan/Juan/47286_AA.txt")

# Para ver los genes implicados, primero ajusto el intervalo del gen segun inserciones y deleciones
# Load BED file with: Gene name, start position, end position
gene_start_end <- readRDS("/home/juan/Github/Proteinas/bedfileMito.RDS")

# Dealing with insertions and deletions on the length of the gene
gene_start_end <- indels_in_gene(mutations, gene_start_end)

# - Si es SNP, analizar si el codon es de STOP. Largar flag y sumar gen posterior tambien
# - Si es INDEL, quilombo:
#    + Buscar el STOP mas cercano
#    + Ver como queda tambien el proximo gen


#' @title Modifies the BED table. Makes a gene longer or shorter depending on the patients variants.
#' @description Depending on the position of the insertion or deletion of the variant, it changes the length of the gene's start and end
#' @param mutations DataFrame of colomns REF-POS-ALT
#' @param start_end_gene BED file - sprintf("%s/MitoRSoftware/bedfileMito.RDS", mitor_files)
#' @return Updated BED table as a DataFrame object
indels_in_gene <- function(mutations, start_end_gene) {
  j <- 1
  for (i in (1:nrow(mutations))) {
    indel <- nchar(mutations[i, "ALT"]) - nchar(mutations[i, "REF"])
    if (indel == 0) {
      # Continues
    } else {
      # Insertion or Deletion
      while(!((mutations[i, "POS"] > start_end_gene[j, "start"]) && (mutations[i, "POS"] < start_end_gene[j, "end"]))){
        j <- j + 1
      }
      start_end_gene[j, "end"] <- start_end_gene[j, "end"] + indel
      if (j < nrow(start_end_gene)){
        start_end_gene[j+1, "start"] <- start_end_gene[j+1, "start"] + indel
      }
    }
  }
  return(start_end_gene)
}

BED_file <- readRDS("/home/juan/Github/Proteinas/bedfileMito.RDS")
polipeptidos <- BED_file[BED_file$name %in% c("ATP8", "ATP6", "CO1", "CO2", "CO3", "CYB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6") ,]

DNA_Seq <- c()
AA_Seq <- c()
for (i in (1:nrow(polipeptidos))) {
  DNA_to_Add <- substr(mito_reference, start = polipeptidos[i, ]$start, stop = polipeptidos[i, ]$end)
  DNA_Seq <- c(DNA_Seq, DNA_to_Add)

  DNA_to_Add <- bioseq::as_dna(DNA_to_Add)
  AA_to_Add <- bioseq::seq_transcribe(DNA_to_Add) %>%
    bioseq::seq_translate(, code = 2) %>% as.character()

  AA_Seq <- c(AA_Seq, AA_to_Add)
}

polipeptidos <- cbind(polipeptidos, DNA_Seq = DNA_Seq, AA_Seq = AA_Seq)

