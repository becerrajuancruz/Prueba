#' @title Add the patient's analysis to the MitoR DB
#' @description Updates your MitoR DB with the new patient's analysis. Fixes the frequencies values, adds the new mutations and saves the date of the last MitoR DB update.
#' @param filtered_DF Patient's SNP & INDEL report
#' @param filtered_DF2 Patient's software report
#' @param patient_ID Patient's ID
#' @return Vector of frequencies from MitoR DB to add to the patient's report

add_to_RDS <- function(filtered_DF, filtered_DF2, patient_ID) {

  mitor_db <- sprintf("%s/mitorDB/DB", Sys.getenv('R_LIBS_USER'))
  new_patient <- filtered_DF
  filter_params <- filtered_DF2["Filter_Params"] # Filter parameters
  date <- filtered_DF2["Date"] # Date of analysis
  mutations <- paste(new_patient[, "REF"], new_patient[, "POS"], new_patient[, "ALT"], sep = "/") # REF/POS/ALT of every mutation

  # Once the data is load, let's begin
  if (!"RDS_DB.rds" %in% list.files(mitor_db)){
    # Add the frequency colomn with 100% and count = 1 and then sort all the colomns
    new_patient <- cbind(new_patient, Freq_MitoR = "100 %", Count = 1)

    new_patient <- new_patient[c("CHROM", "POS", "REF", "ALT", "Freq_HMTVAR", "Freq_MitoR", "Count", "clinVAR",	"dbSNP",	"MitoMap",	"Omim",	"Disease", "Franklin", "VarSome")]

    # Add the previous variables to the DB
    RDS_DB <- list(Freq_MitoR = new_patient, list(Filter_Params = filter_params, Date = format(Sys.Date(), "%d/%m/%Y"),
                                                  Mutations = mutations))
    names(RDS_DB)[length(RDS_DB)] <- patient_ID
    frequency_of_interest <- RDS_DB$Freq_MitoR$Freq_MitoR

    saveRDS(RDS_DB, sprintf("%s/RDS_DB.rds", mitor_db))

  } else {
    RDS_DB <- readRDS(sprintf("%s/RDS_DB.rds", mitor_db))
    # Adding the new patient as a new element of the list.
    if (!RDS_DB[sprintf("%s", patient_ID)] %in% RDS_DB) {
      # Adds the patient information in the list. The element is named as the name of the patient's ID
      RDS_DB[[length(RDS_DB)+1]] <- list(Filter_Params = filter_params, Date = format(Sys.Date(), "%d/%m/%Y"), Mutations = mutations)
      names(RDS_DB)[length(RDS_DB)] <- patient_ID

      results <- addFreq_MitoR(new_patient)
      Freq_MitoR <- results[[1]]
      frequency_of_interest <- results[[2]]

      RDS_DB$Freq_MitoR <- Freq_MitoR
      saveRDS(RDS_DB, sprintf("%s/RDS_DB.rds", mitor_db))

    } else {     # If it was already there, it removes the previous one to avoid having the same analysis twice
      Freq_MitoR <- deleteFreq_MitoR(patient_ID)
      RDS_DB[[sprintf("%s", patient_ID)]] <- NULL

      # After removing the previous analysis, it adds the patient' again with the patient's new values to the DB
      # Adds the patient information in the list. The element is named as the name of the patient's ID
      RDS_DB[[length(RDS_DB)+1]] <- list(Filter_Params = filter_params, Date = format(Sys.Date(), "%d/%m/%Y"), Mutations = mutations)
      names(RDS_DB)[length(RDS_DB)] <- patient_ID

      results <- addFreq_MitoR(new_patient)
      Freq_MitoR <- results[[1]]
      frequency_of_interest <- results[[2]]

      RDS_DB$Freq_MitoR <- Freq_MitoR
      saveRDS(RDS_DB, sprintf("%s/RDS_DB.rds", mitor_db))
    }
  }
  return(frequency_of_interest) # MitoR DB frequencies
}


#' @title Updates the MitoR DB with the new patient
#' @description Fixes the frequencies and count values by going through the patient's DataFrame report and the previous MitoR DB, comparing each mutation and looking for matches.
#' @param patient_ID Patient's ID
#' @return List with the new MitoR DB and the frequencies of the mutations of the new patients on the MitoR DB

addFreq_MitoR <- function(new_patient) {

  #libPath <- dirname(dirname(dirname(system.file(package = "MitoR"))))
  mitor_db <- sprintf("%s/mitorDB/DB", Sys.getenv('R_LIBS_USER'))
  RDS_DB <- readRDS(sprintf("%s/RDS_DB.rds", mitor_db))
  Freq_MitoR <- as.data.frame(RDS_DB[[1]])
  new_patient <- cbind(new_patient, Count = 1, Freq_MitoR = "100 %")
  new_patient <- new_patient[c("CHROM", "POS", "REF", "ALT", "Freq_HMTVAR", "Freq_MitoR", "Count", "clinVAR", "dbSNP", "MitoMap", "Omim", "Disease",  "Franklin", "VarSome")]
  newAmount <- length(RDS_DB) # Don't substract 1 because the new patient will be saved on the RDS at the previous function

  # newMut is a DF with the mutations that were not found in the DB before loading the new patient.
  newMut <- data.frame()
  # frequency_of_interest will be used as a new coloumn in the particular analysis of the patient.
  frequency_of_interest <- c()
  # Initialize to go through the previous DB and the new patients DF
  i <- 1 # Related to the DB rows (Freq_MitoR)
  j <- 1 # Related to the new patient's DF rows (new_patient)
  Freq_MitoR$POS <- as.integer(Freq_MitoR$POS)
  new_patient$POS <- as.integer(new_patient$POS)

  # It goes through the DB and the new patients's DF simultaneously. To simplify we will call them as
  # table A and table B. It goes through them checking the POS of both tables: First it checks the first
  # row of each table. If table A has a higher POS than table B, it will compare the same row of table A
  # with the next row of table B. This will keep happening until table B has a higher position than
  # table A, or they match their positions. In case table B has a higher position, the same will happen
  # the other way around.
  # This stops once one table has come to its last row.

  while ((i <= nrow(Freq_MitoR)) && (j <= nrow(new_patient))){
    if (Freq_MitoR$POS[i] > new_patient$POS[j]){
      # In this case the DB table has increased its rows until its POS is higher than
      # the new patient's one and the mutation was not found.
      # Meaning the mutation of the patient was not in the DB and it has to be loaded as a new one.
      new_patient$Freq_MitoR[j] <- sprintf("%s%s", 100/(newAmount), "%")
      frequency_of_interest <- c(frequency_of_interest, new_patient$Freq_MitoR[j])

      newMut <- rbind(new_patient[j, ], newMut)
      j <- j + 1

    } else if (Freq_MitoR$POS[i] < new_patient$POS[j]){
      # In this case the DB POS is lower tan the new patient's POS. It increases its row trying
      # to match with the next mutation.
      Freq_MitoR$Freq_MitoR[i] <- sprintf("%s%s", (Freq_MitoR$Count[i]*100)/(newAmount), "%")
      i <- i + 1

    } else if (Freq_MitoR$POS[i] == new_patient$POS[j]){ # Checks only their POS at first
      # When the mutation match, the DB values must be changed and then both DF go the their next row.
      if ((Freq_MitoR$ALT[i] == new_patient$ALT[j]) && (Freq_MitoR$REF[i] == new_patient$REF[j])){ # Chequea el resto
        # Updates the allele frequency just in case it has changed in HMTVAR
        Freq_MitoR$All_freq_h[i] <- new_patient$All_freq_h[j]

        # Increases its count and changes the DB frequency
        Freq_MitoR$Count[i] <- (Freq_MitoR$Count[i] + 1)
        Freq_MitoR$Freq_MitoR[i] <- sprintf("%s%s", (Freq_MitoR$Count[i]*100)/(newAmount), "%")
        frequency_of_interest <- c(frequency_of_interest, Freq_MitoR$Freq_MitoR[i])
        j <- j + 1
        i <- i + 1
      }
    }
  }

  if (i < nrow(Freq_MitoR)){
    # In case it has come to the end of the new patient's DF, it means there are no more mutations
    # to compare. In this case it updates the frequency of the uncompared rows of the previous DB.
    Freq_MitoR$Freq_MitoR[i:nrow(Freq_MitoR)] <- sprintf("%s%s", (Freq_MitoR$Count[i:nrow(Freq_MitoR)]*100)/newAmount, "%")
  }

  # On the other hand, if the DB is the one that has come to the end first, all the resting mutations
  # of the new patients table are new to the DB and must be saved in the DB.
  if (j < nrow(new_patient)) {
    new_patient$Freq_MitoR[j:nrow(new_patient)] <- sprintf("%s%s", 100/(newAmount), "%")
    frequency_of_interest <- c(frequency_of_interest, new_patient$Freq_MitoR[j:nrow(new_patient)])

    newMut <- rbind(new_patient[j:nrow(new_patient), ], newMut)
  }

  # Makes sure there are new mutations. If not, just avoid the computational process of adding nothing
  # to a DF
  if (length(newMut) != 0){
    colnames(newMut) <- c("CHROM", "POS", "REF", "ALT", "Freq_HMTVAR", "Freq_MitoR", "Count", "clinVAR", "dbSNP", "MitoMap", "Omim", "Disease", "Franklin", "VarSome")
    Freq_MitoR <- rbind(Freq_MitoR, newMut)
  }

  Freq_MitoR <- Freq_MitoR[order(Freq_MitoR$POS), ]
  Freq_MitoR$Freq_HMTVAR[Freq_MitoR$Freq_HMTVAR == "-"] <- "0"
  return(list(Freq_MitoR, frequency_of_interest))
}

#' @title Deletes a patient from the MitoR DB
#' @description Fixes the frequencies and count values by going through the patient's mutations element of the RDS list. Mutations with count = 0 are deleted.
#' @param patient_ID Patient's ID
#' @return Updated MitoR DB with fixed values of frequencies and counts

deleteFreq_MitoR <- function(patient_ID) {
  mitor_db <- sprintf("%s/mitorDB/DB", Sys.getenv('R_LIBS_USER'))
  RDS_DB <- readRDS(sprintf("%s/RDS_DB.rds", mitor_db))
  Freq_MitoR <- data.frame(RDS_DB[[1]])
  # Split the mutation coloumn in ALT POS REF
  toDelete <- unlist(RDS_DB[[patient_ID]]["Mutations"], use.names = FALSE) %>%
    strsplit("/")
  newAmount <- length(RDS_DB) - 2

  for (i in (1:length(toDelete))){
    Freq_MitoR[(Freq_MitoR$POS == as.integer(toDelete[[i]][2])) & (Freq_MitoR$ALT == toDelete[[i]][3]), ]$Count <- (Freq_MitoR[(Freq_MitoR$POS == as.integer(toDelete[[i]][2])) & (Freq_MitoR$ALT == toDelete[[i]][3]), ]$Count - 1)

    # In case we delete a mutation with Count = 1, just delete the entire row
    if (Freq_MitoR[(Freq_MitoR$POS == as.integer(toDelete[[i]][2])) & (Freq_MitoR$ALT == toDelete[[i]][3]), ]$Count == 0) {
      Freq_MitoR <- Freq_MitoR[!((Freq_MitoR$POS == as.integer(toDelete[[i]][2])) & (Freq_MitoR$ALT == toDelete[[i]][3])), ]
    }
  }
  # Updates the frequencies of the DB
  Freq_MitoR$Freq_MitoR <- sprintf("%s%s", (Freq_MitoR$Count * 100)/newAmount, "%")

  return(Freq_MitoR)
}

