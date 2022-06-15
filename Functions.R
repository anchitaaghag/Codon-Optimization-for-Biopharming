# FUNCTIONS FOR 

# 5 FUNCTIONS 

# install.packages("rentrez")
library("rentrez")
# install.packages("XML")
library("XML")
#BiocManager::install("genbankr")
library("genbankr")

# Input should be a NCBI ID

NCBI <- IDs
Accession <- list()
Start <- list()
End <- list()
Length <- list()

#### Get_Accession_Number() Function ####

# This function will take a NCBI id and return the corresponding GenBank accession number.
# It uses the "rentrez" pacakge functionality 

Get_Accession_Number <- function(ncbi_id) {
  
  # First, using rentrez package's entrez_summary() function to get summary information for the NCBI id provided.
  summary_info <- entrez_summary(db="nucleotide", id=ncbi_id)
  # Then, using the created summary_info list, find the accession version number.
  genbank_accession_number <- summary_info[["accessionversion"]]
  # Finally, return the GenBank accession number.
  return(genbank_accession_number)
  
   }

#### Get_CDS_Range_Start() Function ####

Get_CDS_Range_Start <- function(genbank_accession_number) {
  
  # Next, using the genbank package to retrieve GenBank information by the versioned accession found above.
  # Search for entries that match the search term.
  gba = GBAccession(genbank_accession_number)
  # Read in the genbank information.
  test <- readGenBank(gba, partial=TRUE)
  # Save the range information to a dataframe.
  testing <- as.data.frame(test@cds@ranges)
  # Save the start, end, and length of the CDS sequences to corresponsding lists.
  Range_Start <-  testing[1]
  Range_End <- testing[2]
  Range_Length <- testing[3]
  
  return(Range_Start)
  
}


