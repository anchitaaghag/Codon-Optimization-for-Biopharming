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

#### 01 Get_Accession_Number() Function ####

# This function will take a NCBI id and return the corresponding GenBank accession number.
# It uses the "rentrez" package functions and extracts the relevant line of information from a created list.

Get_Accession_Number <- function(ncbi_id) {
  
  # First, using rentrez package's entrez_summary() function to get summary information for the NCBI id provided.
  summary_info <- entrez_summary(db="nucleotide", id=ncbi_id)
  # Then, using the created summary_info list, find the accession version number.
  genbank_accession_number <- summary_info[["accessionversion"]]
  # Finally, return the GenBank accession number.
  return(genbank_accession_number)
  
  }

#### 02 Get_CDS_Range_Start() Function ####

# This function will take a GenBank accession number and return the corresponding CDS start site from the range.
# It uses the "genbankr" package functions and extracts the relevant line of information from a created object.

Get_CDS_Range_Start <- function(genbank_accession_number) {
  
  # First, using the genbankr package to retrieve GenBank information by the provided versioned accession number.
  # Search for entries that match the search term.
  gba = GBAccession(genbank_accession_number)
  # Read in the genbank information.
  test <- readGenBank(gba, partial=TRUE)
  # Save the range information to a dataframe.
  testing <- as.data.frame(test@cds@ranges)
  # Save the start, end, and length of the CDS sequences to corresponsding lists.
  range_start <-  testing[1]
  # Finally, return the starting point for the range.
  return(range_start) 
  
  }


