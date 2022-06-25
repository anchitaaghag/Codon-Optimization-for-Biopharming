# Functions to Obtain CDS Information from NCBI & GenBank
# 13 June 2022
# Anchitaa Ghag 

#### 01 Install Packages & Load Libraries ####

# install.packages("rentrez")
library("rentrez")
# install.packages("XML")
library("XML")
#BiocManager::install("genbankr")
library("genbankr")

#### 02 Get_Accession_Number() Function ####

# This function will take a NCBI id and return the corresponding GenBank accession number.
# It uses the "rentrez" package functions and extracts the relevant line of information from a created list.

#Get_Accession_Number <- function(ncbi_id) {
  
  # First, using rentrez package's entrez_summary() function to get summary information for the NCBI id provided.
  #summary_info <- entrez_summary(db="nucleotide", id=ncbi_id)
  # Then, using the created summary_info list, find the accession version number.
  #genbank_accession_number <- summary_info[["accessionversion"]]
  # Finally, return the GenBank accession number.
  #return(genbank_accession_number)
  
  #}

#### 03 Get_CDS_Start() Function ####

# This function will take a GenBank accession number and return the corresponding CDS start site from the range.
# It uses the "genbankr" package functions and extracts the relevant line of information from a created object.

Get_CDS_Ranges <- function(genbank_accession_number) {
  
  # To avoid printing several informative messages to the console while running this function, the suppressMessages() function will be used. This can be removed, if messages are desired, by commenting out this function below.
  # First, using the genbankr package to retrieve GenBank information by the provided versioned accession number.
  # Search for entries that match the search term.
  gba = GBAccession(genbank_accession_number)
  # Read in the GenBank information.
 # suppressMessages(gb_info <- readGenBank(gba, partial=TRUE, verbose = FALSE))
  gb_info <- readGenBank(gba, partial=TRUE, verbose = FALSE)
  # Save the range information to a dataframe.
  dfRange <- as.data.frame(gb_info@cds@ranges)
  # Save the start of the coding sequence.
  range_start <- unlist(dfRange[1,1])
  # Save the end of the coding sequence.
  range_end <- unlist(dfRange[1,2])
  # Add the start and stop positions to a new list.
  ss_list <- as.list(c(range_start,range_end))
  # Finally, return the starting point for the range.
  return(ss_list) 
  
}

#### 04 Get_CDS_End() Function ####

# This function will take a GenBank accession number and return the corresponding CDS end site from the range.
# It uses the "genbankr" package functions and extracts the relevant line of information from a created object.

Get_CDS_End <- function(genbank_accession_number) {
  
  # To avoid printing several informative messages to the console while running this function, the suppressMessages() function will be used. This can be removed, if messages are desired, by commenting out this function below.
  # First, using the genbankr package to retrieve GenBank information by the provided versioned accession number.
  # Search for entries that match the search term.
  gba = GBAccession(genbank_accession_number)
  # Read in the GenBank information.
  suppressMessages(gb_info <- readGenBank(gba, partial=TRUE, verbose = FALSE))
  # Save the range information to a dataframe.
  dfRange <- as.data.frame(gb_info@cds@ranges)
  # Save the end of the coding sequence.
  range_end <- unlist(dfRange[1,2])
  # Finally, return the ending point for the range.
  return(range_end) 
  
}

#### 05 Get_CDS_Length() Function ####

# This function will take the start and stopping points for a sequence and return the corresponding sequence length.

#Get_CDS_Length <- function(start,stop) {
  
  # Subtract the start from the stop to get the overall length of the sequence.
  #range_length <- (stop-start)+1
  # Return the length for the sequence.
  #return(range_length) 
  
#}

#### 06 References ####

#http://bioconductor.riken.jp/packages/3.8/bioc/vignettes/genbankr/inst/doc/genbankr.html
# https://bioconductor.org/packages/release/bioc/html/genbankr.html
# https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html#dealing-with-many-records
