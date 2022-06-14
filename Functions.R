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


Get_Accession_Number <- function(ncbi_id) {
  
  # First, using rentrex package to get summary information for the NCBI id provided.
  summary_info <- entrez_summary(db="nucleotide", id=ncbi_id)
  genbank_accession_number <- summary_info[["accessionversion"]]
  # Save the corresponding GenBank Accession Numbers to a list.
  #Accession <- append(Accession, genbank_accession_number)
  
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
  #Start <- append(Start, Range_Start)
  #End <- append(End, Range_End)
  #Length <- append(Length, Range_Length)
  
  return(genbank_accession_number)
  
    }


Get_Accession_Number <- function(ncbi_id) {
  
  # First, using rentrex package to get summary information for the NCBI id provided.
  summary_info <- entrez_summary(db="nucleotide", id=ncbi_id)
  genbank_accession_number <- summary_info[["accessionversion"]]
  # Save the corresponding GenBank Accession Numbers to a list.
  #Accession <- append(Accession, genbank_accession_number)
  
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
  #Start <- append(Start, Range_Start)
  #End <- append(End, Range_End)
  #Length <- append(Length, Range_Length)
  
  return(genbank_accession_number)
  
}

Get_Accession_Number <- function(ncbi_id) {
  
  # First, using rentrex package to get summary information for the NCBI id provided.
  summary_info <- entrez_summary(db="nucleotide", id=ncbi_id)
  genbank_accession_number <- summary_info[["accessionversion"]]
  # Save the corresponding GenBank Accession Numbers to a list.
  #Accession <- append(Accession, genbank_accession_number)
  
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
  #Start <- append(Start, Range_Start)
  #End <- append(End, Range_End)
  #Length <- append(Length, Range_Length)
  
  return(genbank_accession_number)
  
}

Get_Accession_Number <- function(ncbi_id) {
  
  # First, using rentrex package to get summary information for the NCBI id provided.
  summary_info <- entrez_summary(db="nucleotide", id=ncbi_id)
  genbank_accession_number <- summary_info[["accessionversion"]]
  # Save the corresponding GenBank Accession Numbers to a list.
  #Accession <- append(Accession, genbank_accession_number)
  
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
  #Start <- append(Start, Range_Start)
  #End <- append(End, Range_End)
  #Length <- append(Length, Range_Length)
  
  return(genbank_accession_number)
  
}













    
    # Next, get the list of titles.
    nico_summs <- entrez_summary(db="nucleotide", id=nico_search$ids)
    titles <- extract_from_esummary(nico_summs, "title")
    # Save the titles to a list.
    Titles <- unname(titles)
    # Add a column for gene name.
    df_len <- length(Titles)
    Gene.Name <- rep(i,df_len)
    # Create a data frame with ncbi ids and titles.
    dfIDs_and_Titles <- data.frame(IDs,Accession,Titles,Gene.Name,Start,End,Length)
    #Append to dataframe.
    mydf <- rbind(mydf,dfIDs_and_Titles)
    
  }
  
  # If no hits were found then ....
  
  else {
    #Append the name.
    mylist <- append(mylist,i)
  }
  
}
)

# 17 minutes to run For Loop ! NoOoOo! Potentially, change this! Too Slow!

# Write this data to a file.
# write.table(mydf,"Results_From_NCBI.txt",sep="\t",row.names=FALSE)

View(mydf)
View(mylist)