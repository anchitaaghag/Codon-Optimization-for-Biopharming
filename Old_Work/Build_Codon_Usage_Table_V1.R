# Codon Optimization for Biopharming: Building a Codon Usage Table
# Version 1 - Using UniProt & BLAST
# 22 June 2022
# Anchitaa Ghag

#### PERSONAL NOTE: FIXES NEEDED ####

# 1) Currently running two For Loops take ~15-20 minutes each. This is not computationally efficient. As discussed with Dr. AHW, python is more efficient than R when it comes to For Loops. Potentially change these For Loops to a function and lapply or sapply for each entry.
# 2) BLAST can either be run via SHARCNET or remotely using the download_once() function. It is more faster to run blast via ComputeCanada, thus, the instruction for this section will be in a seperate file.

#### 00 OVERVIEW - VERSION 1 ####

# Version 1 of the script uses a set of protein names and/or sequences from UniProt to query NCBI for corresponding coding sequences.
# In addition, this script bypasses the need to run multiple entrez searches and BLAST, since all of the information directly comes from Entrez Programming Utilities (E-utilities).
# This script is either dependent on prior installation of BLAST+ and the megablastr R packages (if blast is being run via this R script) or this script will need BLAST-ing of sequences via sharcnet.
# This script is the better choice if it is necessary to have protein GO annotation information or if it would be more beneficial to go from a set of protein sequences to CDS to a codon usage table (i.e. if there are furture unknown proteins sequences that would benefit from the additional BLAST section of this script).

#### 01 INSTALL PACKAGES & DOWNLOAD LIBRARIES ####

# First, set the working directory by running the following lines.
# setwd()
# getwd()

# This script requires the following packages and libraries.
# If these packages have not yet been installed please remove the "#" and run the following lines.

# install.packages("devtools")
# install.packages("BiocManager")
# install.packages("rentrez")
# install.packages("seqinr")
# install.packages("stringr")
# install.packages("tidyverse")
# install.packages("XML")

# BiocManager::install("Biostrings")
# devtools::install_github("HajkD/metablastr", build_vignettes = TRUE, dependencies = TRUE)

library("Biostrings")
library("rentrez")
library("seqinr")
library("stringr")
library("TeachingDemos") # Using this to label all the outliers in a boxplot.
library("tidyverse")
library("XML") # Using this to parse HTML Kazuza's codon usage table

# Load additional functions that are used in this script.

source("~/Major_Research_project_2022/06_Code/Codon-Optimization-for-Biopharming/Functions.R")
source("https://raw.githubusercontent.com/talgalili/R-code-snippets/master/boxplot.with.outlier.label.r") # Load the function

#### 02 IMPORT EXCEL FILE ####

# A list of the proteins in Nicotiana benthamiana, protein sequences, and corresponding information can be obtained from the UniProt database (Bateman et. al., 2021).
# Data orignially downloaded on 28 May 2022. Re-downloaded on 22 June 2022 from: 

# https://www.uniprot.org/uniprot/?query=organism%3A%22Nicotiana%20benthamiana%20%5B4100%5D%22&columns=id%2Centry%20name%2Creviewed%2Cprotein%20names%2Cgenes(PREFERRED)%2Corganism%2Clength%2Cfragment%2Cdatabase(EMBL)%2Csequence%2Ccitation&sort=score

# Load the excel file with protein sequences downloaded from UniProt.

dfData <- as.data.frame(readxl::read_xlsx("~/Major_Research_project_2022/06_Code/Data_Files/UniProt_Data.xlsx"))

# Reformat column names to ensure there are no spaces.

names(dfData) <- make.names(names(dfData),unique(TRUE))

# View the data frame.

# View(dfData)

#### 03 SUMMARY INFORMATION ####

# Check the type of class and summary information.

class(dfData)

# All character data except for one column specifiying the protein lengths.

summary(dfData)

# Check if all entries are for N. benthamiana.

table(dfData$Organism == "Nicotiana benthamiana") # Yes. All 866 entries are for the study species.

# How many entries have been reviewed?

length(dfData$Gene.names) - sum(dfData$Status == "unreviewed") # Only 16 entries have been reviewed.

# Any missing sequence entries?

sum(is.na(dfData$Sequence)) # 0

# Any missing gene name entries?

sum(is.na(dfData$Gene.names)) # 304 missing gene names.

# Any missing protein name entries?

sum(is.na(dfData$Protein.names)) # 0 

#### 04 CLEANING UP DATA FRAME ####

# Subset the columns we want.
# Remove "Organism" column since it has been verified above.
# Remove" Entry.name" column since the "Entry" column provides similar information.
# Remove "Status" column since there isn't enough of "reviewed" entries that can be subset into its own data frame.

colnames(dfData)

dfData.sub <- dfData[c("Entry", "Protein.names","Gene.names...primary..","Length","Fragment","Cross.reference..EMBL.","Sequence","PubMed.ID")]

# Next, clean up column names.

new_col_names <- c("UniProt_Entry","Protein_Name","Gene_Name","Protein_Sequence_Length", "Is_A_Fragment","EMBL_Accession","Protein_Sequence","PubMed_ID")
colnames(dfData.sub) <- new_col_names
rm(new_col_names)

#### 05 DEALING WITH DUPLICATES ####

# Clean-up the names of proteins.

Protein_Names <- unlist(dfData.sub$Protein_Name) %>%
  str_replace_all("\\(Fragment\\)","") %>% # Remove indications of a fragment of protein sequence since this is already indicated by the  "Is_A_Fragment" column e.g (Fragment)
  str_replace_all("\\(EC.*\\)","") # Remove KEGG Entry codes e.g. (EC 3.2.1.52) from protein names

# Add these back to the column.

dfData.sub["Protein_Name"] <- Protein_Names
# Any duplicate protein names?

length(dfData.sub$Protein_Name) - length(unique(dfData.sub$Protein_Name)) # 173 duplicate entries.

# There appears to be 693 unique gene names.

# Find duplicated entries.

Protein_Names <- unlist(dfData.sub$Protein_Name)
Duplicate_Names <- Protein_Names[ Protein_Names %in% Protein_Names[duplicated(Protein_Names)]]

# https://stackoverflow.com/questions/16905425/find-duplicate-values-in-r

NAs <- is.na(Duplicate_Names) 
dup <- Duplicate_Names[!NAs] 
# https://stackoverflow.com/questions/57832161/how-to-remove-na-in-character-vector-in-r

# Final list of gene names that have multiple entries.

unique_duplicates <- unique(dup)

rm(Duplicate_Names,NAs, dup)

# 88 proteins with multiple entries. Can these be filtered further?
# It would be constrictive to ensure that these associated sequences and information makes sense.

# Subset duplicates into a separate data frame.

dfData_Duplicates <- data.frame()

for (i in unique_duplicates) {
  matching_row <- dfData.sub %>%
    filter(Protein_Name == i)
  dfData_Duplicates <- rbind(dfData_Duplicates, matching_row)
}

rm(unique_duplicates, matching_row, i)
View(dfData_Duplicates)

# Filtering of duplicates. 
# There are several duplicates for a single gene in dfData_Duplicates.
# For a number of these genes, there is one full sequence and several "duplicates" fragmented sequences. e.g. CMT3
# It is expected that the full protein sequence will be longer than any one fragmented sequence.
# In addition, the full sequence will most likely begin with a "M" rather than another amino acid.
# For other entries, there is one sequence longer than the others. e.g. Hsp90 has 3 sequences with the lengths 80.81, and 80.
# In this case, only the longest sequence will be retained since there may be more information present that is valuable to downstream analysis. 
# If the sequence lengths only differ by 1 - 2 amino acids, it is expected that BLAST or aligning of sequences will not be grossly impacted.

#Steps:
# 1. Subset the duplicates from dfData
# 2. Remove these duplicated entries from dfData
# 3. Sort the duplicated entries by sequence length.
# Remove duplicates using function. This will retain the first entry (a.k.a. the longest sequence) and remove the succeeding entries.

# Sort dfData_Duplicates by sequence length.

dfDataSorted <- dfData_Duplicates[order(-dfData_Duplicates$Protein_Sequence_Length),] 

#https://stackoverflow.com/questions/62075537/r-pipe-in-does-not-work-with-stringrs-str-extract-all

terms <- duplicated(dfDataSorted$Protein_Name) %>%
  {str_replace(.,"FALSE","No") %>%
      str_replace(.,"TRUE","Yes")}

dfDataSorted["Is.A.Duplicate"] <- terms
rm(terms)

dfKeep <- subset(dfDataSorted, Is.A.Duplicate == "No") # Subset the column for unique entries.
dfKeep <- dfKeep[1:(length(dfKeep)-1)] # Remove the last column indication is the entry is a duplicate or not.

rm(dfDataSorted)

# Add back to original data frame.
# https://stackoverflow.com/questions/17338411/delete-rows-that-exist-in-another-data-frame

dfNew <- dfData.sub[!(dfData.sub$UniProt_Entry %in% dfData_Duplicates$UniProt_Entry),] # Remove all the duplicates from the complete data frame.

dfFinal <- rbind(dfNew,dfKeep) # Add back only the filtered non-duplicated enties that are in dfKeep

rm(Gene_Names,dfData_Duplicates,dfData_No_LecRK,dfKeep,dfNew)

# Confirm we have unique protein names.

table(duplicated(dfFinal$Protein_Name))

#### 06 EXAMINING THE PROTEIN LENGTHS ####

# https://r-charts.com/distribution/add-points-boxplot/

# Are the protein lengths normally distributed?
# Alternatively, Are there any peptides present or abnormalities in lengths?

hist(dfFinal$Protein_Sequence_Length)
shapiro.test(dfFinal$Protein_Sequence_Length)

# W = 0.82433, p-value < 2.2e-16
# Data is not normally distributed.

boxplot(dfFinal$Protein_Sequence_Length,
        col = "white",
        ylab = "Length of Sequence (# of amino acids)", 
        xlab = "Protein Sequences from Uniprot", 
        main = "Boxplot of Protein Sequence Length (no gene information)")

length(boxplot(dfFinal$Protein_Sequence_Length)$out) # 29 Outliers.

#https://www.r-statistics.com/2011/01/how-to-label-all-the-outliers-in-a-boxplot/
#install.packages("TeachingDemos")
library("TeachingDemos")
source("https://raw.githubusercontent.com/talgalili/R-code-snippets/master/boxplot.with.outlier.label.r") # Load the function

y <- dfFinal$Protein_Sequence_Length
lab_y <- dfFinal$Protein_Name

boxplot.with.outlier.label(y, lab_y, spread_text = F)
set.seed(1212)
stripchart(dfFinal$Protein_Sequence_Length,              # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 4,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

# Retain these outliers since there is no evidence to exclude them at this point.

rm(lab_y,y)

#### 07 IMPORT IDS FROM NCBI ####

# Will be using the "nucleotide" database to search for CDS.

# Reformat all gene names to the search term style of entrez_search() function.
# e.g. "Nicotiana benthamiana[Organism] AND CPP1[Gene]"

Reformat <- function(text, wrapper1, wrapper2) {
  
  solution <- c(wrapper1,text,wrapper2)
  results <- paste(solution, collapse = " ")
}

Pro_Names <- dfFinal$Protein_Name[-64]

# Create an empty data frame.
mydf = data.frame()
# Create an empty list to hold the names of gene which has no hits.
mylist = list()

system.time(
  for (i in Pro_Names) {
    # Create the search term.
    wrapping_txt_begin <- "Nicotiana benthamiana[ORGN] AND" 
    wrapping_txt_end <- "[PROT]"
    wrapped_str <- Reformat(i,wrapping_txt_begin,wrapping_txt_end)
    search_term <- str_trim(wrapped_str)
    # Search for entries that match the search term.
    nico_search <- entrez_search(db="nucleotide", term=search_term, retmax = 2000, use_history=FALSE)
    # Count the number of hits entrez_search found. This will indicate if we can continue or not.
    checkpoint <- nico_search[["count"]]
    
    # If more than 0 hits were found then ....
    
    if (checkpoint > 0) {
      # Save the IDs to a list.
      IDs <- nico_search[["ids"]]
      # Next, get the list of titles.
      nico_summs <- entrez_summary(db="nucleotide", id=nico_search$ids)
      titles <- extract_from_esummary(nico_summs, "title")
      # Save the titles to a list.
      Titles <- unname(titles)
      # Add a column for gene name.
      df_len <- length(Titles)
      Gene.Name <- rep(i,df_len)
      # Create a data frame with ncbi ids and titles.
      dfIDs_and_Titles <- data.frame(IDs,Titles,Gene.Name)
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

# How many hits did we get (assuming one per protein)?

length(unique(mydf$Gene.Name)) # 539 hits (one per protein)

length(dfFinal$Protein_Name) - length(unique(mydf$Gene.Name)) # 154 Proteins left to search

rm(checkpoint,df_len,titles,Titles,nico_summs,nico_search,search_term,wrapped_str,wrapping_txt_begin,wrapping_txt_end,dfIDs_and_Titles)

# Filter the data by information in the title.

Titles_Filtered <- mydf$Titles %>%
  str_extract("Nicotiana benthamiana.*") %>% # Keep only N. benthamiana
  str_extract(".*complete cds") # Keep only complete cds (i.e. remove partial cds)

mydf["Filtered_Title"] <- Titles_Filtered

# Remove NAs

mydf <- na.omit(mydf)

# How many genes have multiple entries?

nrow(mydf) - length(unique(mydf$Protein_Name)) #86

# Find duplicated entries.

Gene_Names <- unlist(mydf$Protein_Name)

Duplicate_Names <- Gene_Names[ Gene_Names %in% Gene_Names[duplicated(Gene_Names)] ]

# https://stackoverflow.com/questions/16905425/find-duplicate-values-in-r

NAs <- is.na(Duplicate_Names) 
dup <- Duplicate_Names[!NAs] 
# https://stackoverflow.com/questions/57832161/how-to-remove-na-in-character-vector-in-r

# Final list of gene names that have multiple entries.

unique_duplicates <- unique(dup)

rm(Duplicate_Names,NAs, dup)

# 13 genes with multiple entries. Can these be filtered further?
# It would be constrictive to ensure that these associated sequences and information makes sense.

# Subset duplicates into a separate data frame.

dfData_Duplicates <- data.frame()

for (i in unique_duplicates) {
  matching_row <- mydf %>%
    filter(Protein_Name == i)
  dfData_Duplicates <- rbind(dfData_Duplicates, matching_row)
}

rm(unique_duplicates, matching_row, i)
View(dfData_Duplicates)

# Filtering of duplicates. 
# Remove duplicates using function. This will retain the first entry (a.k.a. the longest sequence) and remove the succeeding entries.

# Sort dfData_Duplicates by sequence length. 

dfDataSorted <- dfData_Duplicates[order(-dfData_Duplicates$CDS_Length),] 

#https://stackoverflow.com/questions/62075537/r-pipe-in-does-not-work-with-stringrs-str-extract-all

terms <- duplicated(dfDataSorted$Titles) %>%
  {str_replace(.,"FALSE","No") %>%
      str_replace(.,"TRUE","Yes")}

dfDataSorted["Is.A.Duplicate"] <- terms
rm(terms)

dfKeep <- subset(dfDataSorted, Is.A.Duplicate == "No") # Subset the column for unique entries.
dfKeep <- dfKeep[1:(length(dfKeep)-1)] # Remove the last column indication is the entry is a duplicate or not.

rm(dfDataSorted)

# Add back to original data frame.
# https://stackoverflow.com/questions/17338411/delete-rows-that-exist-in-another-data-frame

dfNew <- dfData[!(dfData$NCBI_ID %in% dfData_Duplicates$NCBI_ID),] # Remove all the duplicates from the complete data frame.

dfFinal <- rbind(dfNew,dfKeep) # Add back only the filtered non-duplicated enties that are in dfKeep

rm(dfData_Duplicates,dfKeep,dfNew)

# Confirm we have unique protein names.

table(duplicated(dfFinal$Titles))

#### 08 IMPORT ADDITIONAL IDS FROM NCBI #####

# NEXT, Might look again with free word association :

# Create an empty data frame.
mydf1 = data.frame()
# Create an empty list to hold the names of gene which has no hits.
mylist1 = list()

system.time(
  
  for (i in mylist) {
    # Create the search term.
    wrapping_txt_begin <- "Nicotiana benthamiana[ORGN] AND" 
    wrapping_txt_end <- "[WORD]"
    wrapped_str <- Reformat(i,wrapping_txt_begin,wrapping_txt_end)
    search_term <- str_trim(wrapped_str)
    # Search for entries that match the search term.
    nico_search <- entrez_search(db="nucleotide", term=search_term, retmax = 2000, use_history=FALSE)
    # Count the number of hits entrez_search found. This will indicate if we can continue or not.
    checkpoint <- nico_search[["count"]]
    
    # If more than 0 hits were found then ....
    
    if (checkpoint > 0) {
      # Save the IDs to a list.
      IDs <- nico_search[["ids"]]
      # Next, get the list of titles.
      nico_summs <- entrez_summary(db="nucleotide", id=nico_search$ids)
      titles <- extract_from_esummary(nico_summs, "title")
      # Save the titles to a list.
      Titles <- unname(titles)
      # Add a column for gene name.
      df_len <- length(Titles)
      Gene.Name <- rep(i,df_len)
      # Create a data frame with ncbi ids and titles.
      dfIDs_and_Titles <- data.frame(IDs,Titles,Gene.Name)
      #Append to dataframe.
      mydf1 <- rbind(mydf1,dfIDs_and_Titles)
      
    }
    
    # If no hits were found then ....
    
    else {
      #Append the name.
      mylist1 <- append(mylist1,i)
    }
    
  }
)

# Filter the data by information in the title.

Titles_Filtered <- mydf1$Titles %>%
  str_extract("Nicotiana benthamiana.*") %>% # Keep only N. benthamiana
  str_extract(".*complete cds") # Keep only complete cds (i.e. remove partial cds)

mydf1["Filtered_Title"] <- Titles_Filtered

# Remove NAs

mydf1 <- na.omit(mydf1)

# How many genes have multiple entries?

nrow(mydf1) - length(unique(mydf1$Gene.Name)) #64


#### 09 AFTER NCBI SEARCHES ####

# How many hits were found through NCBI?

length(dfFinal$Protein_Name) - length(mylist1)

# How many protein names are still left to search? 

length(mylist1)

# 109 Proteins left to BLAST.

#### 10 BLAST VIA SHARCNET #####

# FIXME Add a file/link  here to describe the steps to take in SHARCNET

#### 11 OBTAIN CDS RANGE INFORMATION ####

# Convert the IDs to a list.
ids <- dfFinal$NCBI_ID
# Apply the Get_Accession_Number() function to each element in the list.
system.time(Accession <- lapply(ids, Get_Accession_Number))
# Append the corresponding GenBank accession numbers to the final data frame.
dfFinal["GenBank_Accession"] <- unlist(Accession)

# Repeat this process for obtaining the start, end, and length of the coding sequences.
# Convert the GenBank accession numbers to a list.
gba_acc <- dfFinal$GenBank_Accession
# Apply the Get_CDS ... () functions to each element in the list.
system.time(Start <- lapply(gba_acc, Get_CDS_Start))
system.time(End <- lapply(gba_acc, Get_CDS_End))
system.time(Length <- lapply(gba_acc, Get_CDS_Length))
# Append the corresponding information to the final data frame.
dfFinal["Start_of_CDS"] <- unlist(Start)
dfFinal["End_of_CDS"] <- unlist(End)
dfFinal["GenBank_Length_of_CDS"] <- unlist(Length)
rm(gba_acc,Start,End,Length)

#### 12 OBTAIN CODING SEQUENCES ####

# Use ids to obtain fasta file of coding sequences.

#ID Name.of.Protein Start Stop Length CDS.Sequence

nico_retrive <- entrez_fetch(db="nucleotide", id = ids[1:300], rettype = "fasta")
nico_retrive1 <- entrez_fetch(db="nucleotide", id = ids[301:600], rettype = "fasta")
write(nico_retrive, file="nico_retrive.fasta")
write(nico_retrive1, file="nico_retrive1.fasta")

# Obtain the list of CDS and add to the data frame.

#BiocManager::install("Biostrings")
library("Biostrings")

fastaFile <- readDNAStringSet("nico_retrive.fasta")
fastaFile1 <- readDNAStringSet("nico_retrive1.fasta")
seq1 <- paste(fastaFile)
seq2 <- paste(fastaFile1)
sequences = c(seq1,seq2)
dfFinal["CDS.(Untrimmed)"] <- sequences

#### 13 TRIMMING SEQUENCES ####

# To get rid of flanking non-coding regions in the CDS, need to trim the seqences based on the start, stop, and length columns.

start_positions <- dfFinal$Start_of_CDS
end_positions <- dfFinal$End_of_CDS

# https://stackoverflow.com/questions/6827299/r-apply-function-with-multiple-parameters
trimmed_seqs <- mapply(str_sub,string=sequences,start=start_positions,end=end_positions)

# Test if this has been done successfully by comparing the lengths of the trimmed sequences with the expected lengths (i.e. the Length.of.CDS column in the data frame.)

trimmed_lengths <- unlist(lapply(trimmed_seqs,nchar))
dfFinal["Trimmed_Length"] <- trimmed_lengths
identical(dfFinal$GenBank_Length_of_CDS,dfFinal$Trimmed_Length)

# TRUE. 
# Can proceed and add the trimmed CDS to the data frame.

dfFinal["Trimmed_CDS"] <- trimmed_seqs

#### 14 EXAMINING THE CDS LENGTHS ####

# https://r-charts.com/distribution/add-points-boxplot/

# Are the coding sequence lengths normally distributed?
# Alternatively, Are there any abnormalities in lengths?

hist(dfFinal$Trimmed_Length)
shapiro.test(dfFinal$Trimmed_Length)

# W = 0.84524, p-value < 2.2e-16
# Data is not normally distributed.

boxplot(dfFinal$Trimmed_Length,
        col = "white",
        ylab = "Length of Sequence (# of amino acids)", 
        xlab = "Protein Sequences from Uniprot", 
        main = "Boxplot of Protein Sequence Length (no gene information)")

length(boxplot(dfFinal$Trimmed_Length)$out) # 38 Outliers.

#https://www.r-statistics.com/2011/01/how-to-label-all-the-outliers-in-a-boxplot/

y <- dfFinal$Trimmed_Length
lab_y <- dfFinal$Titles

boxplot.with.outlier.label(y, lab_y, spread_text = F)
set.seed(1212)
stripchart(dfFinal$Trimmed_Length,              # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 4,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

# Retain these outliers since there is no evidence to exclude them at this point.

rm(lab_y,y)

# Minimum of 80 (or more codons) to conduct MILC analysis.
# Remove sequences less than 80 codons or (80*3 =) 240 bases in sequence length.

dfCodingSeqs <- subset(dfFinal, Trimmed_Length >= 240)

#Recheck the trimmed sequences. Are these divisible by three?
# If sequence lengths are divisible by three, we can continue to retain. 
# If not, either check the trimming or discard.

#### 15 GENERATE CODON UGAGE TABLE (CUT) ####

#https://bioconductor.riken.jp/packages/3.8/bioc/vignettes/coRdon/inst/doc/coRdon.html

# FIXME Uisng codRon package functions. Add Commenting.

cds <- DNAStringSet(dfFinal$Trimmed_CDS)
CUTs <- codonTable(cds)
CU <- codonCounts(CUTs)
getlen(CUTs)
row.names(CU) <- dfFinal$Titles

Average.CU.All.Genes <- colMeans(CU)

#(colMeans(CU)/colSums(CU))*1000

#### 16 IMPORT EXISTING CUT FROM KAZUZA ####

# code adapted from: https://stackoverflow.com/questions/24546312/vector-of-most-used-codons-from-table-of-codon-usage-in-r

# Table can be found here:
# https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4100&aa=1&style=GCG

# Using the XML package's htmlParse() function to "read" an HTML file and generate an HTML/XMLInternalDocument class object.

Kazuza <- htmlParse('http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4100&aa=1&style=GCG')

# Next, read this object and convert to a dataframe.

dfKazuza <- read.table(text=xpathSApply(Kazuza, "//pre", xmlValue), 
                       header=TRUE, 
                       fill=TRUE)

# Finally, calculate the most used codons (1 for each amino acid residue) from the data frame.

dfKazuza.max <- group_by(dfKazuza, AmAcid) %>% 
  filter(Number==max(Number)) %>% 
  select(AmAcid, Codon)

#### 17 STATISTICS ####

# Currently in progress.

milc <- MILC(CUTs)
head(milc)

#### 18 REFERENCES ####
