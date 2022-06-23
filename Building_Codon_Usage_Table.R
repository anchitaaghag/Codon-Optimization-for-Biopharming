# Codon Optimization for Biopharming: Building a Codon Usage Table
# Version 2 - Using NCBI Only
# 22 June 2022
# Anchitaa Ghag

#### PERSONAL NOTE: FIXES NEEDED ####

# 1) Currently, there are 4 entrez searches being performed since the maximum number of hits I can query per search is ~300. Is there a simpler way to code this without creating 4 search objects, 4 summary objects, 4 titles, and 4 id lists? [i.e perform this search using 1 search object, 1 summary obj ... ]
# 2) At the very beginning, include a section or a way for the user to input ALL search parameters or thresholds (to ensure reproducible code).
# 3) Try to amend the Get_CDS...() functions created to avoid outputting every step (may save time?)
# 4) May want to switch sections for flow. E.g. moving the section of importing of the Kazuza table up so that you can import all information once at the beginning of the script.

#### 00 OVERVIEW - WHY A VERSION 2? ####

# Version 2 of the script features a much faster run time, more easily reproducible code (i.e. less hardcoding of URLs/Files), and a more simple and straightforward approach.
# Compared to the previous version, Version 2 is more efficient because it queries a single database - NCBI, rather than two databases. The overall number of coding sequences is not impacted.
# By omitting the back-and-forth between two databases, a considerable amount of time is saved.
# In addition, this script bypasses the need to run multiple entrez searches and BLAST, since all of the information directly comes from Entrez Programming Utilities (E-utilities).
# One dependency of prior installation of BLAST+ is omitted, as well as a few R packages.

#### 01 INSTALL PACKAGES & DOWNLOAD LIBRARIES ####

# First, set the working directory by running the following lines.
# setwd()
# getwd()

# This script requires the following packages and libraries.
# If these packages have not yet been installed please remove the "#" and run the following lines.

# install.packages("BiocManager")
# install.packages("rentrez")
# install.packages("seqinr")
# install.packages("stringr")
# install.packages("TeachingDemos")
# install.packages("tidyverse")
# install.packages("XML")

# BiocManager::install("Biostrings")
# BiocManager::install("coRdon")

library("Biostrings")
library("coRdon")
library("rentrez")
library("seqinr")
library("stringr")
library("TeachingDemos") # Using this to label all the outliers in a boxplot.
library("tidyverse")
library("XML") # Using this to parse HTML Kazuza's codon usage table

# Load additional functions that are used in this script.

source("~/Major_Research_project_2022/06_Code/Codon-Optimization-for-Biopharming/Functions.R")
source("https://raw.githubusercontent.com/talgalili/R-code-snippets/master/boxplot.with.outlier.label.r") # Load the function

#### 02 IMPORT IDS FROM NCBI ####

# A list of the all complete coding sequences in Nicotiana benthamiana and corresponding information can be obtained from the NCBI database.

# FIXME First, define the organism of interest.

# List of databases available to search using EUtils API. This script is using the "nucleotide" database.

entrez_dbs()

# Once, a database has been chosen define the search field. This script is using the [FKEY] parameter to return only coding sequences.

entrez_db_searchable(db="nucleotide", config = NULL)

nico_search <- entrez_search(db="nucleotide", term="Nicotiana benthamiana[Organism] AND CDS[FKEY]", use_history=FALSE, retmax = 98000)

# Create the summary information for each title from the entrez search result.

nico_summs1 <- entrez_summary(db="nucleotide", id=nico_search$ids[1:300])
nico_summs2 <- entrez_summary(db="nucleotide", id=nico_search$ids[301:600])
nico_summs3 <- entrez_summary(db="nucleotide", id=nico_search$ids[601:900])
nico_summs4 <- entrez_summary(db="nucleotide", id=nico_search$ids[901:1200])

# Extract all the titles.

titles1 <- extract_from_esummary(nico_summs1, "title")
titles2 <- extract_from_esummary(nico_summs2, "title")
titles3 <- extract_from_esummary(nico_summs3, "title")
titles4 <- extract_from_esummary(nico_summs4, "title")

# Extract all the NCBI ids.

# Create an empty list to hold the names of gene which has no hits.
ID1 = list()
ID2 = list()
ID3 = list()
ID4 = list()

LEN1 = list()
LEN2 = list()
LEN3 = list()
LEN4 = list()

for (j in nico_summs1) {
 uid <- j[["uid"]]
 slen <- j[["slen"]]
  #Append the id to the list.
  ID1 <- append(ID1,uid)
  # Append length of sequences to the list.
  LEN1 <- append(LEN1,slen)
}

for (j in nico_summs2) {
  uid <- j[["uid"]]
  slen <- j[["slen"]]
  #Append the id to the list.
  ID2 <- append(ID2,uid)
  # Append length of sequences to the list.
  LEN2 <- append(LEN2,slen)
}

for (j in nico_summs3) {
  uid <- j[["uid"]]
  slen <- j[["slen"]]
  #Append the id to the list.
  ID3 <- append(ID3,uid)
  # Append length of sequences to the list.
  LEN3 <- append(LEN3,slen)
}

for (j in nico_summs4) {
  uid <- j[["uid"]]
  slen <- j[["slen"]]
  #Append the id to the list.
  ID4 <- append(ID4,uid)
  # Append length of sequences to the list.
  LEN4 <- append(LEN4,slen)
}
  
# Save the ids to one list.

NCBI_ID <- unlist(c(ID1,ID2,ID3,ID4))

rm(ID1,ID2,ID3,ID4)

# Save the lengths to one list.

CDS_Length <- unlist(c(LEN1,LEN2,LEN3,LEN4))

rm(LEN1,LEN2,LEN3,LEN4)

# Save the titles to one list.

t<-c(titles1,titles2,titles3,titles4)

Titles <- unname(t)

rm(t,titles1,titles2,titles3,titles4)

# Add the ids and titles into a dataframe.

dfIDs_and_Titles <- data.frame(NCBI_ID,Titles,CDS_Length)

# Filter the data by information in the title.

T_Filtered <- dfIDs_and_Titles$Titles %>%
  str_extract("Nicotiana benthamiana.*") %>% # Keep only N. benthamiana entries
  str_extract(".*complete cds") %>% # Keep only complete cds (i.e. remove partial cds) 
  str_extract(".*mRNA.*") %>% # Keep only mRNA records (i.e. remove any gene or unknown records) 
  str_replace("Nicotiana benthamiana mRNA for hypothetical protein.*", "") %>% # Remove records for hypothetical proteins
  str_replace("Nicotiana benthamiana clone.*", "") # Remove records for any N. benthamiana clones

# Add back to dataframe.

dfIDs_and_Titles["Titles"] <- T_Filtered
  
# Remove "NA"s from the list.

dfFiltered_IDs_Titles <- na.omit(dfIDs_and_Titles)

# How many entires (excluding any duplicates)?

dfData <- subset(dfFiltered_IDs_Titles, Titles != "") # Removing any empty rows.

length(unique(dfData$Titles)) # 585 Entries.

#### 03 REMOVE DUPLICATED TITLES ####

# # Because entrez-searches are not case sensetive, many different enties are possible for the same gene.
# For example, while AGO1b and AGO1B are the same gene, their diffrent spelling leads to 4 entries returned instead of one.

# Find duplicated entries.

Names <- unlist(dfData$Titles)

Duplicate_Names <- Names[ Names %in% Names[duplicated(Names)] ]

# https://stackoverflow.com/questions/16905425/find-duplicate-values-in-r

NAs <- is.na(Duplicate_Names) 
dup <- Duplicate_Names[!NAs] 
# https://stackoverflow.com/questions/57832161/how-to-remove-na-in-character-vector-in-r

# Final list of gene names that have multiple entries.

unique_duplicates <- unique(dup)

#rm(Duplicate_Names,NAs, dup)

length(unique_duplicates)

# 4 proteins with multiple entries. Can these be filtered further?
# It would be constructive to ensure that these associated sequences and information makes sense.

# Subset duplicates into a separate data frame.

dfData_Duplicates <- data.frame()

for (i in unique_duplicates) {
  matching_row <- dfData %>%
    filter(Titles == i)
  dfData_Duplicates <- rbind(dfData_Duplicates, matching_row)
}

#rm(unique_duplicates, matching_row, i)
View(dfData_Duplicates)

# Filtering of duplicates. 
# There are several duplicates for a single gene in dfData_Duplicates.
# Retain the longest sequence. If there are flanking non-coding regions in the sequence, it can be trimmed below.
# In this case, only the longest sequence will be retained since there may be more information present that is valuable to downstream analysis. 
#Steps:
# 1. Subset the duplicates from dfData
# 2. Remove these duplicated entries from dfData
# 3. Sort the duplicated entries by sequence length.
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

#### 04 OBTAIN CDS RANGE INFORMATION ####

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

#### 05 OBTAIN CODING SEQUENCES ####

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

#### 06 TRIMMING SEQUENCES ####

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

#### 07 EXAMINING THE CDS LENGTHS ####

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

#### 08 GENERATE CODON UGAGE TABLE (CUT) ####

#https://bioconductor.riken.jp/packages/3.8/bioc/vignettes/coRdon/inst/doc/coRdon.html

# FIXME Uisng codRon package functions. Add Commenting.

cds <- DNAStringSet(dfFinal$Trimmed_CDS)
CUTs <- codonTable(cds)
CU <- codonCounts(CUTs)
getlen(CUTs)
row.names(CU) <- dfFinal$Titles

Average.CU.All.Genes <- colMeans(CU)

#(colMeans(CU)/colSums(CU))*1000

#### 09 IMPORT EXISTING CUT FROM KAZUZA ####

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

#### 10 STATISTICS ####

# Currently in progress.

milc <- MILC(CUTs)
head(milc)


#### 11 REFERENCES ####








#