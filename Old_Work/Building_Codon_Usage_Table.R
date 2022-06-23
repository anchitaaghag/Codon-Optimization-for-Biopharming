# Codon Optimization for Biopharming: Building a Codon Usage Table
# 22 June 2022
# Anchitaa Ghag

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
library("tidyverse")
library("XML")

#### 02 IMPORT EXCEL FILE ####

# A list of the proteins in Nicotiana benthamiana, protein sequences, and corresponding information can be obtained from the UniProt database (Bateman et. al., 2021).
# Data orignially downloaded on 28 May 2022. Re-downloaded on 22 June 2022 from: 

# https://www.uniprot.org/uniprot/?query=organism%3A%22Nicotiana%20benthamiana%20%5B4100%5D%22&columns=id%2Centry%20name%2Creviewed%2Cprotein%20names%2Cgenes(PREFERRED)%2Corganism%2Clength%2Cfragment%2Cdatabase(EMBL)%2Csequence%2Ccitation&sort=score

# Load the excel file with protein sequences downloaded from UniProt.

dfData <- as.data.frame(readxl::read_xlsx("~/Major_Research_project_2022/06_Code/UniProt_Data.xlsx"))

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

# How many duplicated protein names? 

length(dfData$Protein.names) - length(unique(dfData$Protein.names)) # 171

#### CLEANING UP DATA FRAME ####

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


#### REFERENCES ####






#
