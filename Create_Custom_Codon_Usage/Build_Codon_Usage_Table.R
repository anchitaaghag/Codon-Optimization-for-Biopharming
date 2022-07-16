# Codon Optimization for Biopharming: Building a Codon Usage Table using Coding Sequences
# Anchitaa Ghag

#### PERSONAL NOTE: FIXES NEEDED ####

# At the very beginning, include a section or a way for the user to input ALL search parameters or thresholds (to ensure reproducible code).

#### 01 INSTALL PACKAGES & DOWNLOAD LIBRARIES ####

# First, set the working directory by running the following lines.
# setwd()
# getwd()

# This script requires the following packages and libraries.
# If these packages have not yet been installed please remove the "#" and run the following lines.

# install.packages("ape")
# install.packages("BiocManager")
# install.packages("rentrez")
# install.packages("seqinr")
# install.packages("stringr")
# install.packages("TeachingDemos")
# install.packages("tidyverse")
# install.packages("XML")

# BiocManager::install("Biostrings")
# BiocManager::install("coRdon")

library("ape")
library("Biostrings")
library("coRdon")
library("rentrez")
library("seqinr")
library("stringr")
library("TeachingDemos") # Using this to label all the outliers in a boxplot.
library("tidyverse")
library("XML") # Using this to parse HTML Kazuza's codon usage table

# Load additional functions that are used in this script.

source("https://raw.githubusercontent.com/talgalili/R-code-snippets/master/boxplot.with.outlier.label.r")

#### 02 DATA AQUISITION : IMPORT EXISTING CU FROM KAZUZA ####

# code adapted from: https://stackoverflow.com/questions/24546312/vector-of-most-used-codons-from-table-of-codon-usage-in-r

# The existing codon usage table from Kazuza can be imported from: 
# https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4100&aa=1&style=GCG

# It can also be imported from the "Kazuza_Codon_Usage.txt" file using the following lines of code. File originally created on 24 June 2022.

# Kazuza <- read_table("Kazuza_Codon_Usage.txt")

# Using the XML package's htmlParse() function to "read" an HTML file and generate an HTML/XMLInternalDocument class object.

#Kazuza <- htmlParse('http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4100&aa=1&style=GCG')

#KazuzaTabacum <- htmlParse('http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4097&aa=1&style=GCG')

# Next, read this object and convert to a dataframe.

#dfKazuza <- read.table(text=xpathSApply(Kazuza, "//pre", xmlValue), 
                      # header=TRUE, 
                      # fill=TRUE)

#dfNTabacum <- read.table(text=xpathSApply(KazuzaTabacum, "//pre", xmlValue), 
                      # header=TRUE, 
                      # fill=TRUE)

# List of the codon usage from each CDS that makes up Kazuza CUT is available at: http://www.kazusa.or.jp/codon/current/species/4100
#KazuzaCDS <- read.table("/Users/anchitaa/Major_Research_Project_2022/06_Code/KCC_No_Carot.txt",
                       # header = TRUE)

# Finally, calculate the most used codons (1 for each amino acid residue) from the data frame.

#dfKazuza.max <- group_by(dfKazuza, AmAcid) %>% 
#  filter(Number==max(Number)) %>% 
#  select(AmAcid, Codon)

# Remove all objects no longer needed from the environment.

#rm(Kazuza)

#### 03 DATA AQUISITION : OBTAIN IDS FROM NCBI ####

# A list of the all complete coding sequences in Nicotiana benthamiana and corresponding information can be obtained from the NCBI database.

# FIXME First, define the organism of interest.

# List of databases available to search using EUtils API. This script is using the "nucleotide" database.

entrez_dbs()

# Once, a database has been chosen define the search field. This script is using the [FKEY] parameter to return only coding sequences.

entrez_db_searchable(db="nucleotide")

# Using a web history object. The maximum retmax chaned from defalut 20 to 5000 , can be increaded if it is anticipated that more records will be retrived.

nico_search <- entrez_search(db="nucleotide", 
                             term="Nicotiana benthamiana[Organism] AND CDS[FKEY]", 
                             use_history=TRUE, 
                             retmax = 5000)

# Create the summary information for each title from the entrez search result.
# This script will use version 1.0 or XML (as opposed to the more extensive version 2.0 or JSON). This is because the maximum number of UIDs is 500 for a JSON format output. For longer secuence requests (more than 500), even with a web history object, using a single entrez_summary search will result in a "Too many UIDs in request." error. To ensure that all records can be fetched in one order, regardless of request size, the more limeted summaries of a database record in version 1.0 is being used. This can be changed to version 2.0 and submitted in batches.

nico_summs <- entrez_summary(db="nucleotide", 
                             web_history=nico_search$web_history, 
                             version = "1.0")

# Extract all the titles, NCBI ids, GenBank accession versions, and untrimmed sequence lengths from the esummary as a matrix.

InfoTable <- extract_from_esummary(nico_summs, c("Title","Length","AccessionVersion"))

# To ensure readability, transpose the matrix, to make the "Title","Length","AccessionVersion" fields into columns and convert to a data frame.

dfTemp <- as.data.frame(t(InfoTable))

# Change the column names to be more descriptive.

colnames(dfTemp) <- c("Titles","Untrimmed_Sequence_Length","GenBank_Accession")

# Use the rownames_to_column() from the tibble package to create a NCBI_ID column in the new data frame.

dfTemp <- rownames_to_column(dfTemp, "NCBI_ID")

# Ensure that we have a traditional data frame and not one composed of lists.
# Convert to a traditional data frame format.

dfNCBI <- as.data.frame(lapply(dfTemp, unlist))

# Remove all objects no longer needed from the environment.

rm(InfoTable,dfTemp)

#### 05 DATA AQUISITION : OBTAIN SEQUENCES ####

# Use web history to obtain fasta file of coding sequences.

# Since the entrez_fetch() function

nico_retrive <- entrez_fetch(db="nucleotide", web_history = nico_search$web_history, rettype = "fasta")

write(nico_retrive, file="nico_retrive.fasta")

# Obtain the list of CDS and add to the data frame.

fastaFile <- readDNAStringSet("nico_retrive.fasta")

sequences <- paste(fastaFile)

dfNCBI["Untrimmed_Sequences"] <- sequences

# Remove all objects no longer needed from the environment.

rm(fastaFile, nico_search, nico_summs, nico_retrive, sequences)

#### 04 DATA AQUISITION : OBTAIN ANNOTATIONS ####

# Convert the GenBank accession numbers to a list.
Accessions <- unlist(dfNCBI$GenBank_Accession) 

# Check class. Ensure that the class is a character vector.
class(Accessions) 

# Using the getAnnotationsGenBank() function from the ape package to get the annotations for each sequence from GenBank.

Feature_List <- getAnnotationsGenBank(Accessions, quiet = FALSE) 

# FIXME Add an option to use the feature file for faster script running
# This function takes approximately 9.082817 minutes to run. For a faster script time, the feature list can also be loaded into R from the file using the following lines of code:

# Other options include the biofiles package and the genbankr package functions. However, the ape package function has a faster run time that can also handle large amounts of queries. It also does not require parsing of GenBank files (which may be tedious when dealing with a large number of records).
# Can also be downloaded from NCBI: search "Nicotiana benthamiana"[Organism] AND "CDS"[FKEY] 
# Send to: -> Complete record -> File -> Format: Feature table -> Default order -> Create file

# Convert the list of dataframes to a data frame. 
# https://stackoverflow.com/questions/2851327/combine-a-list-of-data-frames-into-one-data-frame-by-row

dfFeatures <- bind_rows(Feature_List, .id = "GB_Accession")

# write.table(x = dfFeatures, file = "Features_Data_Frame")

# Change the column names to be more descriptive.

colnames(dfFeatures) <- c("GB_Accession","Start_of_CDS","End_of_CDS","Type","Product","Protein_ID","Gene")

# Remove all objects no longer needed from the environment.

rm(Accessions,Feature_List)

#### ADD GENE INFORMATION ####

# The gene name for each entry is distributed over three columns: Product, Protein_ID, and Gene. To ensure that these are formatted correctly, extract the relevant information from each of the columns and merge.

# First, for entries that have a listed gene name, extract the gene names into a character vector.

Gene_Name_List_1 <- dfFeatures$Gene 

# Second, for entries that have a listed gene name in the protein description, extract the gene names into a character vector.

Gene_Name_List_2 <- dfFeatures$Protein_ID %>%
  str_extract("prot_desc.*") %>%
  str_extract(".*;") %>%
  str_replace("prot_desc:", "") %>%
  str_replace(";.*", "") %>%
  str_squish()
    
# Third, for entries that have a listed gene name in the notes, extract the gene names into a character vector.

Gene_Name_List_3 <- dfFeatures$Protein_ID %>%
  str_extract("note.*") %>%
  str_replace("note:", "") %>%
  str_replace(";.*", "") %>%
  str_squish()

# Lastly, for entries that have listed a gene name instead of a protein name, extract the gene names into a character vector.

Gene_Name_List_4 <- dfFeatures$Product

# Lastly, extract names from the entry titles.

df

# Condense these four lists, by order of importance (i.e. as aforementioned above)
  
dfGeneNames["All_Gene_Names"] <- coalesce(x = Gene_Name_List_1, y = Gene_Name_List_2) %>%
  coalesce(y = Gene_Name_List_3) %>%
  coalesce(y = Gene_Name_List_4)

#### SUMMARY ####

# Rename Columns to Nicer Names 
# FIXME Summary Stats here, see previous script work

#### DATA FILTERING : COMBINE ALL INFO INTO ONE DATAFRAME AND SUMMARY INFO HERE ####

# The script now has two data frames as follows:

# 1) dfFeatures - A dataframe holding several features annotated for each sequence.
# 2) dfNCBI - A dataframe holding the untrimmed sequences and version information.

# The number of features can exceed the number of records.
table(duplicated(dfFeatures$GB_Accession)) # FALSE = 332 "extra" features annotated.

# What type of features have been annotated in the records?
table(dfFeatures$Type)

# Several features including coding sequences, mRNA 3' and 5' untranslated regions (UTRs).
# Filter the data frame to retain features marked "CDS", eliminating UTRs, rRNAs, and other non-coding regions.

dfFeatures.sub <- dfFeatures %>%
  filter(Type == "CDS") # Retain only feature "rows" marked as CDS 

# Subsetting dfData by the column names in dfFeatures. This will retain only records marked as "CDS"
  
dfCombined <- subset(dfNCBI, GenBank_Accession %in% dfFeatures.sub$GB_Accession)

# Add the feature columns to the new data frame.

dfData <- cbind(dfCombined, dfFeatures.sub[,2:6])

# Append a new column with the length of the CDS.

dfData["GenBank_Length_of_CDS"] <- ((dfData$End_of_CDS - dfData$Start_of_CDS) + 1)

# Rearrange column names, for readability.

dfData <- dfData[,c("GenBank_Accession", "NCBI_ID", "Titles", "Type", "Untrimmed_Sequence_Length", "Untrimmed_Sequences", "Start_of_CDS", "End_of_CDS", "GenBank_Length_of_CDS", "Protein_ID", "Product")]

# Remove all objects no longer needed from the environment.

rm(dfNCBI, dfFeatures, dfFeatures.sub, dfCombined)

#### 06 DATA FILTERING : FILTER RECORDS ####

# Next, filter the data frame to ensure that only complete coding sequences are retained.
# This is done because the empirical codon counts for each entry necessitates complete sequences.
# If partial sequences are retained, the incomplete codon counts for each sequence would skew the overall frequencies per codon.

Filtered_Titles <- dfData$Titles %>%
  str_replace("Nicotiana benthamiana mRNA for hypothetical protein.*", "") %>% # Remove records for hypothetical proteins
  str_replace("Nicotiana benthamiana clone.*", "") %>% # Remove records for any N. benthamiana clones
  str_extract("Nicotiana benthamiana.*") %>% # Ensure that only N. benthamiana entries are retained
  str_extract(".*mRNA, complete cds.*")  # Keep only complete cds (i.e. remove partial cds) 

# Remove any coding sequences for hypothetical proteins.
# Subsetting dfData by the column names in dfFeatures. This will retain only records marked as "CDS"

dfData.sub <- subset(dfData, Titles %in% Filtered_Titles)

# How many entries (excluding any duplicates)?

length(unique(dfData.sub$Titles)) # 281 Entries.

# Remove all objects no longer needed from the environment.

rm(Filtered_Titles, dfData)

#### 07 DATA FILTERING : REMOVE DUPLICATED TITLES ####

# Find duplicated entries.

Names <- unlist(dfData.sub$Titles)

Duplicate_Names <- Names[ Names %in% Names[duplicated(Names)] ]

# https://stackoverflow.com/questions/16905425/find-duplicate-values-in-r

NAs <- is.na(Duplicate_Names) 
dup <- Duplicate_Names[!NAs] 
# https://stackoverflow.com/questions/57832161/how-to-remove-na-in-character-vector-in-r

# Final list of gene names that have multiple entries.

unique_duplicates <- unique(dup)

#rm(Duplicate_Names,NAs, dup)

length(unique_duplicates)

# 2 records with multiple entries. Can these be filtered further?
# It would be constructive to ensure that these associated sequences and information makes sense.

# Subset duplicates into a separate data frame.

dfData_Duplicates <- data.frame()

for (i in unique_duplicates) {
  matching_row <- dfData.sub %>%
    filter(Titles == i)
  dfData_Duplicates <- rbind(dfData_Duplicates, matching_row)
}

#rm(unique_duplicates, matching_row, i)

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

dfDataSorted <- dfData_Duplicates[order(dfData_Duplicates$Untrimmed_Sequence_Length, decreasing = TRUE),] 

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

dfNew <- dfData.sub[!(dfData.sub$NCBI_ID %in% dfData_Duplicates$NCBI_ID),] # Remove all the duplicates from the complete data frame.

dfFinal <- rbind(dfNew,dfKeep) # Add back only the filtered non-duplicated enties that are in dfKeep

# Confirm we have unique protein names.

table(duplicated(dfFinal$Titles))

# Remove all objects no longer needed from the environment.

rm(dfData_Duplicates,dfKeep,dfNew)

#### 08 DATA FILTERING : TRIM SEQUENCES ####

# To get rid of flanking non-coding regions in the CDS, need to trim the sequences based on the start, stop, and length columns.

# https://stackoverflow.com/questions/6827299/r-apply-function-with-multiple-parameters
Trimmed_Sequences <- mapply(str_sub,
                       string=dfFinal$Untrimmed_Sequences,
                       start=dfFinal$Start_of_CDS,
                       end=dfFinal$End_of_CDS)

# Test if this has been done successfully by comparing the lengths of the trimmed sequences with the expected lengths (i.e. the Length.of.CDS column in the data frame.)

Trimmed_Lengths <- unlist(lapply(Trimmed_Sequences,nchar))
dfFinal["Trimmed_Length"] <- Trimmed_Lengths

identical(dfFinal$GenBank_Length_of_CDS,as.numeric(Trimmed_Lengths))

# TRUE. 
# Can proceed and add the trimmed CDS to the data frame.

dfFinal["Trimmed_CDS"] <- Trimmed_Sequences

# Remove all objects no longer needed from the environment.

rm(Untrimmed_Sequences,Start_Positions,End_Positions,Trimmed_Sequences,Trimmed_Lengths)

#### GET GENE NAMES ###

Gene_Names <- dfFinal$Protein_ID %>%
  str_match("prot_desc:\\s*(.*?)\\s*;") %>%
  str_replace("prot_desc:","") %>%
  str_replace("Nb","") %>%
  str_replace("\\;","") %>%
  str_trim()
Gene_Names[5] <- "GRP7"

GenBank_Protein_ID <- dfFinal$Protein_ID %>%
  str_match(".*gb\\|([^|]+)")

#### 09 DATA FILTERING : EXAMINE THE CDS LENGTHS ####

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

# A minimum of 80 (or more codons) is required to conduct MILC analysis.
# Remove sequences less than 80 codons or (80*3 =) 240 bases in sequence length.

dfCodingSeqs <- subset(dfFinal, Trimmed_Length >= 240)

#Recheck the trimmed sequences. Are these divisible by three?
  # If sequence lengths are divisible by three, we can continue to retain. 
  # If not, either check the trimming or discard.

# Remove all objects no longer needed from the environment.

rm(lab_y,y)

#### 10 GENERATE CODON USAGE TABLES ####

#https://bioconductor.riken.jp/packages/3.8/bioc/vignettes/coRdon/inst/doc/coRdon.html

# For Kazuza

CUT_Kazuza <- codonTable(KazuzaCDS)

# FIXME Uisng codRon package functions. Add Commenting.

cds <- DNAStringSet(dfCodingSeqs$Trimmed_CDS)
CUTs <- codonTable(cds)
CU <- codonCounts(CUTs)
getlen(CUTs)
row.names(CU) <- dfCodingSeqs$Titles

Average.CU.All.Genes <- colMeans(CU)

#(colMeans(CU)/colSums(CU))*1000

write.csv(x=dfCodingSeqs, file="dfCodingSeqs.csv")
write.csv(x=dfKazuza, file="dfKazuza.csv")
write.csv(x=dfNTabacum, file="dfNTabacum.csv")

#### 11 STATISTICS ####

# https://www.bioconductor.org/packages/devel/bioc/vignettes/coRdon/inst/doc/coRdon.html#calculate-cu-bias

# Currently in progress.

milc <- MILC(CUTs)
head(milc)

# "Other CU statistics can be calculated in the same way as MILC(), using one of the functions: B(), MCB(), ENCprime(), ENC() or SCUO(). Note however, that when calculating ENC and SCUO, one doesn’t need to provide a subset of refe"

SCUO(CUTs)
ENC(CUTs)
B(CUTs)
MCB(CUTs)
ENCprime(CUTs)

# Next, compare the CU bias for every coding sequence between the created CUT and Kazuza through visualizations.
# 

library(ggplot2)
xlab <- "MILC distance from sample centroid"
ylab <- "MILC distance from ribosomal genes"

MILC_Created_CUT <- MILC(CUTs, ribosomal = TRUE)
Bplot(x = "self", y = "self", data = MILC_Created_CUT) +
  labs(x = xlab, y = ylab)

# Visualization of Codon Usage - Existing vs. Current

# Code Adapted From: https://www.bioconductor.org/packages/devel/bioc/vignettes/coRdon/inst/doc/coRdon.html#calculate-cu-bias

# Use the intraBplot() function in the coRdon package and our two codonTable objects to plot a  plot of codon usage distances between existing vs. created codon tables.
#Intra-samples Karlin B plot

intraBplot(x = CUT_Kazuza, 
           y = CUTs, 
           names = c("Kazuza", "Current"), 
           variable = "MILC", 
           size = 3, 
           alpha = 1.0) + 
           ggtitle("B Plot of Existing v.s. Current CU Distances")

#hi <- lapply(list, function)

formatted_cds <- s2c((unlist(dfCodingSeqs$Trimmed_CDS[1])))
s2c("hi")
uco(seq=formatted_cds,
    frame = 0,
    index = "rscu",
    #as.data.frame = TRUE,
    NA.rscu = NA)

## Make a coding sequence from this:
(cds <- s2c(paste(words(), collapse = "")))
uco(cds, index = "rscu")


rcds <- read.fasta(file = system.file("sequences/malM.fasta", package = "seqinr"))[[1]]
uco( rcds, index = "freq")
uco( rcds, index = "eff")
uco( rcds, index = "rscu")
uco( rcds, as.data.frame = TRUE)


#### Comparison of all 3 per amino acid ####3


#### 11 REFERENCES ####


