# Codon Optimization for Biopharming: Building an Updated Codon Usage Table for Nicotiana benthaminana
# Anchitaa Ghag

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
library("tidyverse")
library("XML")

#### 02 INPUT PARAMETERS AND THRESHOLDS ####

# This script can also be run for a different organism by changing the following search parameters or thresholds.
# First, specify the species of interest. All coding sequences "CDS" will be searched for the species of interest Nicotiana benthamiana. 

Entrez_Search_Term <- "Nicotiana benthamiana[Organism] AND CDS[FKEY]"

# The maximum retmax has changed from the default 20 to 5000, this can be increased if it is anticipated that more records will be retrieved (i.e with a more well-studied species such as Nicotiana tabacum).
# Currently, for N. benthamiana, the total number of hits found are 1084. Therefore, this threshold has been set to 5000 (i.e. well above 1084).

Maximum_Sequences_To_Retrieve <- 5000

#### 03 DATA AQUISITION ALTERNATIVE: LOAD DATA FROM FILE ####

# If desired, the data can be loaded from a file by running the following lines of code and skipping Sections 04 to 07 below. This may be useful to reproduce the results of the original script. Data originally downloaded on 17 July 2022.
# dfNCBI <- 
# dfFeatures <-

# If not, the data from NCBI can be obtained by running the following script below (Sections 04 to 07).
# The latter may be particularly helpful if additional sequences have been deposited in the GenBank/NCBI database or if an updated codon usage table is required.

#### 04 DATA AQUISITION : OBTAIN IDS FROM NCBI ####

# A list of the all complete coding sequences in Nicotiana benthamiana and corresponding information can be obtained from the NCBI database.
# The full list of databases available to search using EUtils API can be found using the entrez_dbs() function. This script is using the "nucleotide" database.
# The code in this section was adapted from the Rentrez package vignette available from: https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html. (Winter, 2020).

entrez_dbs()

# Once a database has been chosen,the search field needs to be defined. This has been completed above in the user input section. This script is using the [FKEY] parameter to return only coding sequences. (Winter, 2020).

entrez_db_searchable(db="nucleotide")

# Search the NCBI database and store the resulting hits as a web history object. (Winter, 2020).

nico_search <- entrez_search(db="nucleotide", 
                             term= Entrez_Search_Term, 
                             use_history=TRUE, 
                             retmax = Maximum_Sequences_To_Retrieve)

# Create the summary information for each title from the entrez search result. (Winter, 2020).
# This script will use version 1.0 or XML (as opposed to the more extensive version 2.0 or JSON). This is because the maximum number of UIDs is 500 for a JSON format output. For longer sequence requests (more than 500), even with a web history object, using a single entrez_summary search will result in a "Too many UIDs in request." error. To ensure that all records can be fetched in one order, regardless of request size, the more limited summaries of a database record in version 1.0 is being used. This can be changed to version 2.0 and submitted in batches.

nico_summs <- entrez_summary(db="nucleotide", 
                             web_history=nico_search$web_history, 
                             version = "1.0")

# Extract all the titles, NCBI ids, GenBank accession versions, and untrimmed sequence lengths from the esummary as a matrix. (Winter, 2020).

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

rm(InfoTable,dfTemp, Entrez_Search_Term, Maximum_Sequences_To_Retrieve)

#### 05 DATA AQUISITION : OBTAIN SEQUENCES ####

# The code in this section was adapted from the Rentrez package vignette available from: https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html. (Winter, 2020).
# Use web history created above to coding sequences from the nucleotide database in a fasta format. (Winter, 2020).

nico_retrive <- entrez_fetch(db="nucleotide", 
                             web_history = nico_search$web_history, 
                             rettype = "fasta")

# Next, write the fasta sequences to a fasta file. (Winter, 2020).

write(nico_retrive, file="nico_retrive.fasta")

# Obtain the list of coding sequences and add them to the data frame. (Winter, 2020).

dfNCBI["Untrimmed_Sequences"] <- paste(readDNAStringSet("nico_retrive.fasta"))

# Remove all objects no longer needed from the environment.

rm(nico_search, nico_summs, nico_retrive)

#### 06 DATA AQUISITION : OBTAIN ANNOTATIONS ####

# Convert the GenBank accession version numbers to a list.

Accessions <- unlist(dfNCBI$GenBank_Accession) 

# Check class. Ensure that the class is a character vector.

class(Accessions) 

# The following feature data can also be downloaded from NCBI by searching "Nicotiana benthamiana"[Organism] AND "CDS"[FKEY]" on https://www.ncbi.nlm.nih.gov and then clicking on the "Send to:" tab -> "Complete record" -> "File" -> "Format: Feature table" -> "Default order" -> "Create file".
# Using the getAnnotationsGenBank() function from the ape package to get the annotations for each sequence from GenBank. 
# The following line of code was adapted from the APE package documentation available from: https://cran.r-project.org/web/packages/ape/ape.pdf. (Paradis, 2022).

Feature_List <- getAnnotationsGenBank(Accessions, quiet = FALSE) 

# This function takes approximately 9.082817 minutes to run. For a faster script time, the feature list can also be loaded into R from the file using the following lines of code:
# Other options include the biofiles package and the genbankr package functions. However, the ape package function has a faster run time that can also handle large amounts of queries. It also does not require parsing of GenBank files (which may be tedious when dealing with a large number of records).

# Convert the list of dataframes to a data frame. 
# The following line of code was adapted from: https://stackoverflow.com/questions/2851327/combine-a-list-of-data-frames-into-one-data-frame-by-row. (Klieg, 2018).

dfFeatures <- bind_rows(Feature_List, .id = "GB_Accession") 

# Change the column names to be more descriptive.

colnames(dfFeatures) <- c("GB_Accession","Start_of_CDS","End_of_CDS","Type","Product","Protein_ID","Gene")

# Remove all objects no longer needed from the environment.

rm(Accessions,Feature_List)

#### 07 DATA AQUISITION : OBTAIN GENE NAMES ####

# The gene name for each entry is distributed over three columns: Product, Protein_ID, and Gene. To ensure that these are formatted correctly, extract the relevant information from each of the columns and merge.

# First, for entries that have a listed gene name, extract the gene names into a character vector.

Gene_Name_List_1 <- dfFeatures$Gene 

# Second, for entries that have a listed gene name in the protein description, extract the gene names into a character vector.

Gene_Name_List_2 <- dfFeatures$Protein_ID %>%
  str_extract("prot_desc.*") %>%
  str_extract(".*;") %>%
  str_replace("prot_desc:", "") %>%
  str_replace(";.*", "") %>%
  str_replace("autophagy", "") %>%
  str_squish()

# Third, for entries that have a listed gene name in the notes, extract the gene names into a character vector.

Gene_Name_List_3 <- dfFeatures$Protein_ID %>%
  str_extract("note.*") %>%
  str_replace("note:", "") %>%
  str_replace(";.*", "") %>%
  str_squish()

# Condense these four lists, by order of importance (i.e. as aforementioned above)
  
All_Gene_Names <- coalesce(x = Gene_Name_List_1, y = Gene_Name_List_2) %>%
  coalesce(y = Gene_Name_List_3) 

# Add back to the data frame.

dfFeatures["Gene"] <- All_Gene_Names

# Edit the protein id column to neatly display only the GenBank Protein IDs.

New_Protein_IDs <- str_match(dfFeatures$Protein_ID, ".*gb\\|([^|]+)")

dfFeatures["Protein_ID"] <- New_Protein_IDs[,2]

# Export these two data frames as tables to .csv format files (if desired). 

write.table(x = dfFeatures, 
            file = "Features_Data.csv")

write.table(x = dfNCBI, 
            file = "NCBI_Data.csv")

# Remove all objects no longer needed from the environment.

rm(Gene_Name_List_1,Gene_Name_List_2,Gene_Name_List_3,All_Gene_Names,New_Protein_IDs)

#### 08 DATA SUMMARY : NCBI INFORMATION ####

# Check the type of class and summary information.

class(dfNCBI)
class(dfFeatures)

# Both of these are data frames.

# Check the type of information in each data frame.

summary(dfNCBI)
summary(dfFeatures)

# For the dfNCBI, all are character data except for one column specifiying the sequence lengths.
# For the dfFeatures, all are character data except for two columns specifiying the sequence start and stop positions.

# Any missing sequence entries?

sum(is.na(dfNCBI$Untrimmed_Sequences)) # 0

# Any missing gene name entries?

sum(is.na(dfFeatures$Gene)) # 535 missing gene names.

# Any missing protein name entries?

sum(is.na(dfNCBI$Titles)) # 0 

# Check if the number of entries are the same.

length(unique(dfNCBI$GenBank_Accession)) == length(unique(dfFeatures$GB_Accession))

# TRUE. Yes, these are the same.

#### 09 DATA FILTERING : COMBINE INTO ONE DATAFRAME ####

# The script now has two data frames as follows:

# 1) dfFeatures - A dataframe holding several features annotated for each sequence.
# 2) dfNCBI - A dataframe holding the untrimmed sequences and version information.

# The number of features can exceed the number of records.

table(duplicated(dfFeatures$GB_Accession)) # FALSE = 332 "extra" features annotated.

# What type of features have been annotated in the records?

table(dfFeatures$Type)

# Several features including coding sequences, mRNA 3' and 5' untranslated regions (UTRs).
# Filter the data frame to retain features marked "CDS", eliminating UTRs, rRNAs, and other non-coding regions.
# Retain only feature "rows" marked as CDS.

dfFeatures.CDS <- dfFeatures %>%
  filter(Type == "CDS")

# Subsetting dfData by the column names in dfFeatures. This will retain only records marked as "CDS".
  
dfCombined <- subset(dfNCBI, GenBank_Accession %in% dfFeatures.CDS$GB_Accession)

# Add the feature columns to the new data frame.

dfData <- cbind(dfCombined, dfFeatures.CDS[,2:7])

# Append a new column with the length of the CDS.

dfData["GenBank_Length_of_CDS"] <- ((dfData$End_of_CDS - dfData$Start_of_CDS) + 1)

# Rearrange column names and edit row names for readability.

dfData <- dfData[,c("GenBank_Accession", "NCBI_ID", "Titles", "Type", "Untrimmed_Sequence_Length", "Untrimmed_Sequences", "Start_of_CDS", "End_of_CDS", "GenBank_Length_of_CDS", "Gene", "Protein_ID", "Product")]

rownames(dfData) <- NULL

# Remove all objects no longer needed from the environment.

rm(dfNCBI, dfFeatures, dfFeatures.sub, dfCombined)

#### 10 DATA FILTERING : FILTER RECORDS ####

# Next, filter the data frame to ensure that only complete coding sequences are retained.
# This is done because the empirical codon counts for each entry necessitates complete sequences.
# If partial sequences are retained, the incomplete codon counts for each sequence would skew the overall frequencies per codon.

Filtered_Titles <- dfData$Titles %>%
  str_replace("Nicotiana benthamiana mRNA for hypothetical protein.*", "") %>% # Remove records for hypothetical proteins
  str_replace("Nicotiana benthamiana clone.*", "") %>% # Remove records for any N. benthamiana clones
  str_extract("Nicotiana benthamiana.*") %>% # Ensure that only N. benthamiana entries are retained
  str_extract(".*mRNA, complete cds.*")  # Keep only complete cds (i.e. remove partial cds) 

# Remove any coding sequences for hypothetical proteins.
# Subsetting dfData by the filtered protein names in Filtered_Titles.

dfData.sub <- subset(dfData, Titles %in% Filtered_Titles)

# How many entries (excluding any duplicates)?

length(unique(dfData.sub$Titles)) # 281 Entries.

# Remove all objects no longer needed from the environment.

rm(Filtered_Titles, dfData)

#### 11 DATA FILTERING : REMOVE DUPLICATED TITLES ####

# Find duplicated entries.

Names <- unlist(dfData.sub$Titles)

# The following line of code was adapted from: https://stackoverflow.com/questions/16905425/find-duplicate-values-in-r. (David, 2016).

Duplicate_Names <- Names[ Names %in% Names[duplicated(Names)] ]

# The following lines of code were adapted from: https://stackoverflow.com/questions/57832161/how-to-remove-na-in-character-vector-in-r. (Laurent, 2019).

NAs <- is.na(Duplicate_Names) 

All_Duplicates <- Duplicate_Names[!NAs] 

# Final list of gene names that have multiple entries.

Unique_Duplicates <- unique(All_Duplicates)

#rm(Duplicate_Names,NAs, dup)

length(Unique_Duplicates)

# 2 records with multiple entries. Can these be filtered further?
# It would be constructive to ensure that these associated sequences and information makes sense.

# Subset duplicates into a separate data frame.

dfData_Duplicates <- data.frame()

for (i in Unique_Duplicates) {
  matching_row <- dfData.sub %>%
    filter(Titles == i)
  dfData_Duplicates <- rbind(dfData_Duplicates, matching_row)
}

# Filtering of duplicates. 
# There are several duplicates for a single gene in dfData_Duplicates.
# Retain the longest sequence. If there are flanking non-coding regions in the sequence, it can be trimmed below.
# In this case, only the longest sequence will be retained since there may be more information present that is valuable to downstream analysis. 

# Steps:
# 1. Subset the duplicates from dfData
# 2. Remove these duplicated entries from dfData
# 3. Sort the duplicated entries by sequence length.
# Remove duplicates using function. This will retain the first entry (a.k.a. the longest sequence) and remove the succeeding entries.

# Sort dfData_Duplicates by sequence length. 

dfDataSorted <- dfData_Duplicates[order(dfData_Duplicates$Untrimmed_Sequence_Length, decreasing = TRUE),] 

# The following lines of code were fixed using the advice from: https://stackoverflow.com/questions/62075537/r-pipe-in-does-not-work-with-stringrs-str-extract-all. (Kirshna, A. "Akrun" 2020).

terms <- duplicated(dfDataSorted$Titles) %>%
  {str_replace(.,"FALSE","No") %>%
  str_replace(.,"TRUE","Yes")}

dfDataSorted["Is.A.Duplicate"] <- terms
rm(terms)

# Subset the column for unique entries.

dfKeep <- subset(dfDataSorted, Is.A.Duplicate == "No") 

# Remove the last column indication is the entry is a duplicate or not.

dfKeep <- dfKeep[1:(length(dfKeep)-1)] 

rm(dfDataSorted)

# The following line of code was adapted from: https://stackoverflow.com/questions/17338411/delete-rows-that-exist-in-another-data-frame. (Gillespie, 2013).
# Remove all the duplicates from the complete data frame.

dfNew <- dfData.sub[!(dfData.sub$NCBI_ID %in% dfData_Duplicates$NCBI_ID),] 

# Add back only the filtered non-duplicated enties that are in dfKeep

dfFinal <- rbind(dfNew,dfKeep) 

# Confirm we have unique protein names.

table(duplicated(dfFinal$Titles))

# Remove all objects no longer needed from the environment.

rm(dfData_Duplicates,dfKeep,dfNew)

#### 12 DATA FILTERING : TRIM SEQUENCES ####

# To get rid of flanking non-coding regions in the CDS, need to trim the sequences based on the start, stop, and length columns.
# The idea to use mapply function came from: https://stackoverflow.com/questions/6827299/r-apply-function-with-multiple-parameters. (Borochkin, A. A. "Alexander" 2016).

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

rm(Trimmed_Sequences,Trimmed_Lengths)


#### 13 GENERATE CODON USAGE TABLE ####

# https://bioconductor.riken.jp/packages/3.8/bioc/vignettes/coRdon/inst/doc/coRdon.html

# First, store all the coding sequences as a DNAStringSet object.

Coding_Seqs <- DNAStringSet(dfFinal$Trimmed_CDS)

# Then, using the codonCounts() and codonTable() functions from the codRon package create a matrix to view the counts per codon.

Updated_CU <- codonCounts(codonTable(Coding_Seqs))

# Third, calculate the total number of counts per column.

Number <- colSums(Updated_CU)

# Next, calculate the frequencies of each codon (out of 1000) (i.e. how the Kazuza data base reports the frequencies)

`X.1000` <- round(x = ((Number/sum(Updated_CU))*1000), 
                  digits = 2)

# Lastly, create two character vectors for the amino acid and corresponding codons. The amino acids are assigned with the assumtion that the codons are listed from AAA to TTT as defined in the codonTable function.

AmAcid <- c("Lys","Asn","Lys","Asn","Thr","Thr","Thr","Thr","Arg","Ser","Arg","Ser","Ile","Ile","Met","Ile","Gln","His","Gln","His","Pro","Pro","Pro","Pro","Arg","Arg","Arg","Arg","Leu","Leu","Leu","Leu","Glu","Asp","Glu","Asp","Ala","Ala","Ala","Ala","Gly","Gly","Gly","Gly","Val","Val","Val","Val","End","Tyr","End","Tyr","Ser","Ser","Ser","Ser","End","Cys","Trp","Cys","Leu","Phe","Leu","Phe")

Single_Letter_Abbreviation <- c("K","N","K","N","T","T","T","T","R","S","R","S","I","I","M","I","Q","H","Q","H","P","P","P","P","R","R","R","R","L","L","L","L","E","D","E","D","A","A","A","A","G","G","G","G","V","V","V","V","X","Y","X","Y","S","S","S","S","X","C","W","C","L","F","L","F")

Codon <- colnames(Updated_CU)

# Finally, create a data frame using the four vectors created above.
  
dfCodon_Usage_Table <- data.frame(AmAcid, Single_Letter_Abbreviation, Codon, Number, `X.1000`)

row.names(dfCodon_Usage_Table) <- NULL

#### 14 EXPORT CODON USAGE AND ADDITIONAL DATA ####

# Write the updated codon usage data frame to a text file.

write.table(x = dfCodon_Usage_Table, 
            file = "Updated_Codon_Usage.txt",
            quote = FALSE,
            row.names = FALSE)

# Write the updated codon usage data frame to a .csv format file.

write_csv(x = dfCodon_Usage_Table, 
          file = "Updated_Codon_Usage.csv")

# Write the final data frame to a .csv format file.

write_csv(x = dfFinal, 
          file = "Updated_Codon_Usage_Information.csv")

#### 15 REFERENCES ####

# Winter, D. (2020, November 11). Rentrez Tutorial. https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html
# Klieg, J. (2018, February 27). Combine a list of data frames into one data frame by row. Stack Overflow. https://stackoverflow.com/questions/2851327/combine-a-list-of-data-frames-into-one-data-frame-by-row
# Paradis, E. (2022, March 2). Package ‘ape’. Comprehensive R Archive Network (CRAN). https://cran.r-project.org/web/packages/ape/ape.pdf
# David, C. (2016, October 1). Find duplicate values in R [duplicate]. Stack Overflow. https://stackoverflow.com/questions/16905425/find-duplicate-values-in-r
# Laurent, S. (2019, September 7). How to remove NA in character vector in R [duplicate]. Stack Overflow. https://stackoverflow.com/questions/57832161/how-to-remove-na-in-character-vector-in-r
# Kirshna, A. "Akrun" (2020, May 28). R pipe in does not work with stringR’s str_extract_all(). Stack Overflow. https://stackoverflow.com/questions/62075537/r-pipe-in-does-not-work-with-stringrs-str-extract-all
# Gillespie, C. S. (2013, June 27). Delete rows that exist in another data frame? [duplicate]. Stack Overflow. https://stackoverflow.com/questions/17338411/delete-rows-that-exist-in-another-data-frame
# Borochkin, A. A. "Alexander" (2016, November 1). R apply function with multiple parameters. Stack Overflow. https://stackoverflow.com/questions/6827299/r-apply-function-with-multiple-parameters


#