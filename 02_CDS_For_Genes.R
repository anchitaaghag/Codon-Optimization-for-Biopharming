# Importing CDS Sequences for Proteins with Gene Information
# 9 June 2022
# Anchitaa Ghag

#### POTENTIAL CHANGES TO BE MADE ####

# 1) It currently takes around 17 minutes to run the For Loop below. This is too slow. 
# Potentially, this could be converted to a function and then use lapply or sapply to reiterate over all the gene names.

#### 01 INSTALL PACKAGES & DOWNLOAD LIBRARIES ####

# install.packages("seqinr")
library("seqinr")
# install.packages("tidyverse")
library("tidyverse")
# install.packages("rentrez")
library("rentrez")
# install.packages("XML")
library("XML")
#BiocManager::install("genbankr")
library("genbankr")


#### 02 IMPORT FILE ####

# Load the .txt file with gene names.

dfData_For_NCBI <- (read.table("~/Major_Research_project_2022/06_Code/Codon-Optimization-for-Biopharming/Data_For_NCBI.txt"))
dfLecRK <- (read.table("~/Major_Research_project_2022/06_Code/Codon-Optimization-for-Biopharming/Data_For_LecRK.txt"))

# Check the type of class and summary information.

class(dfData_For_NCBI)
summary(dfData_For_NCBI)

class(dfLecRK)
summary(dfLecRK)

# Update the column names.

column_names_1 <- dfData_For_NCBI[1,]
column_names_2 <- dfLecRK[1,]

colnames(dfData_For_NCBI) <- column_names_1
colnames(dfLecRK) <- column_names_2

dfData_For_NCBI <- dfData_For_NCBI[-1,]
dfLecRK <- dfLecRK[-1,]

rm(column_names_1, column_names_2)

# Ensure there are no duplicates.

length(unique(dfData_For_NCBI$Gene.name)) # 480
length(dfData_For_NCBI$Gene.names) # 480

#### 03 CREATE LIST OF GENES ####

Gene_Names <- dfData_For_NCBI$Gene.names

# 480 genes are too large even for a web history object. In addition, one gene search may have multiple entries.
# Thus, we will run the search for one gene at a time.

# Finding the corresponding (CDS) coding sequence using rentrez package, if possible.
# If corresponding CDS are not found, can pass on sequences to BLAST.

# Get a list of all the databases available. 

entrez_dbs(config = NULL)

# Will be using the "nucleotide" database to search for CDS.

# Reformat all gene names to the search term style of entrez_search() function.
# e.g. "Nicotiana benthamiana[Organism] AND CPP1[Gene]"

Reformat <- function(text, wrapper1, wrapper2) {
  
  solution <- c(wrapper1,text,wrapper2)
  results <- paste(solution, collapse = " ")
}

#### 04 For Loop to Find Entries for Each Gene ####

# https://www.dataquest.io/blog/control-structures-in-r-using-loops-and-if-else-statements/

# Create an empty data frame.
mydf = data.frame()
# Create an empty list to hold the names of gene which has no hits.
mylist = list()

system.time (
# For Loop through each gene name in the gene name list.
for (i in Gene_Names) {
  # Create the search term.
  wrapping_txt_begin <- "Nicotiana benthamiana[Organism] AND" 
  wrapping_txt_end <- "[Gene]"
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

rm(checkpoint,df_len,titles,Titles,nico_summs,nico_search,search_term,wrapped_str,wrapping_txt_begin,wrapping_txt_end,dfIDs_and_Titles)

# Write this data to a file.
# write.table(mydf,"Results_From_NCBI.txt",sep="\t",row.names=FALSE)

View(mydf)
View(mylist)

# Filter the data by information in the title.

Titles_Filtered <- mydf$Titles %>%
  str_extract("Nicotiana benthamiana.*") %>% # Keep only N. benthamiana
  str_extract(".*complete cds") # Keep only complete cds (i.e. remove partial cds)

mydf["Filtered_Title"] <- Titles_Filtered

# Remove NAs

mydf <- na.omit(mydf)

# How many genes have multiple entries?

nrow(mydf) - length(unique(mydf$Gene.Name)) #16

#### 05 HANDLING MULTIPLE ENTRIES ####

# Because entrez-searches are not case sensetive, many different enties are possible for the same gene.
# For example, while AGO1b and AGO1B are the same gene, their diffrent spelling leads to 4 entries returned instead of one.

# Find duplicated entries.

Gene_Names <- unlist(mydf$Gene.Name)

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
    filter(Gene.Name == i)
  dfData_Duplicates <- rbind(dfData_Duplicates, matching_row)
}

rm(unique_duplicates, matching_row, i)
View(dfData_Duplicates)

# Filtering of duplicates. 
# There are several duplicates for a single gene in dfData_Duplicates.
# For a number of these genes, there is one mRNA sequence and several "duplicates" gene sequences. e.g. psbS1
# Need to retain the mRNA sequnce and not complete gene information.
# In addition, the full sequence will most likely begin with a "M" rather than another amino acid.
# For other entries, there are some "precursor" or "clone" entries. e.g. PCNA has 2 clone sequences and one complete CDS sequence.
# In this case, only the complete CDS sequence will be retained since there may be more information present that is valuable to downstream analysis. 

# Clean-up dfData_Duplicates gene names.

new_gene_names <- str_replace(dfData_Duplicates$Gene.Name,"AGO1a", "AGO1A") %>%
  {str_replace(.,"AGO1b", "AGO1B")} 

dfData_Duplicates["Gene.Name"] <- new_gene_names
rm(new_gene_names)

# Clean-up dfData_Duplicates titles.

new_titles <- str_replace(dfData_Duplicates$Titles,".*clone.*", "NA") %>%
  {str_replace(.,".*gene, complete .*", "NA") %>%
  str_replace(.,".*AGO1a.*", "NA") %>%
  str_replace(.,".*AGO1b.*", "NA")} 

dfData_Duplicates["Titles"] <- new_titles
rm(new_titles)

dfData_Duplicates <- dfData_Duplicates %>%
  filter(!Titles == "NA")

#https://stackoverflow.com/questions/62075537/r-pipe-in-does-not-work-with-stringrs-str-extract-all

terms <- duplicated(dfData_Duplicates$Gene.Name) %>%
  {str_replace(.,"FALSE", "No") %>%
  str_replace(.,"TRUE", "Yes")}

dfData_Duplicates["Is.A.Duplicate"] <- terms
rm(terms)

dfKeep <- subset(dfData_Duplicates, Is.A.Duplicate == "No")
dfKeep <- dfKeep[1:(length(dfKeep)-1)]

# Add back to original data frame.
# https://stackoverflow.com/questions/17338411/delete-rows-that-exist-in-another-data-frame

dfNew <- mydf[!(mydf$Titles %in% dfData_Duplicates$Titles),]

df_Final <- rbind(dfNew,dfKeep)

rm(Gene_Names,dfData_Duplicates,dfKeep,dfNew)

#### IDs to CDS Information ####

source("~/Major_Research_project_2022/06_Code/Codon-Optimization-for-Biopharming/Functions.R")

# Convert the IDs to a list.
ids <- df_Final$IDs
# Apply the Get_Accession_Number() function to each element in the list.
system.time(Accession <- lapply(ids, Get_Accession_Number))
# Append the corresponding GenBank accession numbers to the final data frame.
df_Final["Accession.Number"] <- unlist(Accession)
rm(ids)

# Repeat this process for obtaining the start, end, and length of the coding sequences.
# Convert the GenBank accession numbers to a list.
gba_acc <- df_Final$Accession.Number
# Apply the Get_CDS ... () functions to each element in the list.
system.time(Start <- lapply(gba_acc, Get_CDS_Start))
system.time(End <- lapply(gba_acc, Get_CDS_End))
system.time(Length <- lapply(gba_acc, Get_CDS_Length))
# Append the corresponding information to the final data frame.
df_Final["Start.of.CDS"] <- unlist(Start)
df_Final["End.of.CDS"] <- unlist(End)
df_Final["Length.of.CDS"] <- unlist(Length)
rm(gba_acc,Start,End,Length)

#### IDs to CDS Sequences #### 

# Use ids to obtain fasta file of coding sequences.

#ID Name.of.Protein Start Stop Length CDS.Sequence

ids = c(df_Final$IDs)
nico_retrive <- entrez_fetch(db="nucleotide", id = ids[1:300], rettype = "fasta")
nico_retrive1 <- entrez_fetch(db="nucleotide", id = ids[301:356], rettype = "fasta")
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
df_Final["CDS.(Untrimmed)"] <- sequences

### TIMMING THE SEQUENCES ####

# To get rid of flanking non-coding regions in the CDS, need to trim the seqences based on the start, stop, and length columns.

start_positions <- df_Final$Start.of.CDS
end_positions <- df_Final$End.of.CDS

# https://stackoverflow.com/questions/6827299/r-apply-function-with-multiple-parameters
trimmed_seqs <- mapply(str_sub,string=sequences,start=start_positions,end=end_positions)

# Test if this has been done successfully by comparing the lengths of the trimmed sequences with the expected lengths (i.e. the Length.of.CDS column in the data frame.)

trimmed_lengths <- unlist(lapply(trimmed_seqs,nchar))
df_Final["Trimmed.Length"] <- trimmed_lengths
identical(df_Final$Length.of.CDS,df_Final$Trimmed.Length)

# TRUE. 
# Can proceed and add the trimmed CDS to the data frame.

df_Final["Trimmed.CDS"] <- trimmed_seqs

# Recheck the timmed sequences. Are these divisible by three?
# If sequence lengths are divisible by three, we can continue to retain. 
# If not, either check the trimming or discard.

# FIXME


#### NO HITS FOUND: DATA FOR BLAST####

# For the genes where no corresponding CDS information was found in Nicotiana benthamiana.

# This will be merged with the dfData_For_BLAST data frame.
# Next, this data set will be run through blast.
# If no hits are found by protein name then:
# Use N. tabacum or use homologs

#### CU Table Generation ####

#https://bioconductor.riken.jp/packages/3.8/bioc/vignettes/coRdon/inst/doc/coRdon.html

#BiocManager::install("coRdon")
library("coRdon")

cds <- DNAStringSet(df_Final$Trimmed.CDS)
CUTs <- codonTable(cds)
CU <- codonCounts(CUTs)


row.names(CU) <- df_Final$Titles

#### Statistical Tests ####

milc <- MILC(CUTs)
head(milc)
  
  
