# Importing and Cleaning UniProt Data
# 9 June 2022
# Anchitaa Ghag

#### 01 INSTALL PACKAGES & DOWNLOAD LIBRARIES ####

# install.packages("seqinr")
library("seqinr")
# install.packages("tidyverse")
library("tidyverse")
# install.packages("rentrez")
library("rentrez")
# install.packages("XML")
library("XML")

#### 02 IMPORT EXCEL FILE ####

# A list of the proteins in Nicotiana benthamiana, protein sequences, and corresponding information can be obtained from the UniProt database (Bateman et. al., 2021).
# https://www.uniprot.org/uniprot/?query=Nicotiana%20benthamiana&columns=id%2Centry%20name%2Creviewed%2Cprotein%20names%2Cgenes%2Corganism%2Clength%2Cfragment%2Ccitation%2Csequence&sort=score

# Load the excel file with protein sequences downloaded from UniProt.

dfUniProt <- as.data.frame(readxl::read_xlsx("~/Major_Research_project_2022/06_Code/UniProt_Data.xlsx"))

# Reformat column names to ensure there are no spaces.

names(dfUniProt) <- make.names(names(dfUniProt),unique(TRUE))

# View the data frame.

View(dfUniProt)

# Not all organisms are N. benthamiana. Subset only the N. benthamiana entries.

dfData <- dfUniProt[dfUniProt$Organism == "Nicotiana benthamiana",]

# Reformat column names to ensure there are no spaces.

names(dfData) <- make.names(names(dfData),unique(TRUE))

#### 03 SUMMARY STATS ####

# Check the type of class and summary information.

class(dfData)

# All character data except for one column specifiying the protein lengths.

summary(dfData)

# How many entries have been reviewed?

length(dfData$Gene.names) - sum(dfData$Status == "unreviewed") 

# Only 16 entries have been reviewed.

# Any missing sequence entries?

sum(is.na(dfData$Sequence)) # 0

# Any missing gene name entries?

sum(is.na(dfData$Gene.names)) # 304 missing gene names.

# Any missing protein name entries?

sum(is.na(dfData$Protein.names)) # 0 

# Since there are missing gene names the data set must be split. 
# Next, we will look at the gene names for the entries that have them.

#### 04 CLEANING UP GENE NAMES ####

# View the gene names.

dfData$Gene.names

# There are some gene names that have been marked with "Nb" at the beginning. 
# Filter this out since we have already subsetted the data above to include only N. benthamiana entries.

Replaced_Gene_Names <- str_trim(str_replace(dfData$Gene.names,"^Nb"," ")) %>% 

# There is one entry for "GIP2 Niben101Scf03191g04002" which needs to be changed to "GIP2".
  {str_trim(str_replace(Replaced_Gene_Names,"GIP2 Niben101Scf03191g04002","GIP2"))}

dfData["Gene.names"] <- Replaced_Gene_Names

rm(Replaced_Gene_Names)

#### 05 DEALING WITH DUPLICATES ####

# Any duplicate gene names?

length(dfData$Gene.names) - length(unique(dfData$Gene.names)) # 340 duplicate entries.

# There appears to be 526 unique gene names.

# View the gene names.

unique(dfData$Gene.names)

# There are several entries of the LecRK gene name because it contains an alpha-numeric string at the end.
# Change the name of these entries using str_replace() function & view again.

unique(str_replace(unique(dfData$Gene.names),"LecRK.*","LecRK")) 

# There are 490 unique gene names.

# There are several LecRK entries corresponding to different clades. Filter these out separately into a data frame. 
# Treat these entries separeately when importing the CDS from NCBI.

dfData["Gene.names"] <- str_trim(str_replace(dfData$Gene.names,"LecRK.*","LecRK")) 

dfLecRK <- subset(dfData, Gene.names == "LecRK") # (Kincaid, 2011). Code adapted from: https://stackoverflow.com/questions/7381455/filtering-a-data-frame-by-values-in-a-column

# Remove LecRK entries from dfData data frame.

dfData_No_LecRK <- dfData[!(dfData$Protein.names %in% dfLecRK$Protein.names),]
# https://stackoverflow.com/questions/16905425/find-duplicate-values-in-r

rm(dfData)

# Find duplicated entries.

Gene_Names <- unlist(dfData_No_LecRK$Gene.names)

Duplicate_Names <- Gene_Names[ Gene_Names %in% Gene_Names[duplicated(Gene_Names)] ]

# https://stackoverflow.com/questions/16905425/find-duplicate-values-in-r

NAs <- is.na(Duplicate_Names) 
dup <- Duplicate_Names[!NAs] 
# https://stackoverflow.com/questions/57832161/how-to-remove-na-in-character-vector-in-r

# Final list of gene names that have multiple entries.

unique_duplicates <- unique(dup)

rm(Duplicate_Names,NAs, dup)

# 22 genes with multiple entries. Can these be filtered further?
# It would be constrictive to ensure that these associated sequences and information makes sense.

# Subset duplicates into a separate data frame.

dfData_Duplicates <- data.frame()

for (i in unique_duplicates) {
  matching_row <- dfData_No_LecRK %>%
    filter(Gene.names == i)
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

dfDataSorted <- dfData_Duplicates[order(-dfData_Duplicates$Length),] 

#https://stackoverflow.com/questions/62075537/r-pipe-in-does-not-work-with-stringrs-str-extract-all

terms <- duplicated(dfDataSorted$Gene.names) %>%
{str_replace(.,"FALSE","No") %>%
str_replace(.,"TRUE","Yes")}

dfDataSorted["Is.A.Duplicate"] <- terms
rm(terms)

dfKeep <- subset(dfDataSorted, Is.A.Duplicate == "No")
dfKeep <- dfKeep[1:(length(dfKeep)-1)]

rm(dfDataSorted)

# Add back to original data frame.
# https://stackoverflow.com/questions/17338411/delete-rows-that-exist-in-another-data-frame

dfNew <- dfData_No_LecRK[!(dfData_No_LecRK$Protein.names %in% dfData_Duplicates$Protein.names),]

df_Final <- rbind(dfNew,dfKeep)

rm(Gene_Names,dfData_Duplicates,dfData_No_LecRK,dfKeep,dfNew)

#### 06 EXAMINING THE PROTEIN LENGTHS ####

# https://r-charts.com/distribution/add-points-boxplot/

# Subset the potential list of entries for BLAST from dfData.

dfData_Selected_For_BLAST <- df_Final %>%
  filter(is.na(Gene.names))

# Are the protein lengths normally distributed?
# Alternatively, Are there any peptides present or abnormalities in lengths?

hist(dfData_Selected_For_BLAST$Length)
shapiro.test(dfData_Selected_For_BLAST$Length)

# W = 0.87271, p-value = 4.144e-15
# Data is not normally distributed.

boxplot(dfData_Selected_For_BLAST$Length,
        col = "white",
        ylab = "Length of Sequence (# of amino acids)", 
        xlab = "Protein Sequences from Uniprot", 
        main = "Boxplot of Protein Sequence Length (no gene information)")

length(boxplot(dfData_Selected_For_BLAST$Length)$out) # 11 Outliers.

#https://www.r-statistics.com/2011/01/how-to-label-all-the-outliers-in-a-boxplot/
#install.packages("TeachingDemos")
library("TeachingDemos")
source("https://raw.githubusercontent.com/talgalili/R-code-snippets/master/boxplot.with.outlier.label.r") # Load the function

y <- dfData_Selected_For_BLAST$Length
lab_y <- dfData_Selected_For_BLAST$Protein.names

boxplot.with.outlier.label(y, lab_y, spread_text = F)
set.seed(1212)
stripchart(dfData_Selected_For_BLAST$Length,              # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 4,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over

# Retain these outliers since there is no evidence to exclude them at this point.

rm(lab_y,y)

#### 07 SUBSET DATA FRAME FOR NEXT STEPS ####

# Subset data frame into two separate data frames depending on gene information presence or absence.

# For proteins with corresponding gene information.

dfData_For_NCBI <- df_Final %>%
  filter(!is.na(Gene.names))

# For proteins without corresponding gene information.

dfData_For_BLAST <- df_Final %>%
  filter(is.na(Gene.names))

# Write both of these newly created data frames to separate .txt files. Saving to a .txt file makes it easier to import in non-R environments (i.e. to open in Excel)
# This line of code was adapted from: https://stackoverflow.com/questions/18514694/how-to-save-a-data-frame-in-a-txt-or-excel-file-separated-by-columns

write.table(dfData_For_NCBI,"Data_For_NCBI.txt",sep="\t",row.names=FALSE)
write.table(dfLecRK,"Data_For_LecRK.txt",sep="\t",row.names=FALSE)
write.table(dfData_For_BLAST,"Data_For_BLAST.txt",sep="\t",row.names=FALSE)


#### 08 REFERENCES ####

# Bateman, A., Martin, M.-J., Orchard, S., Magrane, M., Agivetova, R., Ahmad, S., Alpi, E., Bowler-Barnett, E. H., Britto, R., Bursteinas, B., Bye-A-Jee, H., Coetzee, R., Cukura, A., da Silva, A., Denny, P., Dogan, T., Ebenezer, T., Fan, J., Castro, L. G., … Teodoro, D. (2021). UniProt: the universal protein knowledgebase in 2021. Nucleic Acids Research, 49(D1), D480–D489. https://doi.org/10.1093/nar/gkaa1100/

# Kincaid, D. (2011, September 11). Filtering a data frame by values in a column [duplicate]. Stack Overflow. https://stackoverflow.com/questions/7381455/filtering-a-data-frame-by-values-in-a-column</div>