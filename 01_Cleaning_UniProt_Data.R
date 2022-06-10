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

# Load the excel file with protein sequences downloaded from UniProt.
# FIXME How can we get this data? Is there a package? Can I download through a link??

dfData <- as.data.frame(readxl::read_xlsx("~/Major_Research_project_2022/06_Code/uniprot-nicotiana+benthamiana.xlsx"))

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

# Since there are missing gene names the data set must be split. 
# Next, we will look at the gene names for the entries that have them.

#### 03 DEALING WITH DUPLICATES ####

# Any duplicate gene names?

length(dfData$Gene.names) - length(unique(dfData$Gene.names)) # 339 duplicate entries.

# There appears to be 527 unique gene names.

# View the duplicated names.

unique(dfData$Gene.names)

# There are several entries of the LecRK gene name because it contains an alpha-numeric string at the end.
# Change the name of these entries using str_replace() function & view again.

unique(str_replace(unique(dfData$Gene.names),"LecRK.*","LecRK")) 

Replaced_Gene_Names <- str_replace(dfData$Gene.names,"LecRK.*","LecRK")
dfData["Gene.names"] <- Replaced_Gene_Names

# There are 491 unique gene names.

Duplicate_Names <- unlist(dfData$Gene.names)

duplicates <- Duplicate_Names[ Duplicate_Names %in% Duplicate_Names[duplicated(Duplicate_Names)] ]

# https://stackoverflow.com/questions/16905425/find-duplicate-values-in-r

NAs <- is.na(duplicates) 
dup <- duplicates[!NAs] 
# https://stackoverflow.com/questions/57832161/how-to-remove-na-in-character-vector-in-r

# Final list of gene names that have multiple entries.

unique_duplicates <- unique(dup)

# 22 genes with multiple entries. Can these be filtered further?

# Since these specific protein sequences are not going to be used, keep only unique gene name entries going forward.
# However, it would be constrictive to ensure that these associated sequences and information makes sense.

# Subset duplicates into a separate data frame.

dfData_Duplicates <- data.frame()

for (i in unique_duplicates) {
  matching_row <- dfData %>%
    filter(Gene.names == i)
  dfData_Duplicates <- rbind(dfData_Duplicates, matching_row)
}

View(dfData_Duplicates)

# Manually filter out entries and fragments.
# remove marked (frafment) and keep the longest of the sequences if tie, pick last one
rows_to_keep <- c(1,5,13,14,16,20,22,62,64,65,71,73,75,77,79,81,83,86,88,89,90,92,93)
dfKeep <- dfData_Duplicates[rows_to_keep,]

# Add back to original data frame.

dfNew <- dfData[!(dfData$Protein.names %in% dfData_Duplicates$Protein.names),]

df_Final <- rbind(dfNew,dfKeep)

rm(dup,Duplicate_Names,duplicates)

#### 05 EXAMINING THE PROTEIN LENGTHS ####

# https://r-charts.com/distribution/add-points-boxplot/

# Subset the potential list of entries for BLAST from dfData.

dfData_Selected_For_BLAST <- dfData %>%
  filter(is.na(Gene.names))

# Are the protein lengths normally distributed?
# Alternatively, Are there any peptides present or abnormalities in lengths?

hist(dfData_Selected_For_BLAST$Length)
shapiro.test(dfData_Selected_For_BLAST$Length)

# W = 0.8731, p-value = 3.861e-15
# Data is not normally distributed.

boxplot(dfData_Selected_For_BLAST$Length,
        col = "white",
        ylab = "Length of Sequence (# of amino acids)", 
        xlab = "Protein Sequences from Uniprot", 
        main = "Boxplot of Protein Sequence Length (no gene information)")

length(boxplot(dfData_Selected_For_BLAST$Length)$out) # 9 Outliers.

#https://www.r-statistics.com/2011/01/how-to-label-all-the-outliers-in-a-boxplot/
#install.packages("TeachingDemos")
library("TeachingDemos")
source("https://raw.githubusercontent.com/talgalili/R-code-snippets/master/boxplot.with.outlier.label.r") # Load the function
set.seed(6484)
y <- dfData_Selected_For_BLAST$Length
lab_y <- dfData_Selected_For_BLAST$Protein.names

boxplot.with.outlier.label(y, lab_y, spread_text = F)

stripchart(dfData$Length,              # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 4,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)        # Add it over


# FIXME Interpretation.


#### 06 SUBSET DATA FRAME FOR NEXT STEPS ####

# Subset data frame into two separate data frames depending on gene information presence or absence.

# For proteins with corresponding gene information.

dfData_For_NCBI <- unique(df_Final$Gene.names)

# For proteins without corresponding gene information.

dfData_For_BLAST <- dfData %>%
  filter(is.na(Gene.names))

# Write both of these newly created data frames to separate .txt files. Saving to a .txt file makes it easier to import in non-R environments (i.e. to open in Excel)
# This line of code was adapted from: https://stackoverflow.com/questions/18514694/how-to-save-a-data-frame-in-a-txt-or-excel-file-separated-by-columns

write.table(dfData_For_NCBI,"Data_For_NCBI.txt",sep="\t",row.names=FALSE) 
write.table(dfData_For_BLAST,"Data_For_BLAST.txt",sep="\t",row.names=FALSE) 
