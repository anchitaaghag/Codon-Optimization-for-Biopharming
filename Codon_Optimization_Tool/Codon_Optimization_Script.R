# Codon Optimization Script
# Anchitaa Ghag

# This R code will require three user inputs in section 01 and will then run the Reverse_Translate() function to optimize codons in a chosen organism.
# If running this script for the first time, ensure that the two necessary packages have been installed by uncommenting (i.e. removing the "#") and running lines 33 & 34.

#### 01 USER INPUT ####

# First, ensure you are in the same working directory as your input files by running the following lines (if required).
# getwd()
# setwd()

# Next, please indicate the name of the file with your codon usage table. 
# If you would like to use the updated codon usage of N. benthamiana, leave this as "Updated_Codon_Usage.txt".

Chosen_Codon_Usage <- "Updated_Codon_Usage.txt"

# Then, please indicate the name of the file with your protein sequence(s). 
# If you would like to use the example protein sequences, leave this as "Example_Protein_Sequences.txt".

Protein_Sequence_File <- "Example_Protein_Sequences.txt"

# Finally, please indicate a name for your output results file.
# If you would like to use the example, leave this as "Example_Results.txt".

Result_File_Name <- "Example_Results.txt"

#### 02 LOAD FUNCTIONS ####

# This script requires the following packages and libraries.
# If these packages have not yet been installed please remove the "#" and run the following lines.

#install.packages("dplyr") 
#install.packages("seqinr") 

library("dplyr") 
library("seqinr")

# Load the reverse translate function from the R file.

source("Reverse_Translate_Function.R")

#### 03 REVERSE TRANSLATE AMINO ACID SEQUENCES ####

# Reverse translate the amino acid sequences to DNA sequences based on the user input defined above.

DNA_Sequences <- Reverse_Translate(sequence_file = Protein_Sequence_File, 
                                   codon_usage_table = Chosen_Codon_Usage)

#### 04 VIEW THE RESULTS ####

# View the first sequence that has been translated.

head(DNA_Sequences, 2)

#### 05 OUTPUT RESULTS TO FILE ####

writeLines(text = DNA_Sequences, 
           con = Result_File_Name)
