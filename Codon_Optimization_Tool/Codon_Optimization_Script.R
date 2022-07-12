# Codon Optimization Script
# Anchitaa Ghag

# This R code will require three user inputs in section 01 and will then run the Reverse_Translate() function to optimize codons in a chosen organism. 

#### 01 USER INPUT ####

# First, ensure you are in the same working directory as your input files by running the following lines (if required).
# getwd()
# setwd()

# Next, please indicate if you would like to used the updated codon usage of N. benthamiana ("TRUE") or if you would like to provide your own codon usage ("FALSE"). You will be prompted to enter the name of the file with your custom codon usage in the console if the "FALSE" option is chosen.
# If you would like to use the example, leave this as "TRUE".

Chosen_Option <- TRUE

# Then, please indicate the name of the file with your protein sequence(s). 
# If you would like to use the example protein sequences, leave this as "Example_Protein_Sequences.txt"

Protein_Sequence_File <- "Example_Protein_Sequences.txt"

# Finally, please indicate a name for your output results file.
# If you would like to use the example, leave this as "Example_Results.txt"

Result_File_Name <- "Example_Results.txt"

#### 02 LOAD FUNCTION ####

# Load the reverse translate function.

source("Reverse_Translate_Function.R")

#### 03 REVERSE TRANSLATE AMINO ACID SEQUENCES ####

# Reverse translate the amino acid sequences to DNA sequences based on the user input defined above.

DNA_Sequences <- Reverse_Translate(sequence_file = Protein_Sequence_File, 
                                   default_codon_usage = Chosen_Option)

#### 04 VIEW THE RESULTS ####

# View the first sequence that has been translated.

head(DNA_Sequences, 1)

#### 05 OUTPUT RESULTS TO FILE ####

writeLines(text = DNA_Sequences, 
           con = Result_File_Name)
