# Reverse_Translate Function
# Anchitaa Ghag

# This function requires two inputs. The name of a file containing the amino acid sequences and a logical TRUE/FALSE input for the type of codon "usage" to use.
# The codon usage (i.e. proportions) will be calculated by the function below. However, the user will need to provide a list of empirical codon counts or frequencies per codon.

Reverse_Translate <- function(sequence_file, default_codon_usage) {
  
  # Read in the protein sequence file.

  SeqFile <- readLines(sequence_file)
  
  # Recursively keep only odd numbered lines i.e., the protein names or headers ">"
  
  Names <- SeqFile[c(TRUE,FALSE)] 
  
  # Recursively keep only even numbered lines i.e., the amino acid sequences
  
  All_Sequences <- SeqFile[c(FALSE, TRUE)] 
  
  rm(SeqFile)
  
  # Read in the default codon usage table or provide an option for the user to enter the name of a file with a custom codon usage table.
  
  if (default_codon_usage == TRUE) { 
    
    Codon_Usage_Table <- read.table("Default_Codon_Counts.txt",
                                    header = TRUE)
    
  } else if  (default_codon_usage == FALSE) {
    
    Custom_Codon_Usage_File <- readline(prompt="Please enter the name of the file with your custom counts per codon information: ") # Adapted from: https://stackoverflow.com/questions/60090558/r-language-user-input-if-condition
    Codon_Usage_Table <- read.table(Custom_Codon_Usage_File,
                                    header = TRUE)
    
  } else {
    
    print("A codon usage was not indicated. The default codon usage of Nicotiana benthamiana will be used.")
    Codon_Usage_Table <-  read.table("Default_Codon_Counts.txt",
                                     header = TRUE)
  }
  
  # Create an Amino_Acid_Lookup list that contains the codons and proportions per amino acid. This will be done using the codon usage from above.
  # Ensure that the Codon_Usage_Table is listed alphabetically be codon.
  
  dfAmAcid_Prob <- Codon_Usage_Table[order(Codon_Usage_Table$Codons),]
  rm(Codon_Usage_Table)
  
  # Assign the single letter abbreviation for all amino acids. This is assuming that codons are listed in an alphabetical order ranging from AAA to TTT (as done in the previous step). In addition, "X"s are assigned for the three stop codons
  
  dfAmAcid_Prob["Single_Letter_Abbreviation"] <- c("K","N","K","N","T","T","T","T","R","S","R","S","I","I","M","I","Q","H","Q","H","P","P","P","P","R","R","R","R","L","L","L","L","E","D","E","D","A","A","A","A","G","G","G","G","V","V","V","V","X","Y","X","Y","S","S","S","S","X","C","W","C","L","F","L","F")
  
  Amino_Acids <- list("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","X","Y")
  
  # Create an empty list to hold the proportions of codons per amino acid. (i.e., the codon usage)
  
  Amino_Acid_Lookup <- list()
  
  for (i in Amino_Acids) {
    
    Codon_Names <- dfAmAcid_Prob %>% 
      filter(Single_Letter_Abbreviation == i) %>%
      select(Codons)
    
    Total <- dfAmAcid_Prob %>% 
      filter(Single_Letter_Abbreviation == i) %>%
      select(Counts) %>%
      colSums() 
    
    Prop <- dfAmAcid_Prob %>% 
      filter(Single_Letter_Abbreviation == i) %>%
      select(Counts) %>%
      summarise(Proportions = (Counts/Total)*1) 
    
    Info <- c(Prop,Codon_Names)
    
    Amino_Acid_Lookup <- append(Amino_Acid_Lookup, list(Info))
    
  }
  
  # Rename the lists to match the single letter abbreviations of the amino acids.
  
  names(Amino_Acid_Lookup) <- unlist(Amino_Acids)
  
  rm(dfAmAcid_Prob,Codon_Names, Total, Prop, Info, Amino_Acids)
  
  # Create an empty list to store all the correspoinding DNA sequences.
  
  All_DNA_Seqs <- list()
  
  # Next, reiterate through each sequence in the All_Sequences list and randomly replace each amino acid with a corresponding triplet/codon based on the proportions in the Amino_Acid_Lookup.
  
  for (One_Sequence in All_Sequences) {
    
    Protein_Sequence <- s2c(toupper(One_Sequence))
    
    DNA_Sequence <- list()
    
    for (Amino_Acid in Protein_Sequence) {
      
      if (Amino_Acid == "A") { 
        
        Triplet <- sample(x = Amino_Acid_Lookup$A$Codons,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$A$Proportions)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "C") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$C$Codons,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$C$Proportions)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "D") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$D$Codons,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$D$Proportions)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "E") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$E$Codons,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$E$Proportions)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "F") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$F$Codons,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$F$Proportions)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "G") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$G$Codons,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$G$Proportions)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "H") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$H$Codons,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$H$Proportions)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "I") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$I$Codons,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$I$Proportions)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "K") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$K$Codons,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$K$Proportions)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "L") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$L$Codons,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$L$Proportions)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "M") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$M$Codons,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$M$Proportions)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "N") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$N$Codons,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$N$Proportions)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "P") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$P$Codons,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$P$Proportions)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "Q") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$Q$Codons,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$Q$Proportions)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "R") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$R$Codons,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$R$Proportions)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "S") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$S$Codons,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$S$Proportions)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "T") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$T$Codons,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$T$Proportions)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "V") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$V$Codons,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$V$Proportions)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "W") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$W$Codons,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$W$Proportions)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "Y") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$Y$Codons,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$Y$Proportions)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else {
        
        DNA_Sequence <- append(DNA_Sequence, "NNN") 
        
        # If an unknown amino acid single letter code is encountered, replace with "NNN". This could be any of the 64 triplet.
      }
    }
    
    # Next, add a stop codon to the end of the DNA sequence using the Amino_Acid_Lookup list for "X" or stop codons.
    
    Triplet <- sample(x = Amino_Acid_Lookup$X$Codons,
                      size = 1,
                      replace = FALSE,
                      prob = Amino_Acid_Lookup$X$Proportions)
    
    DNA_Sequence <- append(DNA_Sequence, Triplet)
    
    # Then, collapse the DNA_Sequence "list of lists" into one long string.
    
    Formatted_DNA_Seq <- (paste((unlist(DNA_Sequence)), collapse = ""))
    
    # Finally, append this sequence to the list of All_DNA_Seqs.
    
    All_DNA_Seqs <- append(All_DNA_Seqs,Formatted_DNA_Seq)
    
  }
  
  # Output the results of alternating protein names in the Names list and the DNA sequences in the All_DNA_Seqs list.
  
  Results <- unlist(c(rbind(Names, All_DNA_Seqs)))
  
  return(Results)
  
}

