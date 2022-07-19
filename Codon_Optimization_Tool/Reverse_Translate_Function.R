# Reverse_Translate Function
# Anchitaa Ghag

# The following two packages and associated libraries need to be installed prior to running this function.

#install.packages("dplyr") 
#install.packages("seqinr") 

library("dplyr") 
library("seqinr")

# This function requires two input files. The name of a file containing the amino acid sequences and a file containing the codon usage to use.

Reverse_Translate <- function(sequence_file, codon_usage_table) {
  
  # Read in the protein sequence file.
  
  SeqFile <- readLines(sequence_file)
  
  # Recursively keep only odd numbered lines i.e., the protein names or headers ">"
  
  Names <- SeqFile[c(TRUE,FALSE)] 
  
  # Recursively keep only even numbered lines i.e., the amino acid sequences
  
  All_Sequences <- SeqFile[c(FALSE, TRUE)] 
  
  rm(SeqFile)
  
  # Read in the codon usage table file.
  
  CU_Table <- read.table(codon_usage_table,
                         header = TRUE)
  
  # Create an Amino_Acid_Lookup list that contains the codons and proportions per amino acid. This will be done using the codon usage from above.
  # Ensure that the Codon_Usage_Table is listed alphabetically be codon.
  
  dfAmAcid_Prob <- CU_Table[order(CU_Table$Codon),]
  
  rm(CU_Table)
  
  # Create a list of all amino acids listed as their single letter abbreviation.
  
  Amino_Acids <- list("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","X","Y")
  
  # Create an empty list to hold the proportions of codons per amino acid. (i.e., the codon usage)
  
  Amino_Acid_Lookup <- list()
  
  # Generate the Amino_Acid_Lookup data
  
  for (Each_Amino in Amino_Acids) {
    
    # First, get the codon "name" from the dataframe.
    
    Codon_Names <- dfAmAcid_Prob %>% 
      filter(Single_Letter_Abbreviation == Each_Amino) %>%
      select(Codon)
    
    # Then, calculate the total sum of all codon counts.
    
    Total <- dfAmAcid_Prob %>% 
      filter(Single_Letter_Abbreviation == Each_Amino) %>%
      select(Number) %>%
      colSums() 
    
    # Finally, calculate the probability for the codon to be encountered.
    
    Prop <- dfAmAcid_Prob %>% 
      filter(Single_Letter_Abbreviation == Each_Amino) %>%
      select(Number) %>%
      summarise(Proportion = (Number/Total)*1) 
    
    # Add this information to the Amino_Acid_Lookup list.
    
    Amino_Acid_Lookup <- append(Amino_Acid_Lookup, list(c(Prop,Codon_Names)))
    
  }
  
  # Rename the lists to match the single letter abbreviations of the amino acids.
  
  names(Amino_Acid_Lookup) <- unlist(Amino_Acids)
  
  rm(dfAmAcid_Prob, Codon_Names, Total, Prop, Info, Amino_Acids, Each_Amino)
  
  # Create an empty list to store all the corresponding DNA sequences.
  
  All_DNA_Seqs <- list()
  
  # Next, reiterate through each sequence in the All_Sequences list and randomly replace each amino acid with a corresponding triplet/codon based on the proportions in the Amino_Acid_Lookup.
  
  for (One_Sequence in All_Sequences) {
    
    Protein_Sequence <- s2c(toupper(One_Sequence))
    
    DNA_Sequence <- list()
    
    for (Amino_Acid in Protein_Sequence) {
      
      if (Amino_Acid == "A") { 
        
        Triplet <- sample(x = Amino_Acid_Lookup$A$Codon,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$A$Proportion)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "C") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$C$Codon,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$C$Proportion)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "D") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$D$Codon,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$D$Proportion)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "E") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$E$Codon,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$E$Proportion)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "F") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$F$Codon,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$F$Proportion)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "G") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$G$Codon,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$G$Proportion)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "H") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$H$Codon,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$H$Proportion)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "I") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$I$Codon,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$I$Proportion)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "K") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$K$Codon,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$K$Proportion)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "L") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$L$Codon,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$L$Proportion)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "M") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$M$Codon,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$M$Proportion)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "N") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$N$Codon,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$N$Proportion)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "P") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$P$Codon,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$P$Proportion)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "Q") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$Q$Codon,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$Q$Proportion)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "R") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$R$Codon,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$R$Proportion)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "S") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$S$Codon,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$S$Proportion)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "T") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$T$Codon,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$T$Proportion)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "V") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$V$Codon,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$V$Proportion)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "W") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$W$Codon,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$W$Proportion)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else if  (Amino_Acid == "Y") {
        
        Triplet <- sample(x = Amino_Acid_Lookup$Y$Codon,
                          size = 1,
                          replace = FALSE,
                          prob = Amino_Acid_Lookup$Y$Proportion)
        
        DNA_Sequence <- append(DNA_Sequence, Triplet)
        
      } else {
        
        DNA_Sequence <- append(DNA_Sequence, "NNN") 
        
        # If an unknown amino acid single letter code is encountered, replace with "NNN". This could be any of the 64 triplet.
      }
    }
    
    # Next, add a stop codon to the end of the DNA sequence using the Amino_Acid_Lookup list for "X" or stop codons.
    
    Triplet <- sample(x = Amino_Acid_Lookup$X$Codon,
                      size = 1,
                      replace = FALSE,
                      prob = Amino_Acid_Lookup$X$Proportion)
    
    DNA_Sequence <- append(DNA_Sequence, Triplet)
    
    # Then, collapse the DNA_Sequence "list of lists" into one long string.
    
    Formatted_DNA_Seq <- (paste((unlist(DNA_Sequence)), collapse = ""))
    
    # Finally, append this sequence to the list of All_DNA_Seqs.
    
    All_DNA_Seqs <- append(All_DNA_Seqs,Formatted_DNA_Seq)
    
  }
  
  rm(Amino_Acid_Lookup,DNA_Sequence,All_Sequences,Amino_Acid,Formatted_DNA_Seq,One_Sequence,Protein_Sequence,Triplet)
  
  # Output the results of alternating protein names in the Names list and the DNA sequences in the All_DNA_Seqs list.
  
  Results <- unlist(c(rbind(Names, All_DNA_Seqs)))
  
  rm(All_DNA_Seqs, Names)
  
  return(Results)
  
}

