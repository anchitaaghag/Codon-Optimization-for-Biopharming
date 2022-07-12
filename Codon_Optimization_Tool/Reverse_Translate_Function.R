# Reverse_Translate Function
# Anchitaa Ghag

Reverse_Translate <- function(sequence_file, output_file_name, default_codon_usage) {
  
  # Read in the default codon usage table or a custom codon usage table.
  
  if (default_codon_usage == TRUE) { 
    Codon_Usage_Table <- read.table("Default_CU.txt")
  } else if  (default_codon_usage == FALSE) {
    Custom_Codon_Usage_File <- readline(prompt="Please enter the name of the file with your custom codon usage table: ") # Adapted from: https://stackoverflow.com/questions/60090558/r-language-user-input-if-condition
    Codon_Usage_Table <- read.table(Custom_Codon_Usage_File)
  } else {
    print("A codon usage table was not indicated. The default codon usage of Nicotiana benthamiana will be used.")
    Codon_Usage_Table <-  read.table("Default_CU.txt")
  }
  
  # Create an Amino_Acid_Lookup after reading in the codon usage.
  
  # Note this is based on an alphabetical order ranging from AAA to TTT.
  
  Single_Letter_Abbreviation <- c("K","N","K","N","T","T","T","T","R","S","R","S","I","I","M","I","Q","H","Q","H","P","P","P","P","R","R","R","R","L","L","L","L","E","D","E","D","A","A","A","A","G","G","G","G","V","V","V","V","X","Y","X","Y","S","S","S","S","X","C","W","C","L","F","L","F")
  
  
  dfAmAcid_Prob <- data.frame(dfNew_MeanSD$`Sum of Counts Per Codon`,rownames(dfNew_MeanSD), Single_Letter_Abbreviation)
  colnames(dfAmAcid_Prob) <- c("Counts_Per_Codon", "Codons", "Single_Letter_Abbreviation")
  
  
  # "X" for the three stop codons
  
  
  
  Amino_Acids <- list("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","X","Y")
  
  
  Amino_Acid_Chart <- list()
  
  for (i in Amino_Acids){
    
    Codon_Names <- dfAmAcid_Prob %>% filter(Single_Letter_Abbreviation == i) %>%
      select(Codons)
    
    Total <- dfAmAcid_Prob %>% filter(Single_Letter_Abbreviation == i) %>%
      select(Counts_Per_Codon) %>%
      colSums() 
    
    Prop <- dfAmAcid_Prob %>% filter(Single_Letter_Abbreviation == i) %>%
      select(Counts_Per_Codon) %>%
      summarise(Proportions = (Counts_Per_Codon/Total)*1) 
    
    Info <- c(Prop,Codon_Names)
    
    Amino_Acid_Chart <- append(Amino_Acid_Chart, list(Info))
    
    
    
  }
  
  names(Amino_Acid_Chart) <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","X","Y")
  
  
  
  # Check the sequence file format to ensure the file is in a .txt or .fasta format.
  
  File_Format <- sub(".*\\.", "", Protein_Sequence_File) # https://stackoverflow.com/questions/31774086/extracting-text-after-last-period-in-string
  
  if (File_Format == ".txt") { 
    #  dfNames_Sequences <-
    #} else if  (File_Format == ".fasta") {
    dfNames_Sequences <-
      #} else {
      print("The protein sequence file format is not currently supported. Please provide a text or fasta file.")
  }
  
  
  All_Sequences <- as.list(dfNames_Sequences$Sequences)
  
  # Random replacement of values.
  
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
    
    # Next, add a stop codon.
    
    Triplet <- sample(x = Amino_Acid_Lookup$X$Codons,
                      size = 1,
                      replace = FALSE,
                      prob = Amino_Acid_Lookup$X$Proportions)
    
    DNA_Sequence <- append(DNA_Sequence, Triplet)
    
    
    
    print(paste((unlist(DNA_Sequence)), collapse = ""))
    
    print(c("Your results have been written to:", Output_File_Name))
    return()
    
  }
}
