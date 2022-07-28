# Functions for Statistical Analysis

#### Plotting Function ####

# In this section, we will calculate the average (i.e. mean) codon counts per codon and the standard deviation.
# This step will be performed for each of the 64 amino acids across the three samples (i.e. data frames).

# Code for std dev calculation across columns adapted from : https://www.datasciencemadesimple.com/get-standard-deviation-of-a-column-in-r-2/

Calculate_Means_Std_Dev <- function(matrix) {
  
  # First, for every column (i.e. triplet), calculate the average.
  Averages <- as.data.frame(colMeans(matrix))
  # Then, for every column (i.e. triplet), calculate the standard deviation.
  StdDevs <- as.data.frame(matrix) %>% 
    summarise_if(is.numeric, sd) %>%
    t()
  # Total number of counts per column.
  Sums <- as.data.frame(colSums(matrix))  
  # (Sum of counts of codons / Total number of codons)*1000
  Frequencies <- as.data.frame((Sums/sum(matrix))*1000) 
  # Add all the lists into a data frame.
  Dataframe <- cbind(Averages, StdDevs, Sums, Frequencies)
  # Rename the columns.
  colnames(Dataframe) <- c("Average_Codon_Count","Std_Deviation","Sum of Counts Per Codon","Frequency (Per 1000 Codons)")
  
  return(Dataframe)
}


#### Plot_AmAcid() Function ####

# Takes an argument "amino_acid" which is the name of the amino acid. A data frame containing average and stand daeviation information.

Plot_AmAcid <- function(amino_acid, df) {
  
  if (amino_acid == "Alanine") { 
    Three_Letter_Form <- "Ala"
  } else if (amino_acid == "Arginine") {
    Three_Letter_Form <- "Arg"
  } else if  (amino_acid == "Asparagine") {
    Three_Letter_Form <- "Asn"
  } else if  (amino_acid == "Aspartic acid") {
    Three_Letter_Form <- "Asp"
  } else if  (amino_acid == "Cysteine") {
    Three_Letter_Form <- "Cys"
  } else if  (amino_acid == "Glutamine") {
    Three_Letter_Form <- "Gln"
  } else if  (amino_acid == "Glutamic acid") {
    Three_Letter_Form <- "Glu"
  } else if  (amino_acid == "Glycine") {
    Three_Letter_Form <- "Gly"
  } else if  (amino_acid == "Histidine") {
    Three_Letter_Form <- "His"
  } else if  (amino_acid == "Isoleucine") {
    Three_Letter_Form <- "Ile"
  } else if  (amino_acid == "Leucine") {
    Three_Letter_Form <- "Leu"
  } else if  (amino_acid == "Lysine") {
    Three_Letter_Form <- "Lys"
  } else if  (amino_acid == "Methionine") {
    Three_Letter_Form <- "Met"
  } else if  (amino_acid == "Phenylalanine") {
    Three_Letter_Form <- "Phe"
  } else if  (amino_acid == "Proline") {
    Three_Letter_Form <- "Pro"
  } else if  (amino_acid == "Serine") {
    Three_Letter_Form <- "Ser"
  } else if  (amino_acid == "Threonine") {
    Three_Letter_Form <- "Thr"
  } else if  (amino_acid == "Tryptophan") {
    Three_Letter_Form <- "Trp"
  } else if  (amino_acid == "Tyrosine") {
    Three_Letter_Form <- "Tyr"
  } else {
    Three_Letter_Form <- "Val"
  }
  
  # Creates a data frame containing the codon count averages and std dev from original data frame.
  
  dfAmAcid <- (subset(df, AmAcid == Three_Letter_Form))[,c(1,3,4,5)]
  
  # Generate the plot using ggplot package functions.
  
  Plot <- ggplot(data=dfAmAcid, 
                 aes(x=Codon, y=Avrgs, fill=Sample)) +
    geom_bar(stat="identity", position=position_dodge(), alpha = 0.8) +
    ggtitle(amino_acid) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5, size=10)) + #, legend.position = "none") +
    xlab("") +
    ylab("") +
    geom_errorbar(aes(ymin=Avrgs-SD, ymax=Avrgs+SD), width=.2, position=position_dodge(.9)) +
    scale_y_continuous(breaks = seq(0, 55, by = 10), limits = c(0,55)) + # adapted from: https://stackoverflow.com/questions/37950511/r-ggplot2-setting-tick-mark-interval
    scale_fill_manual(values=c("#FDE725", "#21918C", "#440154")) +
    scale_color_manual(name= "Species",
                       labels = c("Old N.benthamiana (Kazuza)", "New N.benthamiana ", "N.tabacum (Kazuza)"))
  
  return(Plot)
  
}

#### REFERENCES ####
