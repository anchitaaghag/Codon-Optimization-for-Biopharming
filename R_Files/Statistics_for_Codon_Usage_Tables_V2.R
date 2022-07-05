# Statistical Anaylsis of Update Codon Usage in Nicotiana benthaminana
# 27 June 2022
# Anchitaa Ghag

#### 00 PERSONAL NOTES ####

setwd("/Users/anchitaa/Major_Research_Project_2022/06_Code/08_Statistical_Analysis/")

# Still to do:

# Make the axis consistent 0 -> 55  Done. Need to ensure it displays well for AmAcid plots.
# change colors to a color-blind friendly palette like viridis Done. Added transparancy as well to ensure error bars can still be seen.
# Clean up this script, ensure commented throughout
# Email Mr. Muselius & Dr. JGM about the updated Intra B Plot
# Also update Dr. AHW!
# Add another figure perhaps a colored table with the RSCU and ENC values? How to visualize that?

#### 01 INSTALL PACKAGES & DOWNLOAD LIBRARIES ####

# First, set the working directory by running the following lines.
# setwd()
# getwd()

# This script requires the following packages and libraries.
# If these packages have not yet been installed please remove the "#" and run the following lines.

#BiocManager::install("Biostrings")
#install.packages("coRdon")
#install.packages("ggplot2")
#install.packages("gridExtra")
#install.packages("stringr")
#install.packages("tidyverse")
#install.packages("viridis")

library("Biostrings")
library("coRdon")
library("ggplot2")
library("gridExtra")
library("stringr")
library("tidyverse")
library("viridis")

#### 02 DATA AQUISITION : IMPORT CU FROM PREVIOUS SCRIPT ####

dfNew <- read_csv("dfCodingSeqs.csv")[,2:14]
dfOld <- read_csv("dfKazuza.csv")[,2:5]
dfTabacum <- read_csv("dfNTabacum.csv")[,2:5]

#### 03 MINOR NOTE ON AVERAGE CODON COUNTS REPORTED BY KAZUZA ####

# Note that the average number of codon counts per 1000 codons reported by the Kazuza database can differ slightly.
# For example, check the total number of codons (should be 1000) in Kazuza CU for N. benthamiana and N. tabacum.

sum(dfOld$X.1000) # Kazuza CU for N. benthamiana reports 999.99 instead of 1000

sum(dfTabacum$X.1000) # Kazuza CU for N. tabacum reports 1000.04 instead of 1000

# While these are minor departures from the "1000 codon" value, these average codon counts will be recalculated to ensure that the frequencies reported are as exactly equal as possible.
# To do this, need to import the complete set of codon counts for each coding sequence used to generate Kazuza's CU tables.
# These can be found on the individual record in the database under the "List of codon usage for each CDS (format)" link.

# For N. benthamiana: http://www.kazusa.or.jp/codon/current/species/4100
# For N. tabacum: http://www.kazusa.or.jp/codon/current/species/4097

# For easier loading into R, I will be loading previously formatted files (using command line).

Old_CDS <- read.table("N_benthamiana_Codon_Counts_Only.txt",
                        header = TRUE)

Tabacum_CDS <- read.table("N_tabacum_Codon_Counts_Only.txt",
                            header = TRUE)

#### 04 GENERATE CODON USAGE TABLES ####

#https://bioconductor.riken.jp/packages/3.8/bioc/vignettes/coRdon/inst/doc/coRdon.html

# 1) First, convert the trimmed coding sequences to a DNAStringSet object.

cds <- DNAStringSet(dfNew$Trimmed_CDS)

# 2) Second, create a codon usage table object for both existing and current codon usage tables.
# NOTE! The columns in the list of codon counts imported from Kazuza are not sorted alphabetically! 
# These need to be reordered before inputting into the codonTable function or else counts will be assined to the wrong codon!

# Replace U w T

colnames(Old_CDS) <- str_replace_all(colnames(Old_CDS),"U","T") 
colnames(Tabacum_CDS) <- str_replace_all(colnames(Tabacum_CDS),"U","T") 

# Order the columns alphabetically.
Old_CDS_Ordered <- Old_CDS %>%
  select(order(colnames(Old_CDS)))

Tabacum_CDS_Ordered <- Tabacum_CDS %>%
  select(order(colnames(Tabacum_CDS)))

# 

Old_CU <- codonTable(Old_CDS_Ordered)
New_CU <- codonTable(cds)
Tabacum_CU <- codonTable(Tabacum_CDS_Ordered)

# 3) Then, using the codonCounts() function from the codRon package create a matrix to view the counts per codon.

Old_Matrix <- codonCounts(Old_CU)
New_Matrix <- codonCounts(New_CU)
Tabacum_Matrix <- codonCounts(Tabacum_CU)

#### 05 MEAN & STANDARD DEVIATIONS ####

# In this section, we will calculate the average (i.e. mean) codon counts per codon and the standard deviation.
# This step will be performed for each of the 64 amino acids across the three samples (i.e. data frames).

# Code for std dev calculation across columns adapted from : https://www.datasciencemadesimple.com/get-standard-deviation-of-a-column-in-r-2/

# For Kazuza's N. benthamiana coding sequences:

# First, for every column (i.e. triplet), calculate the average.
Averages <- as.data.frame(colMeans(Old_Matrix))
# Then, for every column (i.e. triplet), calculate the standard deviation.
StdDevs <- as.data.frame(Old_Matrix) %>% 
  summarise_if(is.numeric, sd) %>%
  t()
# Finally, transform both the average and the standard deviation.
# Linearly transforming the mean will ensure that all the average codon counts are reported on the same scale. That is, all the data is on the same scale.
# This will also include linearly transforming the standard deviation, which will not be affected by the transformation.
Sums <- as.data.frame(colSums(Old_Matrix))  # Total number of counts per column.
Frequencies <- as.data.frame((Sums/sum(Old_Matrix))*1000) # (Sum of counts of codons / Total number of codons)*1000

dfOld_MeanSD <- cbind(Averages, StdDevs, Sums, Frequencies)
colnames(dfOld_MeanSD) <- c("Average_Codon_Count","Std_Deviation","Sum of Counts Per Codon","Frequency (Per 1000 Codons)")

#dfNew_MeanSD

# For the new  N. benthamiana coding sequences:

# First, for every column (i.e. triplet), calculate the average.
Averages <- as.data.frame(colMeans(New_Matrix))
# Then, for every column (i.e. triplet), calculate the standard deviation.
StdDevs <- as.data.frame(New_Matrix) %>% 
  summarise_if(is.numeric, sd) %>%
  t()
# Finally, transform both the average and the standard deviation.
# Linearly transforming the mean will ensure that all the average codon counts are reported on the same scale. That is, all the data is on the same scale.
# This will also include linearly transforming the standard deviation, which will not be affected by the transformation.
Sums <- as.data.frame(colSums(New_Matrix))  # Total number of counts per column.
Frequencies <- as.data.frame((Sums/sum(New_Matrix))*1000) # (Sum of counts of codons / Total number of codons)*1000

dfNew_MeanSD <- cbind(Averages, StdDevs, Sums, Frequencies)
colnames(dfNew_MeanSD) <- c("Average_Codon_Count","Std_Deviation","Sum of Counts Per Codon","Frequency (Per 1000 Codons)")


#dfTabacum_MeanSD

# For the Kazuza's  N. tabcum coding sequences:

# First, for every column (i.e. triplet), calculate the average.
Averages <- as.data.frame(colMeans(Tabacum_Matrix))
# Then, for every column (i.e. triplet), calculate the standard deviation.
StdDevs <- as.data.frame(Tabacum_Matrix) %>% 
  summarise_if(is.numeric, sd) %>%
  t()
# Finally, transform both the average and the standard deviation.
# Linearly transforming the mean will ensure that all the average codon counts are reported on the same scale. That is, all the data is on the same scale.
# This will also include linearly transforming the standard deviation, which will not be affected by the transformation.
Sums <- as.data.frame(colSums(Tabacum_Matrix))  # Total number of counts per column.
Frequencies <- as.data.frame((Sums/sum(Tabacum_Matrix))*1000) # (Sum of counts of codons / Total number of codons)*1000

dfTabacum_MeanSD <- cbind(Averages, StdDevs, Sums, Frequencies)
colnames(dfTabacum_MeanSD) <- c("Average_Codon_Count","Std_Deviation","Sum of Counts Per Codon","Frequency (Per 1000 Codons)")

#### 06 Histogram #####

#dfMeanSD <-

AmAcid <- rep(rownames(dfOld_MeanSD),3)
Species <- c((rep("Old N. benthamiana (Kazuza)",64)), (rep("New N. benthamiana",64)), (rep("N. tabacum (Kazuza)",64)))
Freq <- c(dfOld_MeanSD$`Frequency (Per 1000 Codons)`, dfNew_MeanSD$`Frequency (Per 1000 Codons)`, dfTabacum_MeanSD$`Frequency (Per 1000 Codons)`)
dfMeanSD <- data.frame(AmAcid,Species,Freq)

# Compare the codon usage frequency using a histogram.

ggplot(data=dfMeanSD, 
      aes(x=AmAcid, y=Freq, fill=Species)) +
  geom_bar(stat="identity", position=position_dodge(), show.legend = TRUE) +
  ggtitle("Histogram of Codon Usage Frequency") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Codons") +
  ylab("Frequency (Per 1000 Codons)") +
  scale_x_discrete(guide = guide_axis(angle = 90)) + #https://stackoverflow.com/questions/1330989/rotating-and-spacing-axis-labels-in-ggplot2
  theme(legend.position = c(0.93, 0.93),
        legend.background = element_rect(fill = "white", color = "black")) + #https://datavizpyr.com/how-to-place-legend-inside-the-plot-with-ggplot2/ 
  ylim(NA,40) +
  scale_fill_viridis(discrete = TRUE, option = "viridis") 

#### 07 CHI SQUARE TEST ####

# Perform a chi-squared test to evaluate if there is a difference between the frequencies of both CU.

chisq.test(x = dfOld_MeanSD$`Frequency (Per 1000 Codons)`,
           y = dfNew_MeanSD$`Frequency (Per 1000 Codons)`)

# Pearson's Chi-squared test

# data:  dfOld_MeanSD$`Frequency (Per 1000 Codons)` and dfNew_MeanSD$`Frequency (Per 1000 Codons)`
# X-squared = 3904, df = 3843, p-value = 0.242

# https://www.bioconductor.org/packages/devel/bioc/vignettes/coRdon/inst/doc/coRdon.html#calculate-cu-bias

# "Other CU statistics can be calculated in the same way as MILC(), using one of the functions: B(), MCB(), ENCprime(), ENC() or SCUO(). Note however, that when calculating ENC and SCUO, one doesn’t need to provide a subset of refe"
# Next, compare the CU bias for every coding sequence between the created CUT and Kazuza through visualizations.

#### 08 Visualization of Codon Usage - Existing vs. Current Karlin B plot ####

# Code Adapted From: https://www.bioconductor.org/packages/devel/bioc/vignettes/coRdon/inst/doc/coRdon.html#calculate-cu-bias

# Use the intraBplot() function in the coRdon package and our two codonTable objects to plot a  plot of codon usage distances between existing vs. created codon tables.
# Intra-samples Karlin B plot

set.seed(1111)

# Randomly sample the same number of "points" (i.e. coding sequences) as the old/exisitng codon usage.

Sample_CU <- sample(New_CU,length(Old_CU)) 

intraBplot(x = Old_CU, 
           y = Sample_CU, 
          names = c("Old", "New"), 
           variable = "MILC", 
           size = 3, 
           alpha = 0.8) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  ggtitle("Karlin B Plot of Existing (Kazuza) v.s.New CU Distances")+
  xlab("MILC Distance from Existing CU") +
  ylab("MILC Distance from Consensus CU") + 
 # scale_fill_manual(name = "Species", labels = c("Old N.benthamiana (Kazuza)", "New N.benthamiana")) +
  theme(legend.position = c(0.93, 0.10),
        legend.background = element_rect(fill = "white", color = "black")) +  #https://datavizpyr.com/how-to-place-legend-inside-the-plot-with-ggplot2/ 
  labs(color="Codon Usage") # https://stackoverflow.com/questions/14622421/how-to-change-legend-title-in-ggplot
  #geom_jitter(width  = 0.05) + # Add noise since the dataset in small to avoid overplotting.
  #geom_smooth(method = "lm", formula = y ~ x, fullrange = FALSE, level = 0.95) # 95% confidence interval
  #geom_smooth(aes(group = 1), formula = y ~ poly(x,2), method = "lm", fullrange = TRUE, level = 0.95) 
  #geom_curve()

# The example dataset in the vignette has ~ 19, 000 + "points". There may just be undersampling/ less data points
# This makes sense, since it is the same species, the codon usage should not differ greatly from the Old CU data set available.



# Compare with N.tabacum 

set.seed(1234)
Sample_CU <- sample(Tabacum_CU,407)

intraBplot(x = New_CU, 
           y = Sample_CU, 
           names = c("Existing", "Consensus"), 
           variable = "MILC", 
           size = 1, 
           alpha = 1.0) + 
  ggtitle("Karlin B Plot of Existing (Kazuza) v.s. Consensus CU Distances")+
  xlab("MILC Distance from Existing CU") +
  ylab("MILC Distance from Consensus CU") + 
#geom_jitter(width  = 0.05) + # Add noise since the dataset in small to avoid overplotting.
#geom_smooth(method = "lm", formula = y ~ x, fullrange = FALSE, level = 0.95) # 95% confidence interval
#geom_smooth(aes(group = 1), formula = y ~ poly(x,2), method = "lm", fullrange = TRUE, level = 0.95) 
#geom_curve()

# Again, not a lot of difference. As expected.
  
#### 09 ADD SECTION THAT COMPARES THE B() TO SUPPLEMENT INTRABPLOT) ####

#### 10 RSCU ####

formatted_cds <- s2c((unlist(dfNew$Trimmed_CDS[1])))

New_RSCU <- uco(seq=formatted_cds,
    frame = 0,
    index = "rscu",
    as.data.frame = TRUE,
    NA.rscu = NA)


table(New_RSCU$RSCU == 1) # RSCU = 1 codon is used as expected by random usage 
table(New_RSCU$RSCU > 1) # RSCU > 1 codon used more frequently than random
table(New_RSCU$RSCU < 1) # RSCU < 1 codon used less frequently than random 

# Are these different from the RSCU for Old N.benthamiana CU?

# Don't have the sequences for Old N.benthamiana CU. May have to import that (if time permits).
# Can report the RSCU measures in a table (perhaps supplementary results?)

#### 11 ENC & Mann Whitney U Test ####

enc <- ENC(New_CU)
Old_ENC <- ENC(Old_CU)

library(ggpubr)
ggqqplot(enc)
ggqqplot(Old_ENC)

shapiro.test(enc) # W = 0.98266, p-value = 8.508e-05
# p-value < 0.05 
# distribution of data ist significantly different from normal distribution
# not normal

shapiro.test(Old_ENC) #W = 0.97915, p-value = 0.114

# p-value > 0.05 
# distribution of data not significantly different from normal distribution

# A non-parametric alternative to the t-test is Mann Whitney U test
wilcox.test(x = enc, 
       y = Old_ENC,
       alternative = "two.sided")

# Wilcoxon rank sum test with continuity correction

# data:  enc and Old_ENC
# W = 18366, p-value = 0.1309
# alternative hypothesis: true location shift is not equal to 0

# p-value < 0.05 significant
# p-value >= 0.05 not significant
# H0 <- mean 1 - mean 2 = 0
# HA <- mean 1 - mean 2 =/= 0

# The result is not significant at p < .05.
# H0 Null Hyp not rejected& 

#http://www.sthda.com/english/wiki/normality-test-in-r

table(enc == 20) # ENC = 20 codon is used as expected by random usage 
table(enc == 61) # ENC = 61 codon used more frequently than random
table(enc < 40) # ENC < 40 low codon usage bias

#### 12 CREATE AVRAGES & STD DEV DATA FRAME ####

Codon <- rep(rownames(dfOld_MeanSD),3)
AmAcid <- as.list((dfOld[order(dfOld$Codon),])[,1])
Sample <- c((rep("Old",64)), (rep("New",64)), (rep("Tabacum",64)))
Avrgs <- c(dfOld_MeanSD$`Frequency (Per 1000 Codons)`, dfNew_MeanSD$`Frequency (Per 1000 Codons)`, dfTabacum_MeanSD$`Frequency (Per 1000 Codons)`)
SD <- c(dfOld_MeanSD$Std_Deviation, dfNew_MeanSD$Std_Deviation, dfTabacum_MeanSD$Std_Deviation)

dfCodon_Mean_SD <- data.frame(Codon,AmAcid,Sample,Avrgs, SD)

#### 13 PLOT AVERAGE COUNTS PER AMINO ACID #####

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
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  xlab("Codon") +
  ylab("Average Codon Count") +
  geom_errorbar(aes(ymin=Avrgs-SD, ymax=Avrgs+SD), width=.2, position=position_dodge(.9)) +
  ylim(NA,55) +
  scale_fill_viridis(discrete = TRUE, option = "viridis") +
  scale_color_manual(name= "Species",
                   labels = c("Old N.benthamiana (Kazuza)", "New N.benthamiana ", "N.tabacum (Kazuza)"))
                   
 return(Plot)

}
                  

# Alanine - Ala

p1 <- Plot_AmAcid(amino_acid = "Alanine", df = dfCodon_Mean_SD)

# Arginine - Arg 

p2 <- Plot_AmAcid(amino_acid = "Arginine", df = dfCodon_Mean_SD)

# Asparagine - Asn

p3 <- Plot_AmAcid(amino_acid = "Asparagine", 
                  df = dfCodon_Mean_SD)

# Aspartic acid - Asp 

p4 <- Plot_AmAcid(amino_acid = "Aspartic acid", 
                  df = dfCodon_Mean_SD)

# Cysteine - Cys

p5 <- Plot_AmAcid(amino_acid = "Cysteine", 
                  df = dfCodon_Mean_SD)

# Glutamine - Gln

p6 <- Plot_AmAcid(amino_acid = "Glutamine", 
                  df = dfCodon_Mean_SD)

# Glutamic acid - Glu

p7 <- Plot_AmAcid(amino_acid = "Glutamic acid", 
                  df = dfCodon_Mean_SD)

# Glycine - Gly 

p8 <- Plot_AmAcid(amino_acid = "Glycine", 
                  df = dfCodon_Mean_SD)

# Histidine - His 

p9 <- Plot_AmAcid(amino_acid = "Histidine", 
                  df = dfCodon_Mean_SD)

# Isoleucine - Ile 

p10 <- Plot_AmAcid(amino_acid = "Isoleucine", 
                  df = dfCodon_Mean_SD)

# Leucine - Leu 

p11 <- Plot_AmAcid(amino_acid = "Leucine", 
                  df = dfCodon_Mean_SD)

# Lysine - Lys 

p12 <- Plot_AmAcid(amino_acid = "Lysine", 
                  df = dfCodon_Mean_SD)

# Methionine - Met

p13 <- Plot_AmAcid(amino_acid = "Methionine", 
                  df = dfCodon_Mean_SD)

# Phenylalanine - Phe 

p14 <- Plot_AmAcid(amino_acid = "Phenylalanine", 
                  df = dfCodon_Mean_SD)

# Proline - Pro 

p15 <- Plot_AmAcid(amino_acid = "Proline", 
                  df = dfCodon_Mean_SD)

# Serine - Ser

p16 <- Plot_AmAcid(amino_acid = "Serine", 
                  df = dfCodon_Mean_SD)

# Threonine - Thr 

p17 <- Plot_AmAcid(amino_acid = "Threonine", 
                  df = dfCodon_Mean_SD)

# Tryptophan - Trp 

p18 <- Plot_AmAcid(amino_acid = "Tryptophan", 
                  df = dfCodon_Mean_SD)

# Tyrosine - Tyr 

p19 <- Plot_AmAcid(amino_acid = "Tyrosine", 
                  df = dfCodon_Mean_SD)

# Valine - Val 

p20 <- Plot_AmAcid(amino_acid = "Valine", 
                  df = dfCodon_Mean_SD)

#### 14 ARRANGE PLOTS IN A GRID ####

# https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html

grid.arrange(p1, p2, p3, p4, nrow = 2)
grid.arrange( p5, p6, p7, p8, nrow = 2)
grid.arrange( p9, p10, p11, p12, nrow = 2)
grid.arrange( p13, p14, p15, p16, nrow = 2)
grid.arrange( p17, p18, p19, p20, nrow = 2)

#### 15 ADITIONAL WORK ####

# Finally, calculate the most used codons (1 for each amino acid residue) from the data frame.

# dfKazuza.max <- group_by(dfKazuza, AmAcid) %>% 
#  filter(Number==max(Number)) %>% 
#  select(AmAcid, Codon)



#getlen(Current_CU)
#row.names(CU) <- dfCodingSeqs$Titles
#(colMeans(CU)/colSums(CU))*1000


#
#### 16 REFERENCES ####