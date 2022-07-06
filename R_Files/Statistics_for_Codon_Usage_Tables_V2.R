# Statistical Analysis of Update Codon Usage in Nicotiana benthaminana
# 27 June 2022
# Anchitaa Ghag

#### 00 PERSONAL NOTES ####

setwd("/Users/anchitaa/Major_Research_Project_2022/06_Code/08_Statistical_Analysis/")

# Still to do:

# Clean up this script, ensure commented throughout
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
#install.packages("ggpubr") # This package is being used for Q-Q Plots, Is there a simpler base R function that can do this?
#install.packages("gridExtra")
#install.packages("seqinr")
#install.packages("tidyverse")
#install.packages("viridis")

library("Biostrings")
library("coRdon")
library("ggplot2")
library("ggpubr")
library("gridExtra")
library("seqinr")
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

dfOld_MeanSD <- Calculate_Means_Std_Dev(matrix = Old_Matrix)

dfNew_MeanSD <- Calculate_Means_Std_Dev(matrix = New_Matrix)

dfTabacum_MeanSD <- Calculate_Means_Std_Dev(matrix = Tabacum_Matrix)

#### 06 HISTOGRAM OF CODON USAGE FREQUENCIES #####

AmAcid <- rep(rownames(dfOld_MeanSD),3)
Species <- c((rep("Old N. benthamiana (Kazuza)",64)), (rep("New N. benthamiana",64)), (rep("N. tabacum (Kazuza)",64)))
Freq <- c(dfOld_MeanSD$`Frequency (Per 1000 Codons)`, dfNew_MeanSD$`Frequency (Per 1000 Codons)`, dfTabacum_MeanSD$`Frequency (Per 1000 Codons)`)
dfMeanSD <- data.frame(AmAcid,Species,Freq)

# Check distribution.

ggqqplot(dfOld_MeanSD$`Frequency (Per 1000 Codons)`)  
ggqqplot(dfNew_MeanSD$`Frequency (Per 1000 Codons)`)
ggqqplot(dfTabacum_MeanSD$`Frequency (Per 1000 Codons)`)

shapiro.test(dfOld_MeanSD$`Frequency (Per 1000 Codons)`) # p-value = 0.01822 Not normal.
shapiro.test(dfNew_MeanSD$`Frequency (Per 1000 Codons)`) # p-value = 0.09106 Normal.
shapiro.test(dfTabacum_MeanSD$`Frequency (Per 1000 Codons)`) # p-value = 0.1056 Normal.

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

#### 08 INTRA-SAMPLES KARLIN B PLOT ####

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
  ylab("MILC Distance from Consensus CU") 
#geom_jitter(width  = 0.05) + # Add noise since the dataset in small to avoid overplotting.
#geom_smooth(method = "lm", formula = y ~ x, fullrange = FALSE, level = 0.95) # 95% confidence interval
#geom_smooth(aes(group = 1), formula = y ~ poly(x,2), method = "lm", fullrange = TRUE, level = 0.95) 
#geom_curve()

# Again, not a lot of difference. As expected.
  
#### 09 ADD SECTION THAT COMPARES THE B() TO SUPPLEMENT INTRABPLOT ####

# Check if the B() values are normally distributed.

# Hypothesis: Data is normally distributed.
# If the p-value =< 0.05, hypothesis rejected, data is not normally distributed.

ggqqplot(B(New_CU)[,1])
shapiro.test(B(New_CU)) #  p-value = 6.007e-10 # Not normal.

ggqqplot(B(Old_CU)[,1])
shapiro.test(B(Old_CU)) #  p-value = 0.0001391 # Not normal.

ggqqplot(B(Tabacum_CU)[,1])
shapiro.test(B(Tabacum_CU)) #  p-value < 2.2e-16 # Not normal.

# Since, the data is not normally distributed, will conduct a non-parametric test to verify.
# A non-parametric alternative to ANOVA (which is parametic) is the Kruskal-Wallis Test.

kruskal.test(list(B(Old_CU),B(New_CU),B(Tabacum_CU)))

kruskal.test(list(B(Old_CU),B(Tabacum_CU)))

# Kruskal-Wallis rank sum test

#data:  list(B(Old_CU), B(New_CU), B(Tabacum_CU))
#Kruskal-Wallis chi-squared = 17.719, df = 2, p-value = 0.0001421

#### 11 ENC & MANN WHITNEY U TEST ####

New_ENC <- ENC(New_CU)
Old_ENC <- ENC(Old_CU)

# Check if the values are normally distributed using a Q-Q plot from the ggpubr package.

ggqqplot(New_ENC)
ggqqplot(Old_ENC)

shapiro.test(New_ENC) # W = 0.98266, p-value = 8.508e-05
# p-value < 0.05 
# distribution of data ist significantly different from normal distribution
# not normal

shapiro.test(Old_ENC) #W = 0.97915, p-value = 0.114

# p-value > 0.05 
# distribution of data not significantly different from normal distribution

# A non-parametric alternative to the t-test is Mann Whitney U test
wilcox.test(x = New_ENC, 
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

#### 12 CREATE AVERAGES & STD DEV DATA FRAME ####

Codon <- rep(rownames(dfOld_MeanSD),3)
AmAcid <- as.list((dfOld[order(dfOld$Codon),])[,1])
Sample <- c((rep("Old",64)), (rep("New",64)), (rep("Tabacum",64)))
Avrgs <- c(dfOld_MeanSD$`Frequency (Per 1000 Codons)`, dfNew_MeanSD$`Frequency (Per 1000 Codons)`, dfTabacum_MeanSD$`Frequency (Per 1000 Codons)`)
SD <- c(dfOld_MeanSD$Std_Deviation, dfNew_MeanSD$Std_Deviation, dfTabacum_MeanSD$Std_Deviation)

dfCodon_Mean_SD <- data.frame(Codon,AmAcid,Sample,Avrgs, SD)

#### ANOVA/ KRUSKAL WALLIS #####

shapiro.test(dfCodon_Mean_SD$Avrgs)
shapiro.test(dfOld_MeanSD$Average_Codon_Count)
shapiro.test(dfNew_MeanSD$Average_Codon_Count)
shapiro.test(dfTabacum_MeanSD$Average_Codon_Count)

# Data follows a normal distribution?

anova <- aov(Avrgs ~ as.factor(Sample), data = dfCodon_Mean_SD)

boxplot(Avrgs ~ Codon, data = dfCodon_Mean_SD)
kwtest <- kruskal.test(Avrgs ~ Sample, data = dfCodon_Mean_SD)

TukeyHSD(x = anova)


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
  scale_y_continuous(breaks = seq(0, 55, by = 5), limits = c(0,55)) + # adapted from: https://stackoverflow.com/questions/37950511/r-ggplot2-setting-tick-mark-interval
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

#### SUPPLEMENTARY TABLE: RSCU VALUES FOR NEW CU ####

# In addition to the ENC, we can also calculate the RSCU for the coding sequences in the new codon usage table.
# First format the coding sequences to a style easily read in to the uco() function in the seqinR package.

Formatted_CDS <- s2c((unlist(dfNew$Trimmed_CDS[1])))

# Calculate the RSCU values for each codon in the coding sequences.

New_RSCU <- uco(seq=Formatted_CDS,
                frame = 0,
                index = "rscu",
                as.data.frame = TRUE,
                NA.rscu = NA)

# View how many codons are used as expected, more than expected, and less than expected.

table(New_RSCU$RSCU == 1) # RSCU = 1 codon is used as expected by random usage 
table(New_RSCU$RSCU > 1) # RSCU > 1 codon used more frequently than random
table(New_RSCU$RSCU < 1) # RSCU < 1 codon used less frequently than random 

# Are these different from the RSCU for Old N.benthamiana CU? At the moment do not have the sequences for Old N.benthamiana CU. May have to import that (if time permits) and calculate the RSCU values.
# For now, can report the RSCU measures in a table (perhaps supplementary results?)

#### 15 REFERENCES ####