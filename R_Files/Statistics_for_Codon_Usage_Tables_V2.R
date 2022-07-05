# Statistical Anaylsis of Update Codon Usage in Nicotiana benthaminana
# 27 June 2022
# Anchitaa Ghag

#### PERSONAL NOTES ####

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

library("ape")
library("Biostrings")
library("coRdon")
library("ggplot2")
library("rentrez")
library("seqinr")
library("stringr")
library("TeachingDemos") # Using this to label all the outliers in a boxplot.
library("tidyverse")
library("XML") # Using this to parse HTML Kazuza's codon usage table
#install.packages("gridExtra")
library("gridExtra")
#install.packages("viridis")
library(viridis)

#### 02 DATA AQUISITION : IMPORT CU FROM PREVIOUS SCRIPT ####

dfNew <- read_csv("dfCodingSeqs.csv")[,2:14]
dfOld <- read_csv("dfKazuza.csv")[,2:5]
dfN.Tabacum <- read_csv("dfNTabacum.csv")[,2:5]

# Finally, calculate the most used codons (1 for each amino acid residue) from the data frame.

#dfKazuza.max <- group_by(dfKazuza, AmAcid) %>% 
#  filter(Number==max(Number)) %>% 
#  select(AmAcid, Codon)

#### MINOR NOTE ON AVERAGE CODON COUNTS REPORTED BY KAZUZA ####

# Note that the average number of codon counts per 1000 codons reported by the Kazuza database can differ slightly.
# For example, check the total number of codons (should be 1000) in Kazuza CU for N. benthamiana and N. tabacum.

sum(dfOld$X.1000) # Kazuza CU for N. benthamiana reports 999.99 instead of 1000

sum(dfN.Tabacum$X.1000) # Kazuza CU for N. tabacum reports 1000.04 instead of 1000

# While these are minor departures from the "1000 codon" value, these average codon counts will be recalculated to ensure that the frequencies reported are as exactly equal as possible.
# To do this, need to import the complete set of codon counts for each coding sequence used to generate Kazuza's CU tables.
# These can be found on the individual record in the database under the "List of codon usage for each CDS (format)" link.

# For N. benthamiana: http://www.kazusa.or.jp/codon/current/species/4100
# For N. tabacum: http://www.kazusa.or.jp/codon/current/species/4097

# For easier loading into R, i will be loading previously formatted files (using command line).

Old_CDS <- read.table("N_benthamiana_Codon_Counts_Only.txt",
                        header = TRUE)

N.tabacum_CDS <- read.table("N_tabacum_Codon_Counts_Only.txt",
                            header = TRUE)

#### GENERATE CODON USAGE TABLES ####

#https://bioconductor.riken.jp/packages/3.8/bioc/vignettes/coRdon/inst/doc/coRdon.html

# 1) First, convert the trimmed coding sequences to a DNAStringSet object.

cds <- DNAStringSet(dfNew$Trimmed_CDS)

# 2) Second, create a codon usage table object for both existing and current codon usage tables.
# NOTE! The columns in the list of codon counts imported from Kazuza are not sorted alphabetically! 
# These need to be reordered before inputting into the codonTable function or else counts will be assined to the wrong codon!

# Replace U w T
colnames(Old_CDS) <- str_replace_all(colnames(Old_CDS),"U","T") 
colnames(N.tabacum_CDS) <- str_replace_all(colnames(Old_CDS),"U","T") 
# Order the columns alphabetically.
Old_CDS_Ordered <- Old_CDS %>%
  select(order(colnames(Old_CDS)))
N.tabacum_CDS_Ordered <- N.tabacum_CDS %>%
  select(order(colnames(N.tabacum_CDS)))


Old_CU <- codonTable(Old_CDS_Ordered)
New_CU <- codonTable(cds)
Tabacum_CU <- codonTable(N.tabacum_CDS_Ordered)

# 3) Then, using the codonCounts() function from the codRon package create a matrix to view the counts per codon.

Old_Matrix <- codonCounts(Old_CU)
New_Matrix <- codonCounts(New_CU)
Tabacum_Matrix <- codonCounts(Tabacum_CU)


#getlen(Current_CU)
#row.names(CU) <- dfCodingSeqs$Titles
#(colMeans(CU)/colSums(CU))*1000

# Finally, create an consensus codon usage table based on the average empirical codon counts of coding sequences.

#Avg_Codon_Counts <- colMeans(New_Matrix)

#Freq_Per_1000 <- list()
#for (i in Avg_Codon_Counts){
#Value_Calc <- (i/sum(Avg_Codon_Counts))*1000
#Freq_Per_1000 <- append(Freq_Per_1000,Value_Calc )
#}

#dfTemp <- data.frame(colnames(New_Matrix),Avg_Codon_Counts,unlist(Freq_Per_1000))

# Reorder the columns to ensure they match.
# https://stackoverflow.com/questions/11977102/order-data-frame-rows-according-to-vector-with-specific-order

#dfTemp.sub <- dfTemp[match(dfExisting$Codon, dfTemp$colnames.CU.),]

# Add to final consensus data frame.

#dfCurrent <- data.frame(dfExisting$AmAcid,dfExisting$Codon,dfTemp.sub$Avg_Codon_Counts,dfTemp.sub$unlist.Freq_Per_1000.)
#colnames(dfCurrent) <- colnames(dfExisting)

#rm(cds,dfTemp,dfTemp.sub)


#### MEAN & STANDARD DEVIATIONS ####

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

##### Histogram #####

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

#### CHI SQUARE TEST ####

# Perform a chi-squared test to evaluate if there is a difference between the frequencies of both CU.

chisq.test(x = dfOld_MeanSD$`Frequency (Per 1000 Codons)`,
           y = dfNew_MeanSD$`Frequency (Per 1000 Codons)`)

# Pearson's Chi-squared test

# data:  dfOld_MeanSD$`Frequency (Per 1000 Codons)` and dfNew_MeanSD$`Frequency (Per 1000 Codons)`
# X-squared = 3904, df = 3843, p-value = 0.242

#######
# https://www.bioconductor.org/packages/devel/bioc/vignettes/coRdon/inst/doc/coRdon.html#calculate-cu-bias

# "Other CU statistics can be calculated in the same way as MILC(), using one of the functions: B(), MCB(), ENCprime(), ENC() or SCUO(). Note however, that when calculating ENC and SCUO, one doesn’t need to provide a subset of refe"
# Next, compare the CU bias for every coding sequence between the created CUT and Kazuza through visualizations.

library(ggplot2)


# Visualization of Codon Usage - Existing vs. Current Karlin B plot ####

# Code Adapted From: https://www.bioconductor.org/packages/devel/bioc/vignettes/coRdon/inst/doc/coRdon.html#calculate-cu-bias

# Use the intraBplot() function in the coRdon package and our two codonTable objects to plot a  plot of codon usage distances between existing vs. created codon tables.
#Intra-samples Karlin B plot

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
  labs(color="Codon Usage")
  #geom_jitter(width  = 0.05) + # Add noise since the dataset in small to avoid overplotting.
  #geom_smooth(method = "lm", formula = y ~ x, fullrange = FALSE, level = 0.95) # 95% confidence interval
  #geom_smooth(aes(group = 1), formula = y ~ poly(x,2), method = "lm", fullrange = TRUE, level = 0.95) 
  #geom_curve()

# The example dataset in the vignette has ~ 19, 000 + "points". There may just be undersampling/ less data points.


# This makes sense, since it is the same species, the codon usage should not differ greatly from the Old CU data set available.

### ADD SECTION THAT COMPARES THE B() TO SUPPLEMENT INTRABPLOT) ####


# Compare with N.tabacum 

set.seed(1111)
Sample_CU <- sample(Tabacum_CU,407)

intraBplot(x = New_CU, 
           y = Sample_CU, 
           names = c("Existing", "Consensus"), 
           variable = "MILC", 
           size = 3, 
           alpha = 1.0) + 
  ggtitle("Karlin B Plot of Existing (Kazuza) v.s. Consensus CU Distances")+
  xlab("MILC Distance from Existing CU") +
  ylab("MILC Distance from Consensus CU") #+ 
#geom_jitter(width  = 0.05) + # Add noise since the dataset in small to avoid overplotting.
#geom_smooth(method = "lm", formula = y ~ x, fullrange = FALSE, level = 0.95) # 95% confidence interval
#geom_smooth(aes(group = 1), formula = y ~ poly(x,2), method = "lm", fullrange = TRUE, level = 0.95) 
#geom_curve()

# Again, not a lot of difference. As expected.


#### RSCU ####

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

#### ENC & Mann Whitney U Test ####

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

table(enc == 20) # RSCU = 1 codon is used as expected by random usage 
table(enc == 61) # RSCU > 1 codon used more frequently than random
table(enc < 40) # RSCU < 1 codon used less frequently than random 

#### The Twenty Amino Acid Distributions ####

#dfCodon_Mean_SD <-

Codon <- rep(rownames(dfOld_MeanSD),3)
AmAcid <- as.list((dfOld[order(dfOld$Codon),])[,1])
Sample <- c((rep("Old",64)), (rep("New",64)), (rep("Tabacum",64)))
Avrgs <- c(dfOld_MeanSD$`Frequency (Per 1000 Codons)`, dfNew_MeanSD$`Frequency (Per 1000 Codons)`, dfTabacum_MeanSD$`Frequency (Per 1000 Codons)`)
SD <- c(dfOld_MeanSD$Std_Deviation, dfNew_MeanSD$Std_Deviation, dfTabacum_MeanSD$Std_Deviation)
dfCodon_Mean_SD <- data.frame(Codon,AmAcid,Sample,Avrgs, SD)

######

library(stringr)

# The number of times X codon appear/is counted in 1000 codons
# Can use one of two options. Codon counts "Number" column in the codon usage table or Frequency (per 1000 codons) "X.1000".
# Codon counts are not unifiorm. For organisms which have more CDS, the total codon count will be higher, and the total numer of counts will be higher for each codon. 
df1 <- dfExisting[,c(2,4)]
df2 <- dfCurrent[,c(2,4)]
df3 <- dfN.Tabacum[,c(2,4)]
dfAmAcid <- cbind(df1,df2,df3)


# Creates a data frame containing the codon count frequencies from each data frame.
dfAmAcid <- rbind(
  (subset(df1, AmAcid == Three_Letter_Form))[,c(2,4)],
  (subset(df2, AmAcid == Three_Letter_Form))[,c(2,4)],
  (subset(df3, AmAcid == Three_Letter_Form))[,c(2,4)]
)


##### Plotting Function #####

# Takes an argument "amino_acid" which is the name of the amino acid. Three data frames in the style of the Kazuza data base.

Plot_AmAcid <- function(amino_acid, three_letter, df) {
  
Three_Letter_Form <- three_letter

# Counts the number of times to repeat the character string.

#Times_To_Repeat <- as.numeric(table(df1$AmAcid == Three_Letter_Form)["TRUE"])

# Creates a data frame containing the codon count averages and std dev from original data frame.
dfAmAcid <- (subset(df, AmAcid == Three_Letter_Form))[,c(1,3,4,5)]

#Sample <- c(rep(x = "Existing", times = Times_To_Repeat),
            #rep(x = "Current", times = Times_To_Repeat),
           # rep(x = "N. tabacum", times = Times_To_Repeat))

#str1 <- (subset(df1, AmAcid == Three_Letter_Form))[,4]
#s1 <-  unlist(str1)

#str2 <- (subset(df2, AmAcid == Three_Letter_Form))[,4]
#s2 <-  unlist(str2)

#str3 <- (subset(df3, AmAcid == Three_Letter_Form))[,4]
#s3 <-  unlist(str3)

# Generate the plot using ggplot package functions.

Plot <- ggplot(data=dfAmAcid, 
       aes(x=Codon, y=Avrgs, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge(),alpha = 0.8) +
  ggtitle(amino_acid) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  xlab("Codon") +
  ylab("Average Codon Count") +
  #scale_fill_discrete(name = "Species", labels = c("Old N.benthamiana (Kazuza)", "New N.benthamiana ", "N.tabacum (Kazuza)"))+
  geom_errorbar(aes(ymin=Avrgs-SD, ymax=Avrgs+SD), width=.2, position=position_dodge(.9)) +
  ylim(NA,55) +
  scale_fill_viridis(discrete = TRUE, option = "viridis") +
  scale_color_manual(name= "Species",
                   labels = c("Old N.benthamiana (Kazuza)", "New N.benthamiana ", "N.tabacum (Kazuza)"))
                   
return(Plot)

}
                  


# Alanine - Ala

p1 <- Plot_AmAcid(amino_acid = "Alanine", 
            three_letter = "Ala", 
            df = dfCodon_Mean_SD)

# Arginine - Arg 

p2 <- Plot_AmAcid(amino_acid = "Arginine", 
            three_letter = "Arg", 
            df = dfCodon_Mean_SD)

# Asparagine - Asn

p3 <- Plot_AmAcid(amino_acid = "Asparagine", 
                  three_letter = "Asn", 
                  df = dfCodon_Mean_SD)

# Aspartic acid - Asp 

p4 <- Plot_AmAcid(amino_acid = "Aspartic acid", 
                  three_letter = "Asp", 
                  df = dfCodon_Mean_SD)

# Cysteine - Cys

p5 <- Plot_AmAcid(amino_acid = "Cysteine", 
                  three_letter = "Cys", 
                  df = dfCodon_Mean_SD)

# Glutamine - Gln

p6 <- Plot_AmAcid(amino_acid = "Glutamine", 
                  three_letter = "Gln", 
                  df = dfCodon_Mean_SD)

# Glutamic acid - Glu

p7 <- Plot_AmAcid(amino_acid = "Glutamic acid", 
                  three_letter = "Glu", 
                  df = dfCodon_Mean_SD)

# Glycine - Gly 

p8 <- Plot_AmAcid(amino_acid = "Glycine", 
                  three_letter = "Gly", 
                  df = dfCodon_Mean_SD)

# Histidine - His 

p9 <- Plot_AmAcid(amino_acid = "Histidine", 
                  three_letter = "His", 
                  df = dfCodon_Mean_SD)

# Isoleucine - Ile 

p10 <- Plot_AmAcid(amino_acid = "Isoleucine", 
                  three_letter = "Ile", 
                  df = dfCodon_Mean_SD)

# Leucine - Leu 

p11 <- Plot_AmAcid(amino_acid = "Leucine", 
                  three_letter = "Leu", 
                  df = dfCodon_Mean_SD)

# Lysine - Lys 

p12 <- Plot_AmAcid(amino_acid = "Lysine", 
                  three_letter = "Lys", 
                  df = dfCodon_Mean_SD)

# Methionine - Met

p13 <- Plot_AmAcid(amino_acid = "Methionine", 
                  three_letter = "Met", 
                  df = dfCodon_Mean_SD)

# Phenylalanine - Phe 

p14 <- Plot_AmAcid(amino_acid = "Phenylalanine", 
                  three_letter = "Phe", 
                  df = dfCodon_Mean_SD)

# Proline - Pro 

p15 <- Plot_AmAcid(amino_acid = "Proline", 
                  three_letter = "Pro", 
                  df = dfCodon_Mean_SD)

# Serine - Ser

p16 <- Plot_AmAcid(amino_acid = "Serine", 
                  three_letter = "Ser", 
                  df = dfCodon_Mean_SD)

# Threonine - Thr 

p17 <- Plot_AmAcid(amino_acid = "Threonine", 
                  three_letter = "Thr", 
                  df = dfCodon_Mean_SD)

# Tryptophan - Trp 

p18 <- Plot_AmAcid(amino_acid = "Tryptophan", 
                  three_letter = "Trp", 
                  df = dfCodon_Mean_SD)

# Tyrosine - Tyr 

p19 <- Plot_AmAcid(amino_acid = "Tyrosine", 
                  three_letter = "Tyr", 
                  df = dfCodon_Mean_SD)

# Valine - Val 

p20 <- Plot_AmAcid(amino_acid = "Valine", 
                  three_letter = "Val", 
                  df = dfCodon_Mean_SD)

#### grid ####

# https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html

grid.arrange(p1, p2, p3, p4, nrow = 2)
grid.arrange( p5, p6, p7, p8, nrow = 2)
grid.arrange( p9, p10, p11, p12, nrow = 2)
grid.arrange( p13, p14, p15, p16, nrow = 2)
grid.arrange( p17, p18, p19, p20, nrow = 2)
