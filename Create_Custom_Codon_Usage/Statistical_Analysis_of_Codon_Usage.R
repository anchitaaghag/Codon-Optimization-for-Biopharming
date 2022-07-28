# Codon Optimization for Biopharming: Statistical Analysis of Updated Codon Usage in Nicotiana benthaminana
# Anchitaa Ghag

setwd("/Users/anchitaa/Major_Research_Project_2022/06_Code/Data/")

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
#install.packages("seqinr")
#install.packages("tidyverse")
#install.packages("viridis")
#install.packages("XML")

library("Biostrings")
library("coRdon")
library("ggplot2")
library("gridExtra")
library("seqinr")
library("tidyverse")
library("viridis")
library("XML")

#### 02 DATA AQUISITION : IMPORT EXISTING CU FROM KAZUZA ####

# The following code was adapted from: https://stackoverflow.com/questions/24546312/vector-of-most-used-codons-from-table-of-codon-usage-in-r

# The existing codon usage table from the Kazuza database can be imported from: 
# https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4100&aa=1&style=GCG

# It can also be imported from the "Kazuza_Codon_Usage.txt" file using the following lines of code. File originally created on 24 June 2022.

# Kazuza <- read_table("Kazuza_Codon_Usage.txt")

# First, use the XML package's htmlParse() function to "read" an HTML file and generate an HTML/XMLInternalDocument class object.

Kazuza <- htmlParse('http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4100&aa=1&style=GCG')

Tabacum <- htmlParse('http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4097&aa=1&style=GCG')

# Next, read this object and convert to a dataframe.

dfOld <- read.table(text=xpathSApply(Kazuza, "//pre", xmlValue), 
                       header=TRUE, 
                       fill=TRUE)

dfTabacum <- read.table(text=xpathSApply(Tabacum, "//pre", xmlValue), 
                        header=TRUE, 
                        fill=TRUE)

# List of the codon usage from each CDS that makes up Kazuza CUT is available at: http://www.kazusa.or.jp/codon/current/species/4100

# KazuzaCDS <- read.table("/Users/anchitaa/Major_Research_Project_2022/06_Code/KCC_No_Carot.txt",
# header = TRUE)

# Remove all objects no longer needed from the environment.

rm(Kazuza, Tabacum)

#### 03 DATA AQUISITION : IMPORT CU FROM PREVIOUS SCRIPT ####

# Import the codon usage (to be compared to Kazuza codon usages) from the previous R file "Build_Codon_Usage_Table.R" or from another source.

dfNew <- read_csv(file = "Updated_Codon_Usage_Information.csv")

#### 04 VERIFY CODON COUNTS REPORTED BY KAZUZA ####

# Note that the average number of codon counts per 1000 codons reported by the Kazuza database can differ slightly.
# For example, check the total number of codons (should be 1000) in Kazuza CU for N. benthamiana and N. tabacum.

# Kazuza CU for N. benthamiana reports 999.99 instead of 1000

sum(dfOld$X.1000) 

# Kazuza CU for N. tabacum reports 1000.04 instead of 1000

sum(dfTabacum$X.1000)

# While these are minor departures from the "1000 codon" value, these average codon counts will be recalculated to ensure that the frequencies reported are as exactly equal as possible.
# To do this, need to import the complete set of codon counts for each coding sequence used to generate Kazuza's CU tables.
# These can be found on the individual record in the database under the "List of codon usage for each CDS (format)" link.

# For N. benthamiana: http://www.kazusa.or.jp/codon/current/species/4100
# For N. tabacum: http://www.kazusa.or.jp/codon/current/species/4097

# For faster loading into R, I will be loading previously formatted files (using command line).

Old_CDS <- read.table("N_benthamiana_Codon_Counts_Only.txt",
                        header = TRUE)

Tabacum_CDS <- read.table("N_tabacum_Codon_Counts_Only.txt",
                            header = TRUE)

#### 05 GENERATE CODON USAGE TABLES ####

#https://bioconductor.riken.jp/packages/3.8/bioc/vignettes/coRdon/inst/doc/coRdon.html

# First, convert the trimmed coding sequences to a DNAStringSet object.

New_CDS <- DNAStringSet(dfNew$Trimmed_CDS)

# Second, create a codon usage table object for both existing and current codon usage tables.
# NOTE! The columns in the list of codon counts imported from Kazuza are not sorted alphabetically! 
# These need to be reordered before inputting into the codonTable function or else counts will be assined to the wrong codon!

# Replace "U" with "T" in the sequences. 

colnames(Old_CDS) <- str_replace_all(colnames(Old_CDS),"U","T") 
colnames(Tabacum_CDS) <- str_replace_all(colnames(Tabacum_CDS),"U","T") 

# Order the columns alphabetically.

Old_CDS_Ordered <- Old_CDS %>%
  select(order(colnames(Old_CDS)))

Tabacum_CDS_Ordered <- Tabacum_CDS %>%
  select(order(colnames(Tabacum_CDS)))

# 

Old_CU <- codonTable(Old_CDS_Ordered)
New_CU <- codonTable(New_CDS)
Tabacum_CU <- codonTable(Tabacum_CDS_Ordered)

# Then, using the codonCounts() function from the codRon package create a matrix to view the counts per codon.

Old_Matrix <- codonCounts(Old_CU)
New_Matrix <- codonCounts(New_CU)
Tabacum_Matrix <- codonCounts(Tabacum_CU)

rm(New_CDS, Old_CDS, Tabacum_CDS)

#### 06 MEAN & STANDARD DEVIATIONS ####

# In this section, we will calculate the average (i.e. mean) codon counts per codon and the standard deviation.
# This step will be performed for each of the 64 amino acids across the three samples (i.e. data frames).

dfOld_MeanSD <- Calculate_Means_Std_Dev(matrix = Old_Matrix)

dfNew_MeanSD <- Calculate_Means_Std_Dev(matrix = New_Matrix)

dfTabacum_MeanSD <- Calculate_Means_Std_Dev(matrix = Tabacum_Matrix)

rm (Old_CDS, Tabacum_CDS, New_Matrix, Old_Matrix, Tabacum_Matrix, Old_CDS_Ordered, Tabacum_CDS_Ordered)

#### 07 CREATE PLOT DATA FRAME ####

Codon <- rep(rownames(dfOld_MeanSD),3)
AmAcid <- (dfOld[order(dfOld$Codon),])[,1]
Sample <- c((rep("Old",64)), (rep("New",64)), (rep("Tabacum",64)))
Avrgs <- c(dfOld_MeanSD$`Frequency (Per 1000 Codons)`, dfNew_MeanSD$`Frequency (Per 1000 Codons)`, dfTabacum_MeanSD$`Frequency (Per 1000 Codons)`)
SD <- c(dfOld_MeanSD$Std_Deviation, dfNew_MeanSD$Std_Deviation, dfTabacum_MeanSD$Std_Deviation)

dfCodon_Mean_SD <- data.frame(Codon,AmAcid,Sample,Avrgs, SD)

rm (Codon, AmAcid, Sample, Avrgs, SD)

#### 08 STATISTICAL TEST : CHI SQUARE TEST ####

# Perform a chi-squared test to evaluate if there is a difference between the frequencies of both CU.

chisq.test(x = dfOld_MeanSD$`Frequency (Per 1000 Codons)`,
           y = dfNew_MeanSD$`Frequency (Per 1000 Codons)`)

# Pearson's Chi-squared test

# data:  dfOld_MeanSD$`Frequency (Per 1000 Codons)` and dfNew_MeanSD$`Frequency (Per 1000 Codons)`
# X-squared = 3904, df = 3843, p-value = 0.242

#### 09 STATISTICAL TEST : CODON USAGE FREQUENCIES & KRUSKAL-WALLIS RANK SUM TEST #####

# Check distribution.

shapiro.test(dfOld_MeanSD$`Frequency (Per 1000 Codons)`) # p-value = 0.01822 Not normal.
shapiro.test(dfNew_MeanSD$`Frequency (Per 1000 Codons)`) # p-value =  0.03754 Not normal.
shapiro.test(dfTabacum_MeanSD$`Frequency (Per 1000 Codons)`) # p-value = 0.1056 Normal.

# Proceed with the non-parametric Kruskal-Wallis Rank Sum test to evaluate if there are significant differences between the frequencies of the three codon usages.

kruskal.test(x = list(dfOld_MeanSD$`Frequency (Per 1000 Codons)`,
                      dfNew_MeanSD$`Frequency (Per 1000 Codons)`,
                      dfTabacum_MeanSD$`Frequency (Per 1000 Codons)`))

# Kruskal-Wallis chi-squared = 0.013975, df = 2, p-value = 0.993. Null hypothesis is not rejected.

#### 10 PLOT : HISTOGRAM OF CODON USAGE FREQUENCIES #####

AmAcid <- rep(rownames(dfOld_MeanSD),3)
Species <- c((rep("Old N. benthamiana (Kazuza)",64)), (rep("New N. benthamiana",64)), (rep("N. tabacum (Kazuza)",64)))
Freq <- c(dfOld_MeanSD$`Frequency (Per 1000 Codons)`, dfNew_MeanSD$`Frequency (Per 1000 Codons)`, dfTabacum_MeanSD$`Frequency (Per 1000 Codons)`)
dfMeanSD <- data.frame(AmAcid,Species,Freq)

# Compare the codon usage frequency using a histogram.

ggplot(data=dfMeanSD, 
      aes(x=AmAcid, y=Freq, fill=Species)) +
  geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8), show.legend = TRUE) + # https://www.learnbyexample.org/r-bar-plot-ggplot2/
  ggtitle("Codon Usage Frequency Among 3 Nicotiana Codon Usages") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Codons") +
  ylab("Frequency (Per 1000 Codons)") +
  scale_x_discrete(guide = guide_axis(angle = 90)) + #https://stackoverflow.com/questions/1330989/rotating-and-spacing-axis-labels-in-ggplot2
  theme(legend.position = c(0.89, 0.90),
        legend.background = element_rect(fill = "white", color = "black")) + #https://datavizpyr.com/how-to-place-legend-inside-the-plot-with-ggplot2/ 
  scale_y_continuous(limits = c(0,40), expand = c(0, 0)) + # https://stackoverflow.com/questions/22945651/remove-space-between-plotted-data-and-the-axes
  scale_fill_viridis(discrete = TRUE, option = "viridis") 
  #scale_fill_manual(values=c("#648FFF", "#DC267F", "#FFB000")) +
  #coord_flip()

rm(AmAcid,Species,Freq)

#### 10 ALTERNATIVE PLOT : CODON USAGE FREQUENCIES #####

# https://stackoverflow.com/questions/61200151/how-to-adjust-relative-transparency-of-ggplot2-points
# https://stackoverflow.com/questions/32423167/ggplot-xy-scatter-how-to-change-alpha-transparency-for-select-points
# https://stackoverflow.com/questions/13035295/overlay-bar-graphs-in-ggplot2
# Inspiration : https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003833

library(ggplot2)
library(reshape)

x <- rownames(dfOld_MeanSD)
y1 <- dfOld_MeanSD$`Frequency (Per 1000 Codons)`
y2 <- dfNew_MeanSD$`Frequency (Per 1000 Codons)`
y3 <- dfTabacum_MeanSD$`Frequency (Per 1000 Codons)`

to_plot <- data.frame(x=x,y1=y1,y2=y2, y3=y3)

melted<-melt(to_plot, id="x")

#### PLOT 1
alpha_vector = rep(0.3, 192)
alpha_vector[c(1:64)] = 1
melted$alpha = alpha_vector

plot1 <- ggplot(melted,aes(x=x,y=value,fill=variable)) + 
  geom_bar(stat="identity", position = "identity", aes(alpha=alpha), show.legend = FALSE) +
  ggtitle(expression('Codon Usage Frequency Among'~italic("Nicotiana")~'Species')) + #https://stackoverflow.com/questions/32555531/how-to-italicize-part-one-or-two-words-of-an-axis-title
  scale_alpha_identity(alpha) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) +
  xlab("") + # https://stackoverflow.com/questions/35090883/remove-all-of-x-axis-labels-in-ggplot
  ylab("Frequency (Per 1000)") +
  scale_x_discrete(labels = NULL, breaks = NULL) + # https://stackoverflow.com/questions/35090883/remove-all-of-x-axis-labels-in-ggplot
  scale_y_continuous(limits = c(0,40), expand = c(0, 0)) + # https://stackoverflow.com/questions/22945651/remove-space-between-plotted-data-and-the-axes
  annotate(geom = "text", label = "N. benthamiana (Kazusa)", x = 57, y = 38) +
  scale_fill_manual(values=c("#FDE725", "#000000", "#000000")) 

#### PLOT 2
alpha_vector = rep(0.3, 192)
alpha_vector[c(65:128)] = 1
melted$alpha = alpha_vector

plot2 <- ggplot(melted,aes(x=x,y=value,fill=variable)) + 
  geom_bar(stat="identity", position = "identity", aes(alpha=alpha), show.legend = FALSE) +
  scale_alpha_identity(alpha) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) +
  xlab("") + # https://stackoverflow.com/questions/35090883/remove-all-of-x-axis-labels-in-ggplot
  ylab("Frequency (Per 1000)") +
  scale_x_discrete(labels = NULL, breaks = NULL) + # https://stackoverflow.com/questions/35090883/remove-all-of-x-axis-labels-in-ggplot
  scale_y_continuous(limits = c(0,40), expand = c(0, 0)) + # https://stackoverflow.com/questions/22945651/remove-space-between-plotted-data-and-the-axes
  annotate(geom = "text", label = "N. benthamiana (Updated)", x = 57, y = 38) +
  scale_fill_manual(values=c("#000000", "#21918C", "#000000")) 

#### PLOT 3 
alpha_vector = rep(0.3, 192)
alpha_vector[c(129:192)] = 1
melted$alpha = alpha_vector

plot3 <- ggplot(melted,aes(x=x,y=value,fill=variable)) + 
  geom_bar(stat="identity", position = "identity", aes(alpha=alpha), show.legend = FALSE) +
  scale_alpha_identity(alpha) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) +
  xlab("Codons") +
  ylab("Frequency (Per 1000)") +
  scale_x_discrete(guide = guide_axis(angle = 90)) + #https://stackoverflow.com/questions/1330989/rotating-and-spacing-axis-labels-in-ggplot2
  scale_y_continuous(limits = c(0,40), expand = c(0, 0)) + # https://stackoverflow.com/questions/22945651/remove-space-between-plotted-data-and-the-axes
  annotate(geom = "text", label = "N. tabacum (Kazusa)", x = 57, y = 38) +
  scale_fill_manual(values=c("#000000", "#000000", "#440154")) 

grid.arrange(plot1, plot2, plot3, nrow=3)

#### 11 STATISTICAL TEST : CODON USAGE BIAS/MILC & KRUSKAL-WALLIS RANK SUM TEST ####

# COMPARE IF CODON USAGE BIAS IS THE DIFFERENT

# A minimum of 80 (or more codons) is required to conduct MILC analysis.
# Remove sequences less than 80 codons or (80*3 =) 240 bases in sequence length.

# Check if the MILC() values are normally distributed.

# Hypothesis: Data is normally distributed.
# If the p-value =< 0.05, hypothesis rejected, data is not normally distributed.

# Create three new objects to hold the MILC values for each group. The sequences will be filtered to exclude sequences that are less than 80 codons in length (as documented in literature).
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1199580/

New_MILC <- MILC(cTobject = New_CU, 
                  filtering = "hard", 
                  len.threshold = 80)

Old_MILC <- MILC(cTobject = Old_CU, 
                  filtering = "hard", 
                  len.threshold = 80)

Tabacum_MILC <- MILC(cTobject = Tabacum_CU, 
                  filtering = "hard", 
                  len.threshold = 80) 



shapiro.test(New_MILC) #  p-value = 2.35e-11 # Not normal.

shapiro.test(Old_MILC) #  p-value = 3.819e-08 # Not normal.

shapiro.test(Tabacum_MILC) #  p-value < 2.2e-16 # Not normal.


# Since, the data is not normally distributed, will conduct a non-parametric test to verify.
# A non-parametric alternative to ANOVA (which is parametic) is the Kruskal-Wallis Test.

kruskal.test(list(New_MILC,Old_MILC,Tabacum_MILC))

# Kruskal-Wallis chi-squared = 6.9591, df = 2, p-value = 0.03082
# Since the Kruskal–Wallis test did not indicate significant differences, a post-hoc analysis was not performed.


# p-value < 0.05 
# distribution of data ist significantly different from normal distribution
# not normal


# p-value > 0.05 
# distribution of data not significantly different from normal distribution
# normal


# p-value < 0.05 significant
# p-value >= 0.05 not significant
# H0 <- mean 1 - mean 2 = 0
# HA <- mean 1 - mean 2 =/= 0

rm(New_MILC,Old_MILC,Tabacum_MILC)

#### 12 STATISTICAL TEST : AVERAGE CODON COUNTS & KRUSKAL-WALLIS RANK SUM TEST #####

shapiro.test(dfOld_MeanSD$Average_Codon_Count) # p-value = 0.01822. Not normal.
shapiro.test(dfNew_MeanSD$Average_Codon_Count) # p-value = 0.03754. Not normal.
shapiro.test(dfTabacum_MeanSD$Average_Codon_Count) # p-value = 0.1056. Normal.

# Data does not follow a normal distribution. Use non-parametric test.

kruskal.test(list(dfNew_MeanSD$Average_Codon_Count,
                  dfOld_MeanSD$Average_Codon_Count,
                  dfTabacum_MeanSD$Average_Codon_Count))

# Kruskal-Wallis chi-squared = 0.86837, df = 2, p-value = 0.6478
# Since the Kruskal–Wallis test did not indicate significant differences, a post-hoc analysis was not performed.

#### 13 PLOT : AVERAGE COUNTS PER AMINO ACID #####

# For each of the 20 amino acids, create a plot comparing the average codon counts and standard deviations using the Plot_AmAcid() function. 

# Alanine - Ala

p1 <- Plot_AmAcid(amino_acid = "Alanine", 
                  df = dfCodon_Mean_SD)

# Arginine - Arg 

p2 <- Plot_AmAcid(amino_acid = "Arginine", 
                  df = dfCodon_Mean_SD)

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

# ARRANGE AMINO ACID PLOTS IN A GRID 

# https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html

grid.arrange(p1, p2, p3, p4, nrow = 2)
grid.arrange( p5, p6, p7, p8, nrow = 2)
grid.arrange( p9, p10, p11, p12, nrow = 2)
grid.arrange( p13, p14, p15, p16, nrow = 2)
grid.arrange( p17, p18, p19, p20, nrow = 2)

# The most used codons (one for each amino acid residue) from the three groups appears to be the same.

rm (p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20)

#### 14 GC% CONTENT ####

# In this last section, the overall G-C content of the coding sequences will be calculated.
# First, concatenate all coding sequences to one character string.

# https://stackoverflow.com/questions/9314328/how-to-collapse-a-list-of-characters-into-a-single-string-in-r

One_Long_Sequence <- toString(dfNew$Trimmed_CDS) 

# For each base, count the number of times that base is encountered and then convert to a percentage.

# https://datacarpentry.org/semester-biology/exercises/Loops-stringr-R/

# A

A_Percentage <- (str_count(One_Long_Sequence, "A")) / str_length(One_Long_Sequence) * 100

# C

C_Percentage <- (str_count(One_Long_Sequence, "C")) / str_length(One_Long_Sequence) * 100

# G

G_Percentage <- (str_count(One_Long_Sequence, "G")) / str_length(One_Long_Sequence) * 100

# T

T_Percentage <- (str_count(One_Long_Sequence, "T")) / str_length(One_Long_Sequence) * 100

# GC content

GC_Content <- (str_count(One_Long_Sequence, "G") + 
                 str_count(One_Long_Sequence, "C")) / str_length(One_Long_Sequence) * 100 

# Output the results.

paste("The ACGT percentages are", 
      round(A_Percentage,2), "%,", 
      round(C_Percentage,2), "%,", 
      round(G_Percentage,2), "%,", 
      round(T_Percentage,2), "%")

# The ACGT percentages are 29.08 %, 19.81 %, 23.91 %, 27.04 %

paste("The overall GC content is" , round(GC_Content,2), "%")

# The overall GC content is 43.72 %

rm(One_Long_Sequence)

#### 15 REFERENCES ####


#
