# Codon Optimization for Biopharming: Statistical Analysis of Updated Codon Usage in Nicotiana benthaminana
# Anchitaa Ghag

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
#install.packages("lemon")
#install.packages("seqinr")
#install.packages("tidyverse")
#install.packages("viridis")
#install.packages("XML")

library("Biostrings")
library("coRdon")
library("ggplot2")
library("gridExtra")
library("lemon")
library("seqinr")
library("tidyverse")
library("viridis")
library("XML")

# Load the reverse translate function from the R file.

source("Functions_For_Stat.R")

#### 02 DATA AQUISITION : IMPORT EXISTING CU FROM KAZUZA ####

# The existing codon usage table can be directly imported from the Kazuza database.

# First, use the XML package's htmlParse() function to "read" an HTML file and generate an HTML/XMLInternalDocument class object.
# The following two lines of code was adapted from: https://stackoverflow.com/questions/24546312/vector-of-most-used-codons-from-table-of-codon-usage-in-r. (shadow, 2014).

Kazuza <- htmlParse('http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4100&aa=1&style=GCG')

Tabacum <- htmlParse('http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4097&aa=1&style=GCG')

# Next, read this object and convert to a dataframe.
# The following lines of code was adapted from: https://stackoverflow.com/questions/24546312/vector-of-most-used-codons-from-table-of-codon-usage-in-r. (shadow, 2014).

dfOld <- read.table(text=xpathSApply(Kazuza, "//pre", xmlValue), 
                       header=TRUE, 
                       fill=TRUE)

dfTabacum <- read.table(text=xpathSApply(Tabacum, "//pre", xmlValue), 
                        header=TRUE, 
                        fill=TRUE)

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

# For faster loading into R, I will be loading previously formatted files (using command line). These files are also available in the Kazuza_Codon_Usage_Data folder.

Old_CDS <- read.table("N_benthamiana_Codon_Counts_Only.txt",
                        header = TRUE)

Tabacum_CDS <- read.table("N_tabacum_Codon_Counts_Only.txt",
                            header = TRUE)

#### 05 GENERATE CODON USAGE TABLES ####

# First, convert the trimmed coding sequences to a DNAStringSet object.

New_CDS <- DNAStringSet(dfNew$Trimmed_CDS)

# Second, create a codon usage table object for both existing and current codon usage tables.
# Note the columns in the list of codon counts imported from Kazuza are not sorted alphabetically. 
# These need to be reordered before inputting into the codonTable function or else counts will be assigned to the wrong codon.

# Replace "U" with "T" in the sequences. 

colnames(Old_CDS) <- str_replace_all(colnames(Old_CDS),"U","T") 

colnames(Tabacum_CDS) <- str_replace_all(colnames(Tabacum_CDS),"U","T") 

# Order the columns alphabetically.

Old_CDS_Ordered <- Old_CDS %>%
  select(order(colnames(Old_CDS)))

Tabacum_CDS_Ordered <- Tabacum_CDS %>%
  select(order(colnames(Tabacum_CDS)))

# The codonTable() and codonCounts() functions were learned from: https://bioconductor.riken.jp/packages/3.8/bioc/vignettes/coRdon/inst/doc/coRdon.html. (Elek, 2019).

Old_CU <- codonTable(Old_CDS_Ordered)

New_CU <- codonTable(New_CDS)

Tabacum_CU <- codonTable(Tabacum_CDS_Ordered)

# Then, using the codonCounts() function from the codRon package create a matrix to view the counts per codon.

Old_Matrix <- codonCounts(Old_CU)

New_Matrix <- codonCounts(New_CU)

Tabacum_Matrix <- codonCounts(Tabacum_CU)

# Remove all objects no longer needed from the environment.

rm(New_CDS, Old_CDS, Tabacum_CDS)

#### 06 MEAN & STANDARD DEVIATIONS ####

# In this section, we will calculate the average (i.e. mean) codon counts per codon and the standard deviation.
# This step will be performed for each of the 64 amino acids across the three samples (i.e. data frames).

dfOld_MeanSD <- Calculate_Means_Std_Dev(matrix = Old_Matrix)

dfNew_MeanSD <- Calculate_Means_Std_Dev(matrix = New_Matrix)

dfTabacum_MeanSD <- Calculate_Means_Std_Dev(matrix = Tabacum_Matrix)

# Remove all objects no longer needed from the environment.

rm (Old_CDS, Tabacum_CDS, New_Matrix, Old_Matrix, Tabacum_Matrix, Old_CDS_Ordered, Tabacum_CDS_Ordered)

#### 07 CREATE PLOT DATA FRAME ####

# To create a data frame with the relevant information to plot figures.
# First, create a list of the names of the species/codon usages.

Codon <- rep(rownames(dfOld_MeanSD),3)

# Then, create a list of the amino acids from the codon usages.

AmAcid <- (dfOld[order(dfOld$Codon),])[,1]

# Third, create a list of the species names.

Species <- c((rep("N. benthamiana Codon Usage (Kazusa)",64)), (rep("N. benthamiana Codon Usage (Updated)",64)), (rep("N. tabacum Codon Usage (Kazusa)",64)))

# Forth, create a list of the averages.

Avrgs <- c(dfOld_MeanSD$`Frequency (Per 1000 Codons)`, dfNew_MeanSD$`Frequency (Per 1000 Codons)`, dfTabacum_MeanSD$`Frequency (Per 1000 Codons)`)

# Fifth, create a list of the standard deviations.

SD <- c(dfOld_MeanSD$Std_Deviation, dfNew_MeanSD$Std_Deviation, dfTabacum_MeanSD$Std_Deviation)

# Finally, put together the lists into one data frame.

dfCodon_Mean_SD <- data.frame(Codon,AmAcid,Species,Avrgs, SD)

# Remove all objects no longer needed from the environment.

rm (Codon, AmAcid, Species, Avrgs, SD)

#### 08 CHI SQUARE TEST ####

# Perform a chi-squared test to evaluate if there is a difference between the frequencies of both CU.

chisq.test(x = dfOld_MeanSD$`Frequency (Per 1000 Codons)`,
           y = dfNew_MeanSD$`Frequency (Per 1000 Codons)`)

# Pearson's Chi-squared test

# data:  dfOld_MeanSD$`Frequency (Per 1000 Codons)` and dfNew_MeanSD$`Frequency (Per 1000 Codons)`
# X-squared = 3904, df = 3843, p-value = 0.242

#### 09 CODON USAGE FREQUENCIES & KRUSKAL-WALLIS RANK SUM TEST #####

# Check if the distribution is normal.

shapiro.test(dfOld_MeanSD$`Frequency (Per 1000 Codons)`) # p-value = 0.01822 Not normal.

shapiro.test(dfNew_MeanSD$`Frequency (Per 1000 Codons)`) # p-value =  0.03754 Not normal.

shapiro.test(dfTabacum_MeanSD$`Frequency (Per 1000 Codons)`) # p-value = 0.1056 Normal.

# Since, the data is not normally distributed, will conduct a non-parametric test to verify.
# A non-parametric alternative to ANOVA (which is parametic) is the Kruskal-Wallis Test.
# Proceed with the non-parametric Kruskal-Wallis Rank Sum test to evaluate if there are significant differences between the frequencies of the three codon usages.

kruskal.test(x = list(dfOld_MeanSD$`Frequency (Per 1000 Codons)`,
                      dfNew_MeanSD$`Frequency (Per 1000 Codons)`,
                      dfTabacum_MeanSD$`Frequency (Per 1000 Codons)`))

# Kruskal-Wallis chi-squared = 0.013975, df = 2, p-value = 0.993. Null hypothesis is not rejected.

#### 10 BAR PLOT OF CODON USAGE FREQUENCIES #####

# To create a data frame with the relevant information to plot figures.
# First, create a list of the amino acids.

AmAcid <- rep(rownames(dfOld_MeanSD),3)

# Then, create a list of the species names.

Species <- c((rep("N. benthamiana (Kazuza)",64)), (rep("N. benthamiana (Updated)",64)), (rep("N. tabacum (Kazuza)",64)))

# Next, create a list of the codon usage frequencies.

Freq <- c(dfOld_MeanSD$`Frequency (Per 1000 Codons)`, dfNew_MeanSD$`Frequency (Per 1000 Codons)`, dfTabacum_MeanSD$`Frequency (Per 1000 Codons)`)

# Finally, put together the lists into one data frame.

dfMeanSD <- data.frame(AmAcid,Species,Freq)

# View summary information.

summary(dfMeanSD$Freq)

# Since the 3rd IQR is 22.1716, most likely the top 10 codon frequencies are above this value.

# Filter the data frame to retain frequencies above 24.

dfTop_Freq <- dfMeanSD %>%
  filter(Freq > 24)

# Are there any missing columns?

table(dfTop_Freq$AmAcid)

# Filter out columns that do not have 3 data points (i.e. 3 groups).

dfTop_10_Freq <- dfTop_Freq %>%
  filter(!AmAcid == "ATG") %>%
  filter(!AmAcid == "GGT") %>%
  filter(!AmAcid == "TTG") %>%
  filter(!AmAcid == "TTT")

# Compare the codon usage frequency using a histogram.

ggplot(data=dfTop_10_Freq, 
       aes(x=AmAcid, y=Freq, fill=Species)) +
  # The following line of code was adapted from: https://www.learnbyexample.org/r-bar-plot-ggplot2/. (LearnByExample, 2022).
  geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8), show.legend = TRUE) + 
  ggtitle("Top 10 Codon Frequencies Among 3 Nicotiana Codon Usages") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
  xlab("Codons") +
  ylab("Frequency (Per 1000 Codons)") +
  # The following line of code was adapted from: https://stackoverflow.com/questions/1330989/rotating-and-spacing-axis-labels-in-ggplot2. (jan-glx, 2020).
  scale_x_discrete(guide = guide_axis(angle = 90)) + 
  # The following line of code was adapted from: https://stackoverflow.com/questions/22945651/remove-space-between-plotted-data-and-the-axes. (Jaap, 2014).
  scale_y_continuous(limits = c(0,40), expand = c(0, 0)) + 
  scale_fill_manual(values=c("#ed7953","#0d0887", "#9c179e")) 

# Remove all objects no longer needed from the environment.

rm(AmAcid,Species,Freq)

#### 11 AVERAGE CODON COUNTS & KRUSKAL-WALLIS RANK SUM TEST #####

shapiro.test(dfOld_MeanSD$Average_Codon_Count) # p-value = 0.01822. Not normal.

shapiro.test(dfNew_MeanSD$Average_Codon_Count) # p-value = 0.03754. Not normal.

shapiro.test(dfTabacum_MeanSD$Average_Codon_Count) # p-value = 0.1056. Normal.

# Data does not follow a normal distribution. Use non-parametric test.

kruskal.test(list(dfNew_MeanSD$Average_Codon_Count,
                  dfOld_MeanSD$Average_Codon_Count,
                  dfTabacum_MeanSD$Average_Codon_Count))

# Kruskal-Wallis chi-squared = 0.86837, df = 2, p-value = 0.6478
# Since the Kruskal–Wallis test did not indicate significant differences, a post-hoc analysis was not performed.

#### 12 PLOT AMINO ACIDS #####

# For each of the 20 amino acids, create a plot comparing the average codon counts and standard deviations using the Plot_AmAcid() function. 

# Alanine - Ala

p1 <- Plot_AmAcid(amino_acid = "Alanine", 
                  df = dfCodon_Mean_SD) +
  theme(legend.position="hidden")

# Arginine - Arg 

p2 <- Plot_AmAcid(amino_acid = "Arginine", 
                  df = dfCodon_Mean_SD) +
  theme(legend.position="hidden") +
  scale_y_continuous(labels = NULL, breaks = NULL)

# Asparagine - Asn

p3 <- Plot_AmAcid(amino_acid = "Asparagine", 
                  df = dfCodon_Mean_SD) +
  theme(legend.position="hidden") +
  scale_y_continuous(labels = NULL, breaks = NULL)

# Aspartic acid - Asp 

p4 <- Plot_AmAcid(amino_acid = "Aspartic acid", 
                  df = dfCodon_Mean_SD) +
  theme(legend.position="hidden") +
  scale_y_continuous(labels = NULL, breaks = NULL)

# Cysteine - Cys

p5 <- Plot_AmAcid(amino_acid = "Cysteine", 
                  df = dfCodon_Mean_SD) +
  theme(legend.position="hidden") +
  scale_y_continuous(labels = NULL, breaks = NULL)

# Glutamine - Gln

p6 <- Plot_AmAcid(amino_acid = "Glutamine", 
                  df = dfCodon_Mean_SD) +
  theme(legend.position="hidden")

# Glutamic acid - Glu

p7 <- Plot_AmAcid(amino_acid = "Glutamic acid", 
                  df = dfCodon_Mean_SD) +
  theme(legend.position="hidden") +
  scale_y_continuous(labels = NULL, breaks = NULL)

# Glycine - Gly 

p8 <- Plot_AmAcid(amino_acid = "Glycine", 
                  df = dfCodon_Mean_SD) +
  theme(legend.position="hidden") +
  scale_y_continuous(labels = NULL, breaks = NULL)

# Histidine - His 

p9 <- Plot_AmAcid(amino_acid = "Histidine", 
                  df = dfCodon_Mean_SD) +
  theme(legend.position="hidden") +
  scale_y_continuous(labels = NULL, breaks = NULL)

# Isoleucine - Ile 

p10 <- Plot_AmAcid(amino_acid = "Isoleucine", 
                  df = dfCodon_Mean_SD) +
  theme(legend.position="hidden") +
  scale_y_continuous(labels = NULL, breaks = NULL)

# Leucine - Leu 

p11 <- Plot_AmAcid(amino_acid = "Leucine", 
                  df = dfCodon_Mean_SD) +
  theme(legend.position="hidden")

# Lysine - Lys 

p12 <- Plot_AmAcid(amino_acid = "Lysine", 
                  df = dfCodon_Mean_SD)  +
  theme(legend.position="hidden") +
  scale_y_continuous(labels = NULL, breaks = NULL)

# Methionine - Met

p13 <- Plot_AmAcid(amino_acid = "Methionine", 
                  df = dfCodon_Mean_SD) +
  theme(legend.position="hidden") +
  scale_y_continuous(labels = NULL, breaks = NULL)

# Phenylalanine - Phe 

p14 <- Plot_AmAcid(amino_acid = "Phenylalanine", 
                  df = dfCodon_Mean_SD) +
  theme(legend.position="hidden") +
  scale_y_continuous(labels = NULL, breaks = NULL)

# Proline - Pro 

p15 <- Plot_AmAcid(amino_acid = "Proline", 
                  df = dfCodon_Mean_SD) +
  theme(legend.position="hidden") +
  scale_y_continuous(labels = NULL, breaks = NULL)

# Serine - Ser

p16 <- Plot_AmAcid(amino_acid = "Serine", 
                  df = dfCodon_Mean_SD) +
  theme(legend.position="hidden")

# Threonine - Thr 

p17 <- Plot_AmAcid(amino_acid = "Threonine", 
                  df = dfCodon_Mean_SD) +
  theme(legend.position="hidden") +
  scale_y_continuous(labels = NULL, breaks = NULL)

# Tryptophan - Trp 

p18 <- Plot_AmAcid(amino_acid = "Tryptophan", 
                  df = dfCodon_Mean_SD) +
  theme(legend.position="hidden") +
  scale_y_continuous(labels = NULL, breaks = NULL)

# Tyrosine - Tyr 

p19 <- Plot_AmAcid(amino_acid = "Tyrosine", 
                  df = dfCodon_Mean_SD) +
  theme(legend.position="hidden") +
  scale_y_continuous(labels = NULL, breaks = NULL)

# Valine - Val 

p20 <- Plot_AmAcid(amino_acid = "Valine", 
                  df = dfCodon_Mean_SD) +
  theme(legend.position="hidden") +
  scale_y_continuous(labels = NULL, breaks = NULL)

#### 13 ARRANGE AMINO ACID PLOTS IN A GRID ####

# Next, using the lemon package, I will be extracting and saving the legend to a separate object.
# The floowing line of code was adapted from https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots (Sturm, 2017).

Legend <- g_legend(p1 + theme(legend.position="bottom"))

# The following line of code was adapted from https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html. (Auguié, 2019).

Combined_Plot <- grid.arrange(p1, p2, p3, p4, p5,
                              p6, p7, p8, p9, p10,
                              p11, p12,p13,p14,p15,
                              p16,p17,p18,p19,p20,
                              nrow=4,
                              ncol=5,
                              top = "Similarities Between Codon Counts Among Nicotiana Codon Usages",
                              bottom = "Codon",
                              left = "Average Count")

# The following line of code was adapted from  https://www.geeksforgeeks.org/add-common-legend-to-combined-ggplot2-plots-in-r/  (Mishra, 2021).

grid.arrange(Combined_Plot, 
             Legend, 
             nrow = 2, 
             heights = c(15, 1))

# The most used codons (one for each amino acid residue) from the three groups appears to be the same.

rm (p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20)

#### 14 GC% CONTENT ####

# In this last section, the overall G-C content of the coding sequences will be calculated.
# First, concatenate all coding sequences to one character string.

# The following line of code was adapted from https://stackoverflow.com/questions/9314328/how-to-collapse-a-list-of-characters-into-a-single-string-in-r (Zacharatos, 2020).

One_Long_Sequence <- toString(dfNew$Trimmed_CDS) 

# For each base, count the number of times that base is encountered and then convert to a percentage.

# The following lines of code (lines 486 to 503) was inspired and adapted from
# https://datacarpentry.org/semester-biology/exercises/Loops-stringr-R/ (White & Brym, n.d.).

# For A

A_Percentage <- (str_count(One_Long_Sequence, "A")) / str_length(One_Long_Sequence) * 100

# For C

C_Percentage <- (str_count(One_Long_Sequence, "C")) / str_length(One_Long_Sequence) * 100

# For G

G_Percentage <- (str_count(One_Long_Sequence, "G")) / str_length(One_Long_Sequence) * 100

# For T

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

# Auguié, B. (2019, July 13). Laying out multiple plots on a page. Comprehensive R Archive Network (CRAN). https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html
# Elek, A. (2019, January 2). Codon usage (CU) analysis in R. Bioconductor.Riken.Jp. https://bioconductor.riken.jp/packages/3.8/bioc/vignettes/coRdon/inst/doc/coRdon.html
# Jaap. (2014, April 8). Remove space between plotted data and the axes. Stack Overflow. https://stackoverflow.com/questions/22945651/remove-space-between-plotted-data-and-the-axes
# jan-glx. (2020, March 12). Rotating and spacing axis labels in ggplot2. Stack Overflow. https://stackoverflow.com/questions/1330989/rotating-and-spacing-axis-labels-in-ggplot2
# LearnByExample. (2022). R Bar Plot – ggplot2. LearnByExample. https://www.learnbyexample.org/r-bar-plot-ggplot2/
# Mishra, P. "mishrapriyank17" (2021, November 28). Add Common Legend to Combined ggplot2 Plots in R. GeeksforGeeks. https://www.geeksforgeeks.org/add-common-legend-to-combined-ggplot2-plots-in-r/
# shadow. (2014, July 13). Vector of most used codons from table of codon usage in R. Stack Overflow. https://stackoverflow.com/questions/24546312/vector-of-most-used-codons-from-table-of-codon-usage-in-r
# Sturm, G. (2017, June 26). Add a common Legend for combined ggplots. Stack Overflow. https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
# White, E., & Brym, Z. (n.d.). stringr (Loops). Data Carpentry for Biologists. https://datacarpentry.org/semester-biology/exercises/Loops-stringr-R/
# Zacharatos, D. (2020, August 4). How to collapse a list of characters into a single string in R. Stack Overflow. https://stackoverflow.com/questions/9314328/how-to-collapse-a-list-of-characters-into-a-single-string-in-r
