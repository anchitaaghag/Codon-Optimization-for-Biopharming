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

# Statistical Tests
# 27 June 2022
# Anchitaa Ghag

#### IMPORT DATA ####
dfCodingSeqs <- read_csv("/Users/anchitaa/Major_Research_Project_2022/06_Code/Codon-Optimization-for-Biopharming/dfCodingSeqs.csv")[,2:14]
dfExisting <- read_csv("/Users/anchitaa/Major_Research_Project_2022/06_Code/Codon-Optimization-for-Biopharming/dfKazuza.csv")[,2:5]
dfN.Tabacum <- read_csv("/Users/anchitaa/Major_Research_Project_2022/06_Code/Codon-Optimization-for-Biopharming/dfNTabacum.csv")[,2:5]

# List of the codon usage from each CDS that makes up Kazuza CUT is available at: http://www.kazusa.or.jp/codon/current/species/4100
KazuzaCDS <- read.table("/Users/anchitaa/Major_Research_Project_2022/06_Code/KCC_No_Carot.txt",
                        header = TRUE)

# Finally, calculate the most used codons (1 for each amino acid residue) from the data frame.

#dfKazuza.max <- group_by(dfKazuza, AmAcid) %>% 
#  filter(Number==max(Number)) %>% 
#  select(AmAcid, Codon)


#### GENERATE CODON USAGE TABLES ####

#https://bioconductor.riken.jp/packages/3.8/bioc/vignettes/coRdon/inst/doc/coRdon.html

# 1) First, convert the trimmed coding sequences to a DNAStringSet object.

cds <- DNAStringSet(dfCodingSeqs$Trimmed_CDS)

# 2) Second, create a codon usage table object for both existing and current codon usage tables.
set.seed(1212)

Existing_CU <- codonTable(KazuzaCDS)
Current_CU <- codonTable(cds)

# 3) Then, 

# FIXME Uisng codRon package functions. Add Commenting.

CU <- codonCounts(Current_CU)
#getlen(Current_CU)
#row.names(CU) <- dfCodingSeqs$Titles
#(colMeans(CU)/colSums(CU))*1000

# Finally, create an consensus codon usage table based on the average empirical codon counts of coding sequences.

Avg_Codon_Counts <- colMeans(CU)

weighted.mean(CU$AAA)

Freq_Per_1000 <- list()
for (i in Avg_Codon_Counts){
Value_Calc <- (i/sum(Avg_Codon_Counts))*1000
Freq_Per_1000 <- append(Freq_Per_1000,Value_Calc )
}

dfTemp <- data.frame(colnames(CU),Avg_Codon_Counts,unlist(Freq_Per_1000))

# Reorder the columns to ensure they match.
# https://stackoverflow.com/questions/11977102/order-data-frame-rows-according-to-vector-with-specific-order

dfTemp.sub <- dfTemp[match(dfExisting$Codon, dfTemp$colnames.CU.),]

# Add to final consensus data frame.

dfCurrent <- data.frame(dfExisting$AmAcid,dfExisting$Codon,dfTemp.sub$Avg_Codon_Counts,dfTemp.sub$unlist.Freq_Per_1000.)
colnames(dfCurrent) <- colnames(dfExisting)

rm(cds,dfTemp,dfTemp.sub)

#### CHI SQUARE TEST: ####

# Perform a chi-squared test to evaluate if there is a difference between the frequencies of both CU.

chisq.test(x = dfExisting$X.1000,
           y = dfCurrent$X.1000)

# Pearson's Chi-squared test

# data:  dfExisting$X.1000 and dfCurrent$X.1000
# X-squared = 3904, df = 3843, p-value = 0.242

#######
# https://www.bioconductor.org/packages/devel/bioc/vignettes/coRdon/inst/doc/coRdon.html#calculate-cu-bias

milc <- MILC(Current_CU)
head(milc)








# "Other CU statistics can be calculated in the same way as MILC(), using one of the functions: B(), MCB(), ENCprime(), ENC() or SCUO(). Note however, that when calculating ENC and SCUO, one doesn’t need to provide a subset of refe"

SCUO(Current_CU)
ENC(Current_CU)
B(Current_CU)
MCB(Current_CU)
ENCprime(Current_CU)

# Next, compare the CU bias for every coding sequence between the created CUT and Kazuza through visualizations.
# 

library(ggplot2)
xlab <- "MILC distance from sample centroid"
ylab <- "MILC distance from ribosomal genes"

MILC_Created_CUT <- MILC(Current_CU, ribosomal = TRUE)
Bplot(x = "self", y = "self", data = MILC_Created_CUT) +
  labs(x = xlab, y = ylab)

# Visualization of Codon Usage - Existing vs. Current Karlin B plot ####

# Code Adapted From: https://www.bioconductor.org/packages/devel/bioc/vignettes/coRdon/inst/doc/coRdon.html#calculate-cu-bias

# Use the intraBplot() function in the coRdon package and our two codonTable objects to plot a  plot of codon usage distances between existing vs. created codon tables.
#Intra-samples Karlin B plot

set.seed(1111)
Sample_CU <- sample(Current_CU,100)


intraBplot(x = Existing_CU, 
           y = Sample_CU, 
           names = c("Existing", "Consensus"), 
           variable = "MILC", 
           size = 3, 
           alpha = 1.0) + 
  ggtitle("Karlin B Plot of Existing (Kazuza) v.s. Consensus CU Distances")+
  xlab("MILC Distance from Existing CU") +
  ylab("MILC Distance from Consensus CU") + 
  geom_jitter(width  = 0.05) + # Add noise since the dataset in small to avoid overplotting.
  #geom_smooth(method = "lm", formula = y ~ x, fullrange = FALSE, level = 0.95) # 95% confidence interval
  geom_smooth(aes(group = 1), formula = y ~ poly(x,2), method = "lm", fullrange = TRUE, level = 0.95) 
  #geom_curve()

# The example dataset in the vignette has ~ 19, 000 + "points". There may just be undersampling/ less data points.




#hi <- lapply(list, function)

formatted_cds <- s2c((unlist(dfCodingSeqs$Trimmed_CDS[1])))

Current_RSCU <- uco(seq=formatted_cds,
    frame = 0,
    index = "rscu",
    #as.data.frame = TRUE,
    NA.rscu = NA)


table(Current_RSCU == 1) # RSCU = 1 codon is used as expected by random usage 
table(Current_RSCU > 1) # RSCU > 1 codon used more frequently than random
table(Current_RSCU < 1) # RSCU < 1 codon used less frequently than random 

Current_RSCU[Current_RSCU < 1]

hist(Current_RSCU[Current_RSCU < 1])

dfCurrent

enc <- ENC(Current_CU)

table(enc == 20) # RSCU = 1 codon is used as expected by random usage 
table(enc == 61) # RSCU > 1 codon used more frequently than random
table(enc < 40) # RSCU < 1 codon used less frequently than random 

enc[enc < 40]

?uco
rcds <- read.fasta(file = system.file("sequences/malM.fasta", package = "seqinr"))[[1]]
uco( rcds, index = "freq")
uco( rcds, index = "eff")
uco( rcds, index = "rscu")
uco( rcds, as.data.frame = TRUE)


#### The Twenty Amino Acid Distributions ####

library(stringr)

# Alanine - Ala

# Takes an argument "amino_acid" which is the name of the amino acid. Three data frames in the style of the Kazuza data base.

Plot_AmAcid <- function(amino_acid, df1, df2, df3) {

# Shortens to three letter format of the amino acid name.
  
Three_Letter_Form <- substr(amino_acid, start = 1, stop = 3)

# Counts the number of times to repeat the character string.

Times_To_Repeat <- as.numeric(table(df1$AmAcid == Three_Letter_Form)["TRUE"])

# Creates a data frame containing the codon count frequencies from each data frame.
dfAmAcid <- rbind(
  (subset(df1, AmAcid == Three_Letter_Form))[,c(2,4)],
  (subset(df2, AmAcid == Three_Letter_Form))[,c(2,4)],
  (subset(df3, AmAcid == Three_Letter_Form))[,c(2,4)]
)


Sample <- c(rep(x = "Existing", times = Times_To_Repeat),
            rep(x = "Current", times = Times_To_Repeat),
            rep(x = "N. tabacum", times = Times_To_Repeat))

str1 <- (subset(df1, AmAcid == Three_Letter_Form))[,4]
s1 <-  unlist(str1)

str2 <- (subset(df2, AmAcid == Three_Letter_Form))[,4]
s2 <-  unlist(str2)

str3 <- (subset(df3, AmAcid == Three_Letter_Form))[,4]
s3 <-  unlist(str3)

Std_Dev <- c(rep(x = sd(s1), times = Times_To_Repeat),
             rep(x = sd(s2), times = Times_To_Repeat),
             rep(x = sd(s3), times = Times_To_Repeat))


dfAmAcid["Sample"] <- Sample
dfAmAcid["Std_Dev"] <- Std_Dev

# Generate the plot using ggplot package functions.

Plot <- ggplot(data=dfAmAcid, 
       aes(x=Codon, y=X.1000, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle(amino_acid) +
  theme_minimal() +
  xlab("Codon") +
  ylab("Frequency (Per 1000)") +
  scale_fill_discrete(name = "CU Table", labels = c("Existing (Kazuza)", "Current", "N. tabacum"))+
  geom_errorbar(aes(ymin=X.1000-Std_Dev, ymax=X.1000+Std_Dev), width=.2, position=position_dodge(.9))

return(Plot)

}

# Arginine - Arg 

a.a. <- "Arg"
No_Time_To_Repeat <- as.numeric(table(dfExisting$AmAcid == a.a.)["TRUE"])

dfAmAcid <- rbind(
  (subset(dfExisting, AmAcid == a.a.))[,c(2,4)],
  (subset(dfCurrent, AmAcid == a.a.))[,c(2,4)],
  (subset(dfN.Tabacum, AmAcid == a.a.))[,c(2,4)]
)

Sample <- c(rep(x = "Existing", times = No_Time_To_Repeat),
            rep(x = "Current", times = No_Time_To_Repeat),
            rep(x = "N. tabacum", times = No_Time_To_Repeat))

str1 <- (subset(dfExisting, AmAcid == a.a.))[,4]
s1 <-  unlist(str1)

str2 <- (subset(dfCurrent, AmAcid == a.a.))[,4]
s2 <-  unlist(str2)

str3 <- (subset(dfN.Tabacum, AmAcid == a.a.))[,4]
s3 <-  unlist(str3)

Std_Dev <- c(rep(x = sd(s1)/2, times = No_Time_To_Repeat),
             rep(x = sd(s2)/2, times = No_Time_To_Repeat),
             rep(x = sd(s3)/2, times = No_Time_To_Repeat))


dfAmAcid["Sample"] <- Sample
dfAmAcid["Std_Dev"] <- Std_Dev


p2 <- ggplot(data=dfAmAcid, 
       aes(x=Codon, y=X.1000, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle(a.a.) +
  theme_minimal() +
  xlab("Codon") +
  ylab("Frequency (Per 1000)") +
  scale_fill_discrete(name = "CU Table", labels = c("Existing (Kazuza)", "Current", "N. tabacum"))+
  geom_errorbar(aes(ymin=X.1000-Std_Dev, ymax=X.1000+Std_Dev), width=.2,
                position=position_dodge(.9))

# Asparagine - Asn

a.a. <- "Asn"
No_Time_To_Repeat <- as.numeric(table(dfExisting$AmAcid == a.a.)["TRUE"])

dfAmAcid <- rbind(
  (subset(dfExisting, AmAcid == a.a.))[,c(2,4)],
  (subset(dfCurrent, AmAcid == a.a.))[,c(2,4)],
  (subset(dfN.Tabacum, AmAcid == a.a.))[,c(2,4)]
)

Sample <- c(rep(x = "Existing", times = No_Time_To_Repeat),
            rep(x = "Current", times = No_Time_To_Repeat),
            rep(x = "N. tabacum", times = No_Time_To_Repeat))

str1 <- (subset(dfExisting, AmAcid == a.a.))[,4]
s1 <-  unlist(str1)

str2 <- (subset(dfCurrent, AmAcid == a.a.))[,4]
s2 <-  unlist(str2)

str3 <- (subset(dfN.Tabacum, AmAcid == a.a.))[,4]
s3 <-  unlist(str3)

Std_Dev <- c(rep(x = sd(s1)/2, times = No_Time_To_Repeat),
             rep(x = sd(s2)/2, times = No_Time_To_Repeat),
             rep(x = sd(s3)/2, times = No_Time_To_Repeat))


dfAmAcid["Sample"] <- Sample
dfAmAcid["Std_Dev"] <- Std_Dev


p3 <- ggplot(data=dfAmAcid, 
       aes(x=Codon, y=X.1000, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle(a.a.) +
  theme_minimal() +
  xlab("Codon") +
  ylab("Frequency (Per 1000)") +
  scale_fill_discrete(name = "CU Table", labels = c("Existing (Kazuza)", "Current", "N. tabacum"))+
  geom_errorbar(aes(ymin=X.1000-Std_Dev, ymax=X.1000+Std_Dev), width=.2,
                position=position_dodge(.9))

# Aspartic acid - Asp 

a.a. <- "Asp"
No_Time_To_Repeat <- as.numeric(table(dfExisting$AmAcid == a.a.)["TRUE"])

dfAmAcid <- rbind(
  (subset(dfExisting, AmAcid == a.a.))[,c(2,4)],
  (subset(dfCurrent, AmAcid == a.a.))[,c(2,4)],
  (subset(dfN.Tabacum, AmAcid == a.a.))[,c(2,4)]
)

Sample <- c(rep(x = "Existing", times = No_Time_To_Repeat),
            rep(x = "Current", times = No_Time_To_Repeat),
            rep(x = "N. tabacum", times = No_Time_To_Repeat))

str1 <- (subset(dfExisting, AmAcid == a.a.))[,4]
s1 <-  unlist(str1)

str2 <- (subset(dfCurrent, AmAcid == a.a.))[,4]
s2 <-  unlist(str2)

str3 <- (subset(dfN.Tabacum, AmAcid == a.a.))[,4]
s3 <-  unlist(str3)

Std_Dev <- c(rep(x = sd(s1)/2, times = No_Time_To_Repeat),
             rep(x = sd(s2)/2, times = No_Time_To_Repeat),
             rep(x = sd(s3)/2, times = No_Time_To_Repeat))


dfAmAcid["Sample"] <- Sample
dfAmAcid["Std_Dev"] <- Std_Dev


p4 <- ggplot(data=dfAmAcid, 
       aes(x=Codon, y=X.1000, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle(a.a.) +
  theme_minimal() +
  xlab("Codon") +
  ylab("Frequency (Per 1000)") +
  scale_fill_discrete(name = "CU Table", labels = c("Existing (Kazuza)", "Current", "N. tabacum"))+
  geom_errorbar(aes(ymin=X.1000-Std_Dev, ymax=X.1000+Std_Dev), width=.2,
                position=position_dodge(.9))

# Cysteine - Cys

a.a. <- "Cys"
No_Time_To_Repeat <- as.numeric(table(dfExisting$AmAcid == a.a.)["TRUE"])

dfAmAcid <- rbind(
  (subset(dfExisting, AmAcid == a.a.))[,c(2,4)],
  (subset(dfCurrent, AmAcid == a.a.))[,c(2,4)],
  (subset(dfN.Tabacum, AmAcid == a.a.))[,c(2,4)]
)

Sample <- c(rep(x = "Existing", times = No_Time_To_Repeat),
            rep(x = "Current", times = No_Time_To_Repeat),
            rep(x = "N. tabacum", times = No_Time_To_Repeat))

str1 <- (subset(dfExisting, AmAcid == a.a.))[,4]
s1 <-  unlist(str1)

str2 <- (subset(dfCurrent, AmAcid == a.a.))[,4]
s2 <-  unlist(str2)

str3 <- (subset(dfN.Tabacum, AmAcid == a.a.))[,4]
s3 <-  unlist(str3)

Std_Dev <- c(rep(x = sd(s1)/2, times = No_Time_To_Repeat),
             rep(x = sd(s2)/2, times = No_Time_To_Repeat),
             rep(x = sd(s3)/2, times = No_Time_To_Repeat))


dfAmAcid["Sample"] <- Sample
dfAmAcid["Std_Dev"] <- Std_Dev


p5 <- ggplot(data=dfAmAcid, 
       aes(x=Codon, y=X.1000, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle(a.a.) +
  theme_minimal() +
  xlab("Codon") +
  ylab("Frequency (Per 1000)") +
  scale_fill_discrete(name = "CU Table", labels = c("Existing (Kazuza)", "Current", "N. tabacum"))+
  geom_errorbar(aes(ymin=X.1000-Std_Dev, ymax=X.1000+Std_Dev), width=.2,
                position=position_dodge(.9))

# Glutamine - Gln

a.a. <- "Gln"
No_Time_To_Repeat <- as.numeric(table(dfExisting$AmAcid == a.a.)["TRUE"])

dfAmAcid <- rbind(
  (subset(dfExisting, AmAcid == a.a.))[,c(2,4)],
  (subset(dfCurrent, AmAcid == a.a.))[,c(2,4)],
  (subset(dfN.Tabacum, AmAcid == a.a.))[,c(2,4)]
)

Sample <- c(rep(x = "Existing", times = No_Time_To_Repeat),
            rep(x = "Current", times = No_Time_To_Repeat),
            rep(x = "N. tabacum", times = No_Time_To_Repeat))

str1 <- (subset(dfExisting, AmAcid == a.a.))[,4]
s1 <-  unlist(str1)

str2 <- (subset(dfCurrent, AmAcid == a.a.))[,4]
s2 <-  unlist(str2)

str3 <- (subset(dfN.Tabacum, AmAcid == a.a.))[,4]
s3 <-  unlist(str3)

Std_Dev <- c(rep(x = sd(s1)/2, times = No_Time_To_Repeat),
             rep(x = sd(s2)/2, times = No_Time_To_Repeat),
             rep(x = sd(s3)/2, times = No_Time_To_Repeat))


dfAmAcid["Sample"] <- Sample
dfAmAcid["Std_Dev"] <- Std_Dev


p6 <- ggplot(data=dfAmAcid, 
       aes(x=Codon, y=X.1000, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle(a.a.) +
  theme_minimal() +
  xlab("Codon") +
  ylab("Frequency (Per 1000)") +
  scale_fill_discrete(name = "CU Table", labels = c("Existing (Kazuza)", "Current", "N. tabacum"))+
  geom_errorbar(aes(ymin=X.1000-Std_Dev, ymax=X.1000+Std_Dev), width=.2,
                position=position_dodge(.9))

# Glutamic acid - Glu

a.a. <- "Glu"
No_Time_To_Repeat <- as.numeric(table(dfExisting$AmAcid == a.a.)["TRUE"])

dfAmAcid <- rbind(
  (subset(dfExisting, AmAcid == a.a.))[,c(2,4)],
  (subset(dfCurrent, AmAcid == a.a.))[,c(2,4)],
  (subset(dfN.Tabacum, AmAcid == a.a.))[,c(2,4)]
)

Sample <- c(rep(x = "Existing", times = No_Time_To_Repeat),
            rep(x = "Current", times = No_Time_To_Repeat),
            rep(x = "N. tabacum", times = No_Time_To_Repeat))

str1 <- (subset(dfExisting, AmAcid == a.a.))[,4]
s1 <-  unlist(str1)

str2 <- (subset(dfCurrent, AmAcid == a.a.))[,4]
s2 <-  unlist(str2)

str3 <- (subset(dfN.Tabacum, AmAcid == a.a.))[,4]
s3 <-  unlist(str3)

Std_Dev <- c(rep(x = sd(s1)/2, times = No_Time_To_Repeat),
             rep(x = sd(s2)/2, times = No_Time_To_Repeat),
             rep(x = sd(s3)/2, times = No_Time_To_Repeat))


dfAmAcid["Sample"] <- Sample
dfAmAcid["Std_Dev"] <- Std_Dev


p7 <- ggplot(data=dfAmAcid, 
       aes(x=Codon, y=X.1000, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle(a.a.) +
  theme_minimal() +
  xlab("Codon") +
  ylab("Frequency (Per 1000)") +
  scale_fill_discrete(name = "CU Table", labels = c("Existing (Kazuza)", "Current", "N. tabacum"))+
  geom_errorbar(aes(ymin=X.1000-Std_Dev, ymax=X.1000+Std_Dev), width=.2,
                position=position_dodge(.9))

# Glycine - Gly 

a.a. <- "Gly"
No_Time_To_Repeat <- as.numeric(table(dfExisting$AmAcid == a.a.)["TRUE"])

dfAmAcid <- rbind(
  (subset(dfExisting, AmAcid == a.a.))[,c(2,4)],
  (subset(dfCurrent, AmAcid == a.a.))[,c(2,4)],
  (subset(dfN.Tabacum, AmAcid == a.a.))[,c(2,4)]
)

Sample <- c(rep(x = "Existing", times = No_Time_To_Repeat),
            rep(x = "Current", times = No_Time_To_Repeat),
            rep(x = "N. tabacum", times = No_Time_To_Repeat))

str1 <- (subset(dfExisting, AmAcid == a.a.))[,4]
s1 <-  unlist(str1)

str2 <- (subset(dfCurrent, AmAcid == a.a.))[,4]
s2 <-  unlist(str2)

str3 <- (subset(dfN.Tabacum, AmAcid == a.a.))[,4]
s3 <-  unlist(str3)

Std_Dev <- c(rep(x = sd(s1)/2, times = No_Time_To_Repeat),
             rep(x = sd(s2)/2, times = No_Time_To_Repeat),
             rep(x = sd(s3)/2, times = No_Time_To_Repeat))


dfAmAcid["Sample"] <- Sample
dfAmAcid["Std_Dev"] <- Std_Dev


p8 <- ggplot(data=dfAmAcid, 
       aes(x=Codon, y=X.1000, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle(a.a.) +
  theme_minimal() +
  xlab("Codon") +
  ylab("Frequency (Per 1000)") +
  scale_fill_discrete(name = "CU Table", labels = c("Existing (Kazuza)", "Current", "N. tabacum"))+
  geom_errorbar(aes(ymin=X.1000-Std_Dev, ymax=X.1000+Std_Dev), width=.2,
                position=position_dodge(.9))

# Histidine - His 

a.a. <- "His"
No_Time_To_Repeat <- as.numeric(table(dfExisting$AmAcid == a.a.)["TRUE"])

dfAmAcid <- rbind(
  (subset(dfExisting, AmAcid == a.a.))[,c(2,4)],
  (subset(dfCurrent, AmAcid == a.a.))[,c(2,4)],
  (subset(dfN.Tabacum, AmAcid == a.a.))[,c(2,4)]
)

Sample <- c(rep(x = "Existing", times = No_Time_To_Repeat),
            rep(x = "Current", times = No_Time_To_Repeat),
            rep(x = "N. tabacum", times = No_Time_To_Repeat))

str1 <- (subset(dfExisting, AmAcid == a.a.))[,4]
s1 <-  unlist(str1)

str2 <- (subset(dfCurrent, AmAcid == a.a.))[,4]
s2 <-  unlist(str2)

str3 <- (subset(dfN.Tabacum, AmAcid == a.a.))[,4]
s3 <-  unlist(str3)

Std_Dev <- c(rep(x = sd(s1)/2, times = No_Time_To_Repeat),
             rep(x = sd(s2)/2, times = No_Time_To_Repeat),
             rep(x = sd(s3)/2, times = No_Time_To_Repeat))


dfAmAcid["Sample"] <- Sample
dfAmAcid["Std_Dev"] <- Std_Dev


p9<- ggplot(data=dfAmAcid, 
       aes(x=Codon, y=X.1000, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle(a.a.) +
  theme_minimal() +
  xlab("Codon") +
  ylab("Frequency (Per 1000)") +
  scale_fill_discrete(name = "CU Table", labels = c("Existing (Kazuza)", "Current", "N. tabacum"))+
  geom_errorbar(aes(ymin=X.1000-Std_Dev, ymax=X.1000+Std_Dev), width=.2,
                position=position_dodge(.9))

# Isoleucine - Ile 

a.a. <- "Ile"
No_Time_To_Repeat <- as.numeric(table(dfExisting$AmAcid == a.a.)["TRUE"])

dfAmAcid <- rbind(
  (subset(dfExisting, AmAcid == a.a.))[,c(2,4)],
  (subset(dfCurrent, AmAcid == a.a.))[,c(2,4)],
  (subset(dfN.Tabacum, AmAcid == a.a.))[,c(2,4)]
)

Sample <- c(rep(x = "Existing", times = No_Time_To_Repeat),
            rep(x = "Current", times = No_Time_To_Repeat),
            rep(x = "N. tabacum", times = No_Time_To_Repeat))

str1 <- (subset(dfExisting, AmAcid == a.a.))[,4]
s1 <-  unlist(str1)

str2 <- (subset(dfCurrent, AmAcid == a.a.))[,4]
s2 <-  unlist(str2)

str3 <- (subset(dfN.Tabacum, AmAcid == a.a.))[,4]
s3 <-  unlist(str3)

Std_Dev <- c(rep(x = sd(s1)/2, times = No_Time_To_Repeat),
             rep(x = sd(s2)/2, times = No_Time_To_Repeat),
             rep(x = sd(s3)/2, times = No_Time_To_Repeat))


dfAmAcid["Sample"] <- Sample
dfAmAcid["Std_Dev"] <- Std_Dev


p10 <- ggplot(data=dfAmAcid, 
       aes(x=Codon, y=X.1000, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle(a.a.) +
  theme_minimal() +
  xlab("Codon") +
  ylab("Frequency (Per 1000)") +
  scale_fill_discrete(name = "CU Table", labels = c("Existing (Kazuza)", "Current", "N. tabacum"))+
  geom_errorbar(aes(ymin=X.1000-Std_Dev, ymax=X.1000+Std_Dev), width=.2,
                position=position_dodge(.9))

# Leucine - Leu 

a.a. <- "Leu"
No_Time_To_Repeat <- as.numeric(table(dfExisting$AmAcid == a.a.)["TRUE"])

dfAmAcid <- rbind(
  (subset(dfExisting, AmAcid == a.a.))[,c(2,4)],
  (subset(dfCurrent, AmAcid == a.a.))[,c(2,4)],
  (subset(dfN.Tabacum, AmAcid == a.a.))[,c(2,4)]
)

Sample <- c(rep(x = "Existing", times = No_Time_To_Repeat),
            rep(x = "Current", times = No_Time_To_Repeat),
            rep(x = "N. tabacum", times = No_Time_To_Repeat))

str1 <- (subset(dfExisting, AmAcid == a.a.))[,4]
s1 <-  unlist(str1)

str2 <- (subset(dfCurrent, AmAcid == a.a.))[,4]
s2 <-  unlist(str2)

str3 <- (subset(dfN.Tabacum, AmAcid == a.a.))[,4]
s3 <-  unlist(str3)

Std_Dev <- c(rep(x = sd(s1)/2, times = No_Time_To_Repeat),
             rep(x = sd(s2)/2, times = No_Time_To_Repeat),
             rep(x = sd(s3)/2, times = No_Time_To_Repeat))


dfAmAcid["Sample"] <- Sample
dfAmAcid["Std_Dev"] <- Std_Dev


p11<- ggplot(data=dfAmAcid, 
       aes(x=Codon, y=X.1000, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle(a.a.) +
  theme_minimal() +
  xlab("Codon") +
  ylab("Frequency (Per 1000)") +
  scale_fill_discrete(name = "CU Table", labels = c("Existing (Kazuza)", "Current", "N. tabacum"))+
  geom_errorbar(aes(ymin=X.1000-Std_Dev, ymax=X.1000+Std_Dev), width=.2,
                position=position_dodge(.9))

# Lysine - Lys 

a.a. <- "Lys"
No_Time_To_Repeat <- as.numeric(table(dfExisting$AmAcid == a.a.)["TRUE"])

dfAmAcid <- rbind(
  (subset(dfExisting, AmAcid == a.a.))[,c(2,4)],
  (subset(dfCurrent, AmAcid == a.a.))[,c(2,4)],
  (subset(dfN.Tabacum, AmAcid == a.a.))[,c(2,4)]
)

Sample <- c(rep(x = "Existing", times = No_Time_To_Repeat),
            rep(x = "Current", times = No_Time_To_Repeat),
            rep(x = "N. tabacum", times = No_Time_To_Repeat))

str1 <- (subset(dfExisting, AmAcid == a.a.))[,4]
s1 <-  unlist(str1)

str2 <- (subset(dfCurrent, AmAcid == a.a.))[,4]
s2 <-  unlist(str2)

str3 <- (subset(dfN.Tabacum, AmAcid == a.a.))[,4]
s3 <-  unlist(str3)

Std_Dev <- c(rep(x = sd(s1)/2, times = No_Time_To_Repeat),
             rep(x = sd(s2)/2, times = No_Time_To_Repeat),
             rep(x = sd(s3)/2, times = No_Time_To_Repeat))


dfAmAcid["Sample"] <- Sample
dfAmAcid["Std_Dev"] <- Std_Dev


p12 <- ggplot(data=dfAmAcid, 
       aes(x=Codon, y=X.1000, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle(a.a.) +
  theme_minimal() +
  xlab("Codon") +
  ylab("Frequency (Per 1000)") +
  scale_fill_discrete(name = "CU Table", labels = c("Existing (Kazuza)", "Current", "N. tabacum"))+
  geom_errorbar(aes(ymin=X.1000-Std_Dev, ymax=X.1000+Std_Dev), width=.2,
                position=position_dodge(.9))

# Methionine - Met

a.a. <- "Met"
No_Time_To_Repeat <- as.numeric(table(dfExisting$AmAcid == a.a.)["TRUE"])

dfAmAcid <- rbind(
  (subset(dfExisting, AmAcid == a.a.))[,c(2,4)],
  (subset(dfCurrent, AmAcid == a.a.))[,c(2,4)],
  (subset(dfN.Tabacum, AmAcid == a.a.))[,c(2,4)]
)

Sample <- c(rep(x = "Existing", times = No_Time_To_Repeat),
            rep(x = "Current", times = No_Time_To_Repeat),
            rep(x = "N. tabacum", times = No_Time_To_Repeat))

str1 <- (subset(dfExisting, AmAcid == a.a.))[,4]
s1 <-  unlist(str1)

str2 <- (subset(dfCurrent, AmAcid == a.a.))[,4]
s2 <-  unlist(str2)

str3 <- (subset(dfN.Tabacum, AmAcid == a.a.))[,4]
s3 <-  unlist(str3)

Std_Dev <- c(rep(x = sd(s1)/2, times = No_Time_To_Repeat),
             rep(x = sd(s2)/2, times = No_Time_To_Repeat),
             rep(x = sd(s3)/2, times = No_Time_To_Repeat))


dfAmAcid["Sample"] <- Sample
dfAmAcid["Std_Dev"] <- Std_Dev


p13<- ggplot(data=dfAmAcid, 
       aes(x=Codon, y=X.1000, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle(a.a.) +
  theme_minimal() +
  xlab("Codon") +
  ylab("Frequency (Per 1000)") +
  scale_fill_discrete(name = "CU Table", labels = c("Existing (Kazuza)", "Current", "N. tabacum"))+
  geom_errorbar(aes(ymin=X.1000-Std_Dev, ymax=X.1000+Std_Dev), width=.2,
                position=position_dodge(.9))

# Phenylalanine - Phe 

a.a. <- "Phe"
No_Time_To_Repeat <- as.numeric(table(dfExisting$AmAcid == a.a.)["TRUE"])

dfAmAcid <- rbind(
  (subset(dfExisting, AmAcid == a.a.))[,c(2,4)],
  (subset(dfCurrent, AmAcid == a.a.))[,c(2,4)],
  (subset(dfN.Tabacum, AmAcid == a.a.))[,c(2,4)]
)

Sample <- c(rep(x = "Existing", times = No_Time_To_Repeat),
            rep(x = "Current", times = No_Time_To_Repeat),
            rep(x = "N. tabacum", times = No_Time_To_Repeat))

str1 <- (subset(dfExisting, AmAcid == a.a.))[,4]
s1 <-  unlist(str1)

str2 <- (subset(dfCurrent, AmAcid == a.a.))[,4]
s2 <-  unlist(str2)

str3 <- (subset(dfN.Tabacum, AmAcid == a.a.))[,4]
s3 <-  unlist(str3)

Std_Dev <- c(rep(x = sd(s1)/2, times = No_Time_To_Repeat),
             rep(x = sd(s2)/2, times = No_Time_To_Repeat),
             rep(x = sd(s3)/2, times = No_Time_To_Repeat))


dfAmAcid["Sample"] <- Sample
dfAmAcid["Std_Dev"] <- Std_Dev


p14 <- ggplot(data=dfAmAcid, 
       aes(x=Codon, y=X.1000, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle(a.a.) +
  theme_minimal() +
  xlab("Codon") +
  ylab("Frequency (Per 1000)") +
  scale_fill_discrete(name = "CU Table", labels = c("Existing (Kazuza)", "Current", "N. tabacum"))+
  geom_errorbar(aes(ymin=X.1000-Std_Dev, ymax=X.1000+Std_Dev), width=.2,
                position=position_dodge(.9))

# Proline - Pro 

a.a. <- "Pro"
No_Time_To_Repeat <- as.numeric(table(dfExisting$AmAcid == a.a.)["TRUE"])

dfAmAcid <- rbind(
  (subset(dfExisting, AmAcid == a.a.))[,c(2,4)],
  (subset(dfCurrent, AmAcid == a.a.))[,c(2,4)],
  (subset(dfN.Tabacum, AmAcid == a.a.))[,c(2,4)]
)

Sample <- c(rep(x = "Existing", times = No_Time_To_Repeat),
            rep(x = "Current", times = No_Time_To_Repeat),
            rep(x = "N. tabacum", times = No_Time_To_Repeat))

str1 <- (subset(dfExisting, AmAcid == a.a.))[,4]
s1 <-  unlist(str1)

str2 <- (subset(dfCurrent, AmAcid == a.a.))[,4]
s2 <-  unlist(str2)

str3 <- (subset(dfN.Tabacum, AmAcid == a.a.))[,4]
s3 <-  unlist(str3)

Std_Dev <- c(rep(x = sd(s1)/2, times = No_Time_To_Repeat),
             rep(x = sd(s2)/2, times = No_Time_To_Repeat),
             rep(x = sd(s3)/2, times = No_Time_To_Repeat))


dfAmAcid["Sample"] <- Sample
dfAmAcid["Std_Dev"] <- Std_Dev


p15 <- ggplot(data=dfAmAcid, 
       aes(x=Codon, y=X.1000, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle(a.a.) +
  theme_minimal() +
  xlab("Codon") +
  ylab("Frequency (Per 1000)") +
  scale_fill_discrete(name = "CU Table", labels = c("Existing (Kazuza)", "Current", "N. tabacum"))+
  geom_errorbar(aes(ymin=X.1000-Std_Dev, ymax=X.1000+Std_Dev), width=.2,
                position=position_dodge(.9))

# Serine - Ser

a.a. <- "Ser"
No_Time_To_Repeat <- as.numeric(table(dfExisting$AmAcid == a.a.)["TRUE"])

dfAmAcid <- rbind(
  (subset(dfExisting, AmAcid == a.a.))[,c(2,4)],
  (subset(dfCurrent, AmAcid == a.a.))[,c(2,4)],
  (subset(dfN.Tabacum, AmAcid == a.a.))[,c(2,4)]
)

Sample <- c(rep(x = "Existing", times = No_Time_To_Repeat),
            rep(x = "Current", times = No_Time_To_Repeat),
            rep(x = "N. tabacum", times = No_Time_To_Repeat))

str1 <- (subset(dfExisting, AmAcid == a.a.))[,4]
s1 <-  unlist(str1)

str2 <- (subset(dfCurrent, AmAcid == a.a.))[,4]
s2 <-  unlist(str2)

str3 <- (subset(dfN.Tabacum, AmAcid == a.a.))[,4]
s3 <-  unlist(str3)

Std_Dev <- c(rep(x = sd(s1)/2, times = No_Time_To_Repeat),
             rep(x = sd(s2)/2, times = No_Time_To_Repeat),
             rep(x = sd(s3)/2, times = No_Time_To_Repeat))


dfAmAcid["Sample"] <- Sample
dfAmAcid["Std_Dev"] <- Std_Dev


p16 <- ggplot(data=dfAmAcid, 
       aes(x=Codon, y=X.1000, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle(a.a.) +
  theme_minimal() +
  xlab("Codon") +
  ylab("Frequency (Per 1000)") +
  scale_fill_discrete(name = "CU Table", labels = c("Existing (Kazuza)", "Current", "N. tabacum"))+
  geom_errorbar(aes(ymin=X.1000-Std_Dev, ymax=X.1000+Std_Dev), width=.2,
                position=position_dodge(.9))

# Threonine - Thr 

a.a. <- "Thr"
No_Time_To_Repeat <- as.numeric(table(dfExisting$AmAcid == a.a.)["TRUE"])

dfAmAcid <- rbind(
  (subset(dfExisting, AmAcid == a.a.))[,c(2,4)],
  (subset(dfCurrent, AmAcid == a.a.))[,c(2,4)],
  (subset(dfN.Tabacum, AmAcid == a.a.))[,c(2,4)]
)

Sample <- c(rep(x = "Existing", times = No_Time_To_Repeat),
            rep(x = "Current", times = No_Time_To_Repeat),
            rep(x = "N. tabacum", times = No_Time_To_Repeat))

str1 <- (subset(dfExisting, AmAcid == a.a.))[,4]
s1 <-  unlist(str1)

str2 <- (subset(dfCurrent, AmAcid == a.a.))[,4]
s2 <-  unlist(str2)

str3 <- (subset(dfN.Tabacum, AmAcid == a.a.))[,4]
s3 <-  unlist(str3)

Std_Dev <- c(rep(x = sd(s1)/2, times = No_Time_To_Repeat),
             rep(x = sd(s2)/2, times = No_Time_To_Repeat),
             rep(x = sd(s3)/2, times = No_Time_To_Repeat))


dfAmAcid["Sample"] <- Sample
dfAmAcid["Std_Dev"] <- Std_Dev


p17 <- ggplot(data=dfAmAcid, 
       aes(x=Codon, y=X.1000, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle(a.a.) +
  theme_minimal() +
  xlab("Codon") +
  ylab("Frequency (Per 1000)") +
  scale_fill_discrete(name = "CU Table", labels = c("Existing (Kazuza)", "Current", "N. tabacum"))+
  geom_errorbar(aes(ymin=X.1000-Std_Dev, ymax=X.1000+Std_Dev), width=.2,
                position=position_dodge(.9))

# Tryptophan - Trp 

a.a. <- "Trp"
No_Time_To_Repeat <- as.numeric(table(dfExisting$AmAcid == a.a.)["TRUE"])

dfAmAcid <- rbind(
  (subset(dfExisting, AmAcid == a.a.))[,c(2,4)],
  (subset(dfCurrent, AmAcid == a.a.))[,c(2,4)],
  (subset(dfN.Tabacum, AmAcid == a.a.))[,c(2,4)]
)

Sample <- c(rep(x = "Existing", times = No_Time_To_Repeat),
            rep(x = "Current", times = No_Time_To_Repeat),
            rep(x = "N. tabacum", times = No_Time_To_Repeat))

str1 <- (subset(dfExisting, AmAcid == a.a.))[,4]
s1 <-  unlist(str1)

str2 <- (subset(dfCurrent, AmAcid == a.a.))[,4]
s2 <-  unlist(str2)

str3 <- (subset(dfN.Tabacum, AmAcid == a.a.))[,4]
s3 <-  unlist(str3)

Std_Dev <- c(rep(x = sd(s1)/2, times = No_Time_To_Repeat),
             rep(x = sd(s2)/2, times = No_Time_To_Repeat),
             rep(x = sd(s3)/2, times = No_Time_To_Repeat))


dfAmAcid["Sample"] <- Sample
dfAmAcid["Std_Dev"] <- Std_Dev


p18<-ggplot(data=dfAmAcid, 
       aes(x=Codon, y=X.1000, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle(a.a.) +
  theme_minimal() +
  xlab("Codon") +
  ylab("Frequency (Per 1000)") +
  scale_fill_discrete(name = "CU Table", labels = c("Existing (Kazuza)", "Current", "N. tabacum"))+
  geom_errorbar(aes(ymin=X.1000-Std_Dev, ymax=X.1000+Std_Dev), width=.2,
                position=position_dodge(.9))

# Tyrosine - Tyr 

a.a. <- "Tyr"
No_Time_To_Repeat <- as.numeric(table(dfExisting$AmAcid == a.a.)["TRUE"])

dfAmAcid <- rbind(
  (subset(dfExisting, AmAcid == a.a.))[,c(2,4)],
  (subset(dfCurrent, AmAcid == a.a.))[,c(2,4)],
  (subset(dfN.Tabacum, AmAcid == a.a.))[,c(2,4)]
)

Sample <- c(rep(x = "Existing", times = No_Time_To_Repeat),
            rep(x = "Current", times = No_Time_To_Repeat),
            rep(x = "N. tabacum", times = No_Time_To_Repeat))

str1 <- (subset(dfExisting, AmAcid == a.a.))[,4]
s1 <-  unlist(str1)

str2 <- (subset(dfCurrent, AmAcid == a.a.))[,4]
s2 <-  unlist(str2)

str3 <- (subset(dfN.Tabacum, AmAcid == a.a.))[,4]
s3 <-  unlist(str3)

Std_Dev <- c(rep(x = sd(s1)/2, times = No_Time_To_Repeat),
             rep(x = sd(s2)/2, times = No_Time_To_Repeat),
             rep(x = sd(s3)/2, times = No_Time_To_Repeat))


dfAmAcid["Sample"] <- Sample
dfAmAcid["Std_Dev"] <- Std_Dev


p19 <- ggplot(data=dfAmAcid, 
       aes(x=Codon, y=X.1000, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle(a.a.) +
  theme_minimal() +
  xlab("Codon") +
  ylab("Frequency (Per 1000)") +
  scale_fill_discrete(name = "CU Table", labels = c("Existing (Kazuza)", "Current", "N. tabacum"))+
  geom_errorbar(aes(ymin=X.1000-Std_Dev, ymax=X.1000+Std_Dev), width=.2,
                position=position_dodge(.9))

# Valine - Val 

a.a. <- "Val"
No_Time_To_Repeat <- as.numeric(table(dfExisting$AmAcid == a.a.)["TRUE"])

dfAmAcid <- rbind(
  (subset(dfExisting, AmAcid == a.a.))[,c(2,4)],
  (subset(dfCurrent, AmAcid == a.a.))[,c(2,4)],
  (subset(dfN.Tabacum, AmAcid == a.a.))[,c(2,4)]
)

Sample <- c(rep(x = "Existing", times = No_Time_To_Repeat),
            rep(x = "Current", times = No_Time_To_Repeat),
            rep(x = "N. tabacum", times = No_Time_To_Repeat))

str1 <- (subset(dfExisting, AmAcid == a.a.))[,4]
s1 <-  unlist(str1)

str2 <- (subset(dfCurrent, AmAcid == a.a.))[,4]
s2 <-  unlist(str2)

str3 <- (subset(dfN.Tabacum, AmAcid == a.a.))[,4]
s3 <-  unlist(str3)

Std_Dev <- c(rep(x = sd(s1)/2, times = No_Time_To_Repeat),
             rep(x = sd(s2)/2, times = No_Time_To_Repeat),
             rep(x = sd(s3)/2, times = No_Time_To_Repeat))


dfAmAcid["Sample"] <- Sample
dfAmAcid["Std_Dev"] <- Std_Dev


p20 <- ggplot(data=dfAmAcid, 
       aes(x=Codon, y=X.1000, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle(a.a.) +
  theme_minimal() +
  xlab("Codon") +
  ylab("Frequency (Per 1000)") +
  scale_fill_discrete(name = "CU Table", labels = c("Existing (Kazuza)", "Current", "N. tabacum"))+
  geom_errorbar(aes(ymin=X.1000-Std_Dev, ymax=X.1000+Std_Dev), width=.2,
                position=position_dodge(.9))

rm(dfAmAcid,Sample, Std_Dev, s1, s2, s3, str1, str2, str3, No_Time_To_Repeat, a.a.)
#### grid ####

# https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html

grid.arrange(p1, p2, p3, p4, nrow = 2)
grid.arrange( p5, p6, p7, p8, nrow = 2)
grid.arrange( p9, p10, p11, p12, nrow = 2)
grid.arrange( p13, p14, p15, p16, nrow = 2)
grid.arrange( p17, p18, p19, p20, nrow = 2)

#