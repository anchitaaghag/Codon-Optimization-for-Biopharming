# Codon Optimization for Biopharming: Codon Optimization Tool Testing
# Anchitaa Ghag
 
#### 01 INSTALL PACKAGES & DOWNLOAD LIBRARIES ####

# First, set the working directory by running the following lines.
# setwd()
# getwd()

# This script requires the following package and library.
# If this package has not yet been installed please remove the "#" and run the following lines.

#install.packages("dplyr") 

library("dplyr") 

#### 02 DATA AQUISITION: LOAD DATA FROM FILES ####

# Read in the file with the CAI and GC content values obtained from SolGenomics and CAIcal server.

CAI_Data <- read.table(file = "CAI_Values.txt",
           header = TRUE)

GC_Data <- read.table(file = "GC_Values.txt",
                       header = TRUE)

#### 03 DATA SUMMARY ####

# Check if the data has been loaded correctly.

summary(CAI_Data)
summary(GC_Data)

# There 12 columns in total. 10 runs of the optimized sequences. 1 run for the original sequences. And 1 column for the protein ids.

#### 04 NORMALITY & WICOXON SIGNED RANK TEST FOR CAI VALUES ####

# Check if the distribution is normal.

shapiro.test(CAI_Data$Original_CAI_Values)
shapiro.test(CAI_Data$Run_1_CAI_Values)
shapiro.test(CAI_Data$Run_2_CAI_Values)
shapiro.test(CAI_Data$Run_3_CAI_Values)
shapiro.test(CAI_Data$Run_4_CAI_Values)
shapiro.test(CAI_Data$Run_5_CAI_Values)
shapiro.test(CAI_Data$Run_6_CAI_Values)
shapiro.test(CAI_Data$Run_7_CAI_Values)
shapiro.test(CAI_Data$Run_8_CAI_Values)
shapiro.test(CAI_Data$Run_9_CAI_Values) # Not normal.
shapiro.test(CAI_Data$Run_10_CAI_Values) # Not normal.

# Perform a series of paired t-tests comparing the original CAI with the optimized CAI.

Run_1 <- wilcox.test(x = CAI_Data$Run_1_CAI_Values,
            y = CAI_Data$Original_CAI_Values,
            paired = TRUE,
            conf.level = 0.95,
            alternative = "two.sided")
Run_2 <- wilcox.test(x = CAI_Data$Run_2_CAI_Values,
                     y = CAI_Data$Original_CAI_Values,
                     paired = TRUE,
                     conf.level = 0.95,
                     alternative = "two.sided")
Run_3 <- wilcox.test(x = CAI_Data$Run_3_CAI_Values,
                     y = CAI_Data$Original_CAI_Values,
                     paired = TRUE,
                     conf.level = 0.95,
                     alternative = "two.sided")
Run_4 <- wilcox.test(x = CAI_Data$Run_4_CAI_Values,
                     y = CAI_Data$Original_CAI_Values,
                     paired = TRUE,
                     conf.level = 0.95,
                     alternative = "two.sided")
Run_5 <- wilcox.test(x = CAI_Data$Run_5_CAI_Values,
                     y = CAI_Data$Original_CAI_Values,
                     paired = TRUE,
                     conf.level = 0.95,
                     alternative = "two.sided")
Run_6 <- wilcox.test(x = CAI_Data$Run_6_CAI_Values,
                     y = CAI_Data$Original_CAI_Values,
                     paired = TRUE,
                     conf.level = 0.95,
                     alternative = "two.sided")
Run_7 <- wilcox.test(x = CAI_Data$Run_7_CAI_Values,
                     y = CAI_Data$Original_CAI_Values,
                     paired = TRUE,
                     conf.level = 0.95,
                     alternative = "two.sided")
Run_8 <- wilcox.test(x = CAI_Data$Run_8_CAI_Values,
                     y = CAI_Data$Original_CAI_Values,
                     paired = TRUE,
                     conf.level = 0.95,
                     alternative = "two.sided")
Run_9 <- wilcox.test(x = CAI_Data$Run_9_CAI_Values,
                     y = CAI_Data$Original_CAI_Values,
                     paired = TRUE,
                     conf.level = 0.95,
                     alternative = "two.sided")
Run_10 <- wilcox.test(x = CAI_Data$Run_10_CAI_Values,
                     y = CAI_Data$Original_CAI_Values,
                     paired = TRUE,
                     conf.level = 0.95,
                     alternative = "two.sided")

P_Val <- c(Run_1$p.value,Run_2$p.value,Run_3$p.value,Run_4$p.value,Run_5$p.value,Run_6$p.value,Run_7$p.value,Run_8$p.value,Run_9$p.value,Run_10$p.value)

# Are there significant differences between original and optimized? P ≤ 0.05?
# If true, then can reject the null hypothesis.

P_Val < 0.05
P_Val == 0.05

# None of the replicates had a p-value ≤ 0.05. No significant differences between original and optimized CAI values.

#### 05 NORMALITY & WICOXON SIGNED RANK TEST FOR GC VALUES ####

# Check if the distribution is normal.

shapiro.test(GC_Data$Original_GC._Values)
shapiro.test(GC_Data$Run_1_GC_Values)
shapiro.test(GC_Data$Run_2_GC_Values) # Not normal.
shapiro.test(GC_Data$Run_3_GC_Values) # Not normal.
shapiro.test(GC_Data$Run_4_GC_Values) # Not normal.
shapiro.test(GC_Data$Run_5_GC_Values)
shapiro.test(GC_Data$Run_6_GC_Values) # Not normal.
shapiro.test(GC_Data$Run_7_GC_Values) # Not normal.
shapiro.test(GC_Data$Run_8_GC_Values)
shapiro.test(GC_Data$Run_9_GC_Values) # Not normal.
shapiro.test(GC_Data$Run_10_GC_Values) # Not normal.

# Perform a series of paired t-tests comparing the original CAI with the optimized CAI.

Run_1 <- wilcox.test(x = GC_Data$Run_1_GC_Values,
                     y = GC_Data$Original_GC._Values,
                     paired = TRUE,
                     conf.level = 0.95,
                     alternative = "two.sided")
Run_2 <- wilcox.test(x = GC_Data$Run_2_GC_Values,
                     y = GC_Data$Original_GC._Values,
                     paired = TRUE,
                     conf.level = 0.95,
                     alternative = "two.sided")
Run_3 <- wilcox.test(x = GC_Data$Run_3_GC_Values,
                     y = GC_Data$Original_GC._Values,
                     paired = TRUE,
                     conf.level = 0.95,
                     alternative = "two.sided")
Run_4 <- wilcox.test(x = GC_Data$Run_4_GC_Values,
                     y = GC_Data$Original_GC._Values,
                     paired = TRUE,
                     conf.level = 0.95,
                     alternative = "two.sided")
Run_5 <- wilcox.test(x = GC_Data$Run_5_GC_Values,
                     y = GC_Data$Original_GC._Values,
                     paired = TRUE,
                     conf.level = 0.95,
                     alternative = "two.sided")
Run_6 <- wilcox.test(x = GC_Data$Run_6_GC_Values,
                     y = GC_Data$Original_GC._Values,
                     paired = TRUE,
                     conf.level = 0.95,
                     alternative = "two.sided")
Run_7 <- wilcox.test(x = GC_Data$Run_7_GC_Values,
                     y = GC_Data$Original_GC._Values,
                     paired = TRUE,
                     conf.level = 0.95,
                     alternative = "two.sided")
Run_8 <- wilcox.test(x = GC_Data$Run_8_GC_Values,
                     y = GC_Data$Original_GC._Values,
                     paired = TRUE,
                     conf.level = 0.95,
                     alternative = "two.sided")
Run_9 <- wilcox.test(x = GC_Data$Run_9_GC_Values,
                     y = GC_Data$Original_GC._Values,
                     paired = TRUE,
                     conf.level = 0.95,
                     alternative = "two.sided")
Run_10 <- wilcox.test(x = GC_Data$Run_10_GC_Values,
                      y = GC_Data$Original_GC._Values,
                      paired = TRUE,
                      conf.level = 0.95,
                      alternative = "two.sided")

P_Val_GC <- c(Run_1$p.value,Run_2$p.value,Run_3$p.value,Run_4$p.value,Run_5$p.value,Run_6$p.value,Run_7$p.value,Run_8$p.value,Run_9$p.value,Run_10$p.value)

# Are there significant differences between original and optimized? P ≤ 0.05?
# If true, then can reject the null hypothesis.

P_Val_GC < 0.05
P_Val_GC == 0.05

# One of the replicates had a p-value ≤ 0.05. No significant differences between original and optimized GC values.

