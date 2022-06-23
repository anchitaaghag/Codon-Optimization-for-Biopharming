# Importing CDS Sequences for Proteins without Gene Information
# 10 June 2022
# Anchitaa Ghag

#### 01 INSTALL PACKAGES & DOWNLOAD LIBRARIES ####

#install.packages("devtools")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#devtools::install_github("mhahsler/rBLAST")
library("rBLAST")

#### 02 DATA AQUISITION ####

# Read in the file with the protein names created in 01_Data_

#read_file("")