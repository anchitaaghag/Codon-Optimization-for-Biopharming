# Codon Optimization Tool for Biopharming

In recent years, the technique of codon optimization has gained popularity in the biopharmaceutical industry to rapidly create and deploy high-yields of commercially valuable protein products. The tobbaco plant, *Nicotiana benthamiana*, is of particular interest due to the plant's rapid growth rate and ease of genetic manipulation. 

Currently, a range of external software tools exists to back-translate given protein sequences to DNA sequences with a user-specified codon usage table. However, the sole *N. benthamiana* [codon usage table currently existing](https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4100) has not been updated in 15 years and is based on a limited number of genes (100 coding sequences). 

The objective of this project is to create a coding script to translate a given protein sequence into a DNA sequence based on a newly constructed codon use table for the target heterologous host, *N. benthamiana*, within 3 months.

## What Can I Do With This?

The code in this project can be used ...

1) To optimize the codon of input protein sequences by using the codon optimization tool [here](https://github.com/anchitaaghag/Codon-Optimization-for-Biopharming/blob/7642f9d98bfbd7849428a40c6a0da9663a253adf/Codon_Optimization_Tool/Codon_Optimization_Script.R).

2) To translate given protein sequence(s) into corresponding DNA sequence(s) via the the Reverse_Translate() function. For example by loading the [Reverse_Translate.R function file](https://github.com/anchitaaghag/Codon-Optimization-for-Biopharming/blob/7642f9d98bfbd7849428a40c6a0da9663a253adf/Codon_Optimization_Tool/Reverse_Translate_Function.R) using the following line of code:

```
source("Reverse_Translate_Function.R")
```
3) To create a codon usage table from coding sequences in NCBI database. For example by adapting the code [here](https://github.com/anchitaaghag/Codon-Optimization-for-Biopharming/blob/7642f9d98bfbd7849428a40c6a0da9663a253adf/Create_Custom_Codon_Usage/Build_Codon_Usage_Table.R) to your species of interest. 

## Prerequisites/Dependencies

Prior to running the codon optimization script, you will need to install the following R packages and download the associated libraries.

```
install.packages("dplyr") 
install.packages("seqinr") 

library("dplyr") 
library("seqinr")
```

### Input File Parameters

You will need to input two files to run the codon optimization tool. First, you will need a text file with the input codon usage table for the species of interest. If you wish to use an updated codon usage for *N. benthamiana* please use [this file](https://github.com/anchitaaghag/Codon-Optimization-for-Biopharming/blob/63c7b9552fed638b1154c84d438cb1eaad6d5a9e/Codon_Optimization_Tool/Updated_Codon_Usage.txt). Please ensure that there are at least three columns in your codon usage file - "Single_Letter_Abbreviation", "Codon", and "Number" - with the single letter abbreviations for each amino acid, the corresponding codons, and the frequencies or empirical codon counts, respectively. 

Then, you will need a text file with the protein sequences you wish to back-translate into a DNA sequence. Please ensure that your protein sequences file has a header or title followed by the corresponding sequence on the following line. An example of the protein sequence file format can be seen [here](https://github.com/anchitaaghag/Codon-Optimization-for-Biopharming/blob/b9df43a085aed205f94423bedfe418a9dd3a2dda/Codon_Optimization_Tool/Example_Protein_Sequences.txt).

## How to Use the Codon Optimization Tool

A step by step walkthrough that shows you all the fields that require user input & how to use the [codon optimization script](https://github.com/anchitaaghag/Codon-Optimization-for-Biopharming/blob/63c7b9552fed638b1154c84d438cb1eaad6d5a9e/Codon_Optimization_Tool/Codon_Optimization_Script.R)

1) Change the name of the file with your codon usage table. If using the [updated codon usage table](https://github.com/anchitaaghag/Codon-Optimization-for-Biopharming/blob/63c7b9552fed638b1154c84d438cb1eaad6d5a9e/Codon_Optimization_Tool/Updated_Codon_Usage.txt) for *N. benthamiana*, leave this field as follows:

```
Chosen_Codon_Usage <- "Updated_Codon_Usage.txt"
```

2) Change the *Your_File_Name_Here* to the file containing your protein sequence(s).

```
Chosen_Codon_Usage <- "Your_File_Name_Here.txt"
```

3) Change the *Your_Desired_Results_File_Name_Here* to indicate a name for your output results file.

```
Result_File_Name <- "Your_Desired_Results_File_Name_Here.txt"
```

4) Run the file. The results will be outputted as a file to your current working directory.

## Example Usage of the Codon Optimization Tool

This is an example of the [input codon usage](https://github.com/anchitaaghag/Codon-Optimization-for-Biopharming/blob/63c7b9552fed638b1154c84d438cb1eaad6d5a9e/Codon_Optimization_Tool/Updated_Codon_Usage.txt).

<img width="685" alt="Screen Shot 2022-08-07 at 2 28 39 AM" src="https://user-images.githubusercontent.com/92746188/183278301-483b6aa4-817d-4bcd-bb8c-8531a77f83d6.png">

This is an example of the [input protein sequences](https://github.com/anchitaaghag/Codon-Optimization-for-Biopharming/blob/b9df43a085aed205f94423bedfe418a9dd3a2dda/Codon_Optimization_Tool/Example_Protein_Sequences.txt). These were two protein sequences taken from [UniProt](https://www.uniprot.org/uniprotkb?query=nicotiana%20benthamiana).

<img width="684" alt="Screen Shot 2022-08-07 at 2 27 46 AM" src="https://user-images.githubusercontent.com/92746188/183278311-759281a7-8855-4712-84db-ad384727ea9e.png">

Next, the following were the user inputs made to the file prior to running.

```
Chosen_Codon_Usage <- "Updated_Codon_Usage.txt"

Protein_Sequence_File <- "Example_Protein_Sequences.txt"

Result_File_Name <- "Example_Results.txt"
```

Then, the file was run to generate the optimized DNA sequences results. 

The output results were as follows: 

<img width="684" alt="Screen Shot 2022-08-07 at 2 27 00 AM" src="https://user-images.githubusercontent.com/92746188/183278317-8c13363a-abd9-4786-970b-db70b67bfb01.png">


These output files are also available for viewing [here](https://github.com/anchitaaghag/Codon-Optimization-for-Biopharming/blob/b9df43a085aed205f94423bedfe418a9dd3a2dda/Codon_Optimization_Tool/Example_Results.txt).

## Authors

* **Anchitaa Ghag** - [anchitaaghag](https://github.com/anchitaaghag)

## Advisors

* **Dr. Andrew Hamilton-Wright** - [andrewhw](https://github.com/andrewhw) Associate Professor, Department of Computer Science, University of Guelph, Guelph, ON N1G 2W1, Canada. *(Computational/Informatics Expertise)*
* **Dr. Doug Cossar** - VP Research, PlantForm Corporation Canada, Toronto, ON M4S 3E2, Canada. *(External/Industry Advisor)*
* **Dr. Jennifer Geddes-McAlister** - Assistant Professor, Department of Molecular and Cellular Biology, University of Guelph, Guelph, ON N1G 2W1, Canada. *(Biological Expertise)*

## Acknowledgments

* This project is part of my major reasearch project (BINF 6999) to fullfill the requirements of the Master of Bioinformatics degree (University of Guelph, ON, Canada). Project Duration: May 2022 - August 2022. 
* I would like to thank my advisors, Drs. Cossar, Geddes-McAlister, and Hamilton-Wright for their support and advice. See also the "Advisors" section above.
* Dr. Jason McAlister - [jmcalist](https://github.com/jmcalist) for assistance with technical aspects (i.e. viewing GitHub repositories )
* Ben Muselius - for help in interpreting figures and conducting statistical analysis on codon usages
* Geddes-McAlister Lab Members - for their feedback and advice during practice presentations

## Template & License Credit

This file was created by adapting the README.md file template by Billie Thompson [PurpleBooth](https://github.com/PurpleBooth). The template can be found [here](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2).

The following section is based on the [cc-licenses repository](https://github.com/santisoler/cc-licenses) by Santiago Soler [santisoler](https://github.com/santisoler). The following section is derived from the README.md file [here](https://github.com/santisoler/cc-licenses#cc-attribution-sharealike-40-international). In addition, the corresponding full license text file is also availiable via the cc-licenses repository [here](https://github.com/santisoler/cc-licenses/blob/8887424b2a1f1a78fca7efbcc2cd5fd4b1998812/LICENSE-CC-BY-SA).

## License

Shield: [![CC BY-SA 4.0][cc-by-sa-shield]][cc-by-sa]

This work is licensed under a
[Creative Commons Attribution-ShareAlike 4.0 International License][cc-by-sa].

[![CC BY-SA 4.0][cc-by-sa-image]][cc-by-sa]

[cc-by-sa]: http://creativecommons.org/licenses/by-sa/4.0/
[cc-by-sa-image]: https://licensebuttons.net/l/by-sa/4.0/88x31.png
[cc-by-sa-shield]: https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg
