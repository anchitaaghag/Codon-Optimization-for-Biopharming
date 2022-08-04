# Codon Optimization Tool for Biopharming

In recent years, the technique of codon optimization has gained popularity in the biopharmaceutical industry to rapidly create and deploy high-yields of commercially valuable protein products. The tobbaco plant, *Nicotiana benthamiana*, is of particular interest due to the plant's rapid growth rate and ease of genetic manipulation. 

Currently, a range of external software tools exists to back-translate given protein sequences to DNA sequences with a user-specified codon usage table. However, the sole *N. benthamiana* [codon usage table currently existing](https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4100) has not been updated in 15 years and is based on a limited number of genes (100 coding sequences). 

The objective of this project is to create a coding script to translate a given protein sequence into a DNA sequence based on a newly constructed codon use table for the target heterologous host, *N. benthamiana*, within 3 months.

## What Can I Do With This?

The code in this project can be used ...

1) To optimize the codon of input protein sequences by using the codon optimization tool here.

2) To translate given protein sequence(s) into corresponding DNA sequence(s) via the the Reverse_Translate() function. For example by loading the Reverse_Translate.R function file using the following line of code:

```
source("Reverse_Translate_Function.R")
```
3) To create a codon usage table from coding sequencies in NCBI database. For example by adapting the code here to you species of interest. 

## Prerequisites/Dependencies

Prior to running the codon optimization script, you will need to install the following R packages and download the associated libraries.

```
install.packages("dplyr") 
install.packages("seqinr") 

library("dplyr") 
library("seqinr")
```

### Input File Parameters

You will need to input two files to run the codon optimization tool. First, you will need a text file with the input codon usage table for the species of intreset. If you wish to use an updated codon usage for N. benthamiana please use the this file. Please ensure that there are at least three columns in your codon usage file. Then, you will need a text file with the protein sequences you wish to back-translate into a DNA sequence. Please ensure that your protein sequences file has a header or title followed by the corresponding sequence on the following line. An example of the protein sequence file format can be seen here.

## How to Use the Codon Optimization Tool

A step by step series of examples that show you all the field that require user input in the file.

1) Change (if using the updated codon usage table for N. benthamiana, leave this field as follows:

```
Chosen_Codon_Usage <- "Updated_Codon_Usage.txt"
```

2) Change the Your_File_Name_Here to the file containing your 

```
Chosen_Codon_Usage <- "Your_File_Name_Here"

```

3)

```
Result_File_Name <- "Your_Desired_Results_File_Name_Here"
```

4) Run

## Example Usage of the Codon Optimization Tool

This is an example of two protein sequences taken from UniProt and the optimized DNA sequences results.

Here is an example of the input codon usage.

Here is an example of the input protein sequences.

These were the user inputs made to the file prior to running.

```
Chosen_Codon_Usage <- "Updated_Codon_Usage.txt"

Protein_Sequence_File <- "Example_Protein_Sequences.txt"

Result_File_Name <- "Example_Results.txt"
```

And then, the file was run and results were saved in the current working direcotry. 

The output results:



These output files are also available for viewing here.

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
