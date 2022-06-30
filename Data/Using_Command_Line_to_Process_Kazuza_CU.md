# Processing Kazuza's Codon Counts for each Coding Sequence

A step by step series of commands that walk through how codon count data from the Kazuza database were formatted for use in the R script. 

###### _Note: the files used in the R script have already been processed and are ready to be used. This is an additional workflow to potentially incorporate other codon use tables from the database (in the future) or to understand how the aformentioned files were generated._

## Prerequisites

The following commands were run using GNU bash, version 3.2.57(1)-release on macOS Monterey version 12.4.

## Navigating the Kazuza Database

For any given species on [the Kazuza database](https://www.kazusa.or.jp/codon/), a list of codon usages for each coding sequence (CDS) used to build the overall codon usage table is available at the bottom of the record. For example, for the _N. benthamiana_ record:

<img width="575" alt="Screen Shot 2022-06-29 at 9 08 51 PM" src="https://user-images.githubusercontent.com/92746188/176571755-2616124e-3fbc-472f-851f-9d0bc544822a.png">

The resulting screen will display each entry with a ">" symbol and descriptive information. This is followed by another line listing the codon counts in a specific order.

<img width="1434" alt="Screen Shot 2022-06-29 at 9 29 19 PM" src="https://user-images.githubusercontent.com/92746188/176573573-0f477cad-448d-4ede-bbaf-e95b99c52a4e.png">

The [format used by the database](http://www.kazusa.or.jp/codon/current/CODON_LABEL) to label each entry of codon counts is as follows:

```
>LOCUS#CDS\ACCESSION\nt..nt\PID(length)\organism\title\descriptions for the CDS
CGA CGC CGG CGU AGA AGG CUA CUC CUG CUU UUA UUG UCA UCC UCG UCU AGC AGU ACA ACC ACG ACU CCA CCC CCG CCU GCA GCC GCG GCU GGA GGC GGG GGU GUA GUC GUG GUU AAA AAG AAC AAU CAA CAG CAC CAU GAA GAG GAC GAU UAC UAU UGC UGU UUC UUU AUA AUC AUU AUG UGG UAA UAG UGA
```

## Files Needed 

For this R script, the list of codon usages for each CDS for [_N. benthamiana_](http://www.kazusa.or.jp/codon/current/species/4100) and [_N. tabacum_](http://www.kazusa.or.jp/codon/current/species/4097) was copied and pasted in a .txt file format. Both .txt files are availabile in the Data folder.

You may also need the file 

The format in this file has simply been modified from the [original format used by the database](http://www.kazusa.or.jp/codon/current/CODON_LABEL) to replace "\" symbols with "---". This is done to ensure that the format of the backslashes does not get mistaken for "\n" newline characters. This can happen for example with a string such as: "1\Nicotiana". The file contains this modified format for simpler copy-and-pasting when in terminal.

## Detailed Walk-through

First, use the **cd** command to change your directory to the where the codon count files are located. You can also use the **ls** command to view all the files present in the directory. For example:

```
cd Major_Research_Project_2022/06_Code/08_Statistical_Analysis
ls
```

(Optional) View the contents of the file using the **cat** command.

```
cat Kazuza_CU_for_each_CDS_Format.txt 
```

(Optinal) The output of the command above:

```
>LOCUS#CDS---ACCESSION---nt---nt---PID(length)---organism---title---descriptions for the CDS
CGA CGC CGG CGU AGA AGG CUA CUC CUG CUU UUA UUG UCA UCC UCG UCU AGC AGU ACA ACC ACG ACU CCA CCC CCG CCU GCA GCC GCG GCU GGA GGC GGG GGU GUA GUC GUG GUU AAA AAG AAC AAU CAA CAG CAC CAU GAA GAG GAC GAU UAC UAU UGC UGU UUC UUU AUA AUC AUU AUG UGG UAA UAG UGA% 
```

Copy and paste the lines in the following command. Ensure that the "\n" newline character is added to the string. 

Use the **sed** command to insert these lines at the beginning of both the files. This will be eventually be our file's headers. 

```
sed -i '' "1s/^/>LOCUS#CDS---ACCESSION---nt---nt---PID(length)---organism---title---descriptions for the CDS \nCGA CGC CGG CGU AGA AGG CUA CUC CUG CUU UUA UUG UCA UCC UCG UCU AGC AGU ACA ACC ACG ACU CCA CCC CCG CCU GCA GCC GCG GCU GGA GGC GGG GGU GUA GUC GUG GUU AAA AAG AAC AAU CAA CAG CAC CAU GAA GAG GAC GAU UAC UAU UGC UGU UUC UUU AUA AUC AUU AUG UGG UAA UAG UGA\n/" Kazuza_CU_for_each_CDS_in_N_benthamiana.txt

sed -i '' "1s/^/>LOCUS#CDS---ACCESSION---nt---nt---PID(length)---organism---title---descriptions for the CDS \nCGA CGC CGG CGU AGA AGG CUA CUC CUG CUU UUA UUG UCA UCC UCG UCU AGC AGU ACA ACC ACG ACU CCA CCC CCG CCU GCA GCC GCG GCU GGA GGC GGG GGU GUA GUC GUG GUU AAA AAG AAC AAU CAA CAG CAC CAU GAA GAG GAC GAU UAC UAU UGC UGU UUC UUU AUA AUC AUU AUG UGG UAA UAG UGA\n/" Kazuza_CU_for_each_CDS_in_N_tabacum.txt 
```

Finally, use the **awk** command to extract every second line from the created file into a new file that contains the codon counts only.

```
awk 'NR % 2 == 0' Kazuza_CU_for_each_CDS_in_N_benthamiana.txt > N_benthamiana_Codon_Counts_Only.txt

awk 'NR % 2 == 0' Kazuza_CU_for_each_CDS_in_N_tabacum.txt > N_tabacum_Codon_Counts_Only.txt
```

## References

* The sed and awk commands have been adapted from the code found here:
* https://stackoverflow.com/questions/9533679/how-to-insert-a-text-at-the-beginning-of-a-file
* https://unix.stackexchange.com/questions/369181/printing-every-nth-line-out-of-a-large-file-into-a-new-file

