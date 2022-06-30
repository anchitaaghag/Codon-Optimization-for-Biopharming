# Processing Kazuza's Codon Counts for each Coding Sequence

A step by step series of commands that tell you how to get codon count data from the Kazuza database are formatted for use in the R script. 

###### _Note: the files used in the R script have already been processed and are ready to be used. This is an additional workflow to potentially incorporate other codon use tables from the database (in the future) or to understand how the aformentioned files were generated._

## Prerequisites

The following commands were run using GNU bash, version 3.2.57(1)-release on macOS Monterey version 12.4.

## Getting Started: Files you will need

You will need your data
This script uses the and

You may also need the file

The [format used by the database](http://www.kazusa.or.jp/codon/current/CODON_LABEL) to label each entry of codon counts is as follows:

```
>LOCUS#CDS\ACCESSION\nt..nt\PID(length)\organism\title\descriptions for the CDS
CGA CGC CGG CGU AGA AGG CUA CUC CUG CUU UUA UUG UCA UCC UCG UCU AGC AGU ACA ACC ACG ACU CCA CCC CCG CCU GCA GCC GCG GCU GGA GGC GGG GGU GUA GUC GUG GUU AAA AAG AAC AAU CAA CAG CAC CAU GAA GAG GAC GAU UAC UAU UGC UGU UUC UUU AUA AUC AUU AUG UGG UAA UAG UGA
```
This format has simply been modified to replace Replace "\" symbols with "---". This is done to ensure that the format of the backslashes does not get mistaken for "\n" newline characters. This can happen for example with a string such as: "1\Nicotiana". The file contains this modified format for simpler copy-and-pasting when in terminal.

## Deatiled Walk-through

First, use the cd command to change your directory to the where the codon count files are located. You can also use the ls command to view all the files present in the directory. For example:

```
Last login: Sat Jun 25 17:36:25 on ttys000
anchitaa@Anchitaas-MacBook-Air ~ % cd Major_Research_Project_2022/06_Code/08_Statistical_Analysis
ls
```

View the contents of the file using the cat command.

```
anchitaa@Anchitaas-MacBook-Air 08_Statistical_Analysis % cat Kazuza_CU_for_each_CDS_Format.txt 
>LOCUS#CDS---ACCESSION---nt---nt---PID(length)---organism---title---descriptions for the CDS
CGA CGC CGG CGU AGA AGG CUA CUC CUG CUU UUA UUG UCA UCC UCG UCU AGC AGU ACA ACC ACG ACU CCA CCC CCG CCU GCA GCC GCG GCU GGA GGC GGG GGU GUA GUC GUG GUU AAA AAG AAC AAU CAA CAG CAC CAU GAA GAG GAC GAU UAC UAU UGC UGU UUC UUU AUA AUC AUU AUG UGG UAA UAG UGA% 
```

Use the sed command to insert this format line as the first line in both the files. This will be eventually be our file's headers.

```
MacBook-Air 08_Statistical_Analysis % sed -i '' "1s/^/>LOCUS#CDS---ACCESSION---nt---nt---PID(length)---organism---title---descriptions for the CDS \nCGA CGC CGG CGU AGA AGG CUA CUC CUG CUU UUA UUG UCA UCC UCG UCU AGC AGU ACA ACC ACG ACU CCA CCC CCG CCU GCA GCC GCG GCU GGA GGC GGG GGU GUA GUC GUG GUU AAA AAG AAC AAU CAA CAG CAC CAU GAA GAG GAC GAU UAC UAU UGC UGU UUC UUU AUA AUC AUU AUG UGG UAA UAG UGA\n/" Kazuza_CU_for_each_CDS_in_N_benthamiana.txt

anchitaa@Anchitaas-MacBook-Air 08_Statistical_Analysis % sed -i '' "1s/^/>LOCUS#CDS---ACCESSION---nt---nt---PID(length)---organism---title---descriptions for the CDS \nCGA CGC CGG CGU AGA AGG CUA CUC CUG CUU UUA UUG UCA UCC UCG UCU AGC AGU ACA ACC ACG ACU CCA CCC CCG CCU GCA GCC GCG GCU GGA GGC GGG GGU GUA GUC GUG GUU AAA AAG AAC AAU CAA CAG CAC CAU GAA GAG GAC GAU UAC UAU UGC UGU UUC UUU AUA AUC AUU AUG UGG UAA UAG UGA\n/" Kazuza_CU_for_each_CDS_in_N_tabacum.txt 
```

Finally, use the awk command to extract every second line from the created file into a new file that contains the codon counts only.

```
anchitaa@Anchitaas-MacBook-Air 08_Statistical_Analysis % awk 'NR % 2 == 0' Kazuza_CU_for_each_CDS_in_N_benthamiana.txt > N_benthamiana_Codon_Counts_Only.txt

anchitaa@Anchitaas-MacBook-Air 08_Statistical_Analysis % awk 'NR % 2 == 0' Kazuza_CU_for_each_CDS_in_N_tabacum.txt > N_tabacum_Codon_Counts_Only.txt
```

## References

* https://stackoverflow.com/questions/9533679/how-to-insert-a-text-at-the-beginning-of-a-file
* https://unix.stackexchange.com/questions/369181/printing-every-nth-line-out-of-a-large-file-into-a-new-file

