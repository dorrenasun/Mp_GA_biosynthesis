# Notes for integrating protein sequences from OneKP project
- This is a script to parse and combine peptide sequences from [the 1000 Plants (OneKP or 1KP) project](https://sites.google.com/a/ualberta.ca/onekp/)
	References:
	- One Thousand Plant Transcriptomes Initiative. Oct 2019. One thousand plant transcriptomes and the phylogenomics of green plants. Nature 574: 679-685
	-Eric Carpenter, … Gane KS Wong. Oct 2019. Access to RNA-sequencing data from 1173 plant species: The 1000 Plant transcriptomes initiative (1KP). GigaScience 8: giz126

- Run "00-Download1kp-2022.r" to automatically download sequences from the OneKP database. The clades to download and merge are included in "Prefixes.csv", which will also provide prefixes for the parsed sequence names. Downloaded sequence files will be in a folder called "Download" in this directory.

- Run "01-CombineOneKP.R" to merge sequences. The outputs will be in a folder called "Extracted" in this directory.