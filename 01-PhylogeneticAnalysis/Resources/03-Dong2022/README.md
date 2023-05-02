# Notes for integrating transcriptomes from Dong et al., 2022
- This is a script to parse and combine transcriptomes from Dong et al., (2022);
Dong, Shanshan, et al. "Phylotranscriptomics of liverworts: revisiting the backbone phylogeny and ancestral gene duplications." *Annals of Botany* 130.7 (2022): 951-964.

- To automatically download from [the original data depository](https://figshare.com/articles/dataset/Transcriptome-based_phylogenomic_studies_of_liverworts/14452662) and combine the files, run "01-CombineDong.R". The output folder for a merged fasta file and IDs will be "Combined" in the current directory.

- The transcriptomes included for analysis is defined by the excel file "List_of_Dong42.xlsx". Please modify it if you want to include/exclude sequences from a particular species.

- Please run these scripts after processing Genome and Onekp data to exclude the sequences from other sources.
