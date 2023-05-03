# Source codes for phylogenetic analysis of KAO (CYP88)
This directory include scripts used for the phylogenetic analysis of KAO (CYP88) homologs.

Working pipeline (All scripts should be run from the 00-Codes folder):

- Step 1: Run "Autotest-1.sh" to generate a initial tree (output file: 04-Candidates-all_FFTNS2_80.treefile.newick)
	1. blastp to search in different local databases. Threshold: e<1e-3. Top 50 genes from each transcriptome/genome were kept for each query sequences.
	2. MAFFT alignment using the FFTNS2 algorithm. Columns with >80% gaps were masked from the alignment.
	3. Phylogenetic inference with IQ-TREE 2 and Ultrafast bootstrap.
- Step 2: Manually select a subtree from the initial tree (easily done with Figtree), and save it as a nexus file "05-Subtree.nex" in the output folder. Chose outgroup sequences and sequences to exclude and modify the corresponding files in 01-Input.
- Step 3: Run "Autotest-2.sh" to test if the sequence selection works okay.
	1. MAFFT alignment using the LINSI or EINSI algorithm. Columns with >98% gaps were masked from the alignment.
	2. Phylogenetic inference with IQ-TREE 2 and Ultrafast bootstrap
- Step 4: Run standard non-parametric analysis with IQ-TREE 2
```
	iqtree -nt AUTO -ntmax 16 -s ./06-Subtree_LINSI.98.phy -m TEST -mset phyml -b 1000 --prefix ./07-Subtree_LINSI_98_np
```
- Step 5: Adjust the tree branch order manually and reroot the tree. Save the adjusted tree as a newick file "07-Subtree_LINSI_80_np.treefile.ordered.newick"
- Step 6: Visualize the tree with "07-Visualization.R". Pay attention to the input files or directory in the script.
- Step 7 (optional): Translate the names of sequence alignments etc. using "08-TranslateAlignment.R". Pay attention to the input files or directory in the script.


The files used in creating the phylogenetic tree of the manuscript is in the directory "ForPublication".
	- KAO-all_seqs.fasta: fasta file for candidate sequences of the final tree
	- KAO-LINSI.fasta: file for sequence alignment, before masking
	- KAO-LINSI-mask80.fasta: file for sequence alignment, masked & used in phylogenetic inference
	- KAO-tree.newick: file for the final phylogenetic tree, in a machine-readable form