# RNA-seq analysis

Workflow for RNA-seq analysis:

1. Download genome and prepare Salmon index with "Prepare_Annotations.R"
2. Salmon quasi-mapping for the fastq files.
3. Use R scripts in "DESeq2" to perform differential gene expression analysis, GO enrichment and visualization.
